
#Authors: Maria Dolores Riesgo (IEO-CSIC)
#Research paper: "Trust your model"
#General objective: fit the simulated data with Random Forest models 

#R version 4.4.2

#If necessary 

rm(list=ls(all=TRUE)) 

library(gridExtra)
library(ggplot2)
library(openxlsx)
library(INLA)
library(GGally)
library(corrplot)
library(tidyverse)
library(randomForest)
library(caret)
library(pROC)
library(ranger)    
library(viridis)
library(gghalves)
library(patchwork)

#Set the seed

set.seed(123456789)
print(.Random.seed[1:4])

# load data of the simulation  -----------------------------------------

head(df_simulacion_fit) #df_muestra
head(df_simulacion_predict) #df_ocup

# RANDOM FOREST GRID SEARCH FOR HYPERPARAMETERS  ------------------------------------------------

df_simulacion_fit <- df_simulacion_fit %>%
  dplyr::select(x, y, temp, time, bathy, pres) #ignore the spatial structure

glimpse(df_simulacion_fit)

df_simulacion_fit$pres <- as.factor(df_simulacion_fit$pres)

df_simulacion_fit <- na.omit(df_simulacion_fit) 

minority_size <- min(table(df_simulacion_fit$pres))

# Train model 
rf_model1 <- randomForest(pres ~ ., 
                          data = df_simulacion_fit, 
                          importance = TRUE, 
                          ntree = 250, 
                          sampsize = c(minority_size, minority_size),
                          mtry        <- floor(sqrt(ncol(df_simulacion_fit) - 1)))


print(rf_model1)

error_oob <- rf_model1$err.rate[, 1]

# Plot OBB error
plot(error_oob, type = "l", 
     main = "Error OOB vs. Number of trees",
     xlab = "Number of trees", ylab = "Error OOB")

min_trees <- which.min(error_oob)

abline(v = min_trees, col = "red", lty = 2)

trees <- c(100, 200, 250, 300) 

errors <- sapply(trees, function(n) {
  rf_tmp <- randomForest(pres ~ ., data = df_simulacion_fit, 
                         importance = TRUE, 
                         ntree = 1000, 
                         sampsize = c(minority_size, minority_size),
                         mtry = sqrt(ncol(df_simulacion_fit) - 1))
  rf_tmp$err.rate[n, "OOB"]
})

plot(trees, errors, type = "b", xlab = "Number of trees", ylab = "Error OOB")

min_error <- min(errors)
best_n_tree <- trees[which.min(errors)]

cat("The minor error", min_error, "occurs with", best_n_tree, "trees.\n")

# Assessing model performance with bootstrap resampling -----------------

df <- df_simulacion_fit
df$pres <- factor(df$pres, levels = c(0,1))
df$pres_num <- as.numeric(as.character(df$pres))

idx1 <- which(df$pres == 1)
idx0 <- which(df$pres == 0)
test1 <- sample(idx1, floor(0.3 * length(idx1)))
test0 <- sample(idx0, floor(0.3 * length(idx0)))
test_idx <- sort(c(test1, test0))

test_df  <- df[test_idx, ]
pool_df  <- df[-test_idx, ]  

idx_pres_pool <- which(pool_df$pres == 1)
idx_abs_pool  <- which(pool_df$pres == 0)

#fix parameters 

prevalences <- c(0.5, 0.20,  0.10,  0.01)
n_pos_fixed <- 642
n_iter      <- 100
ntree       <- 200
mtry        <- floor(sqrt(ncol(pool_df) - 1))

res <- list()

with_seed(123456789, {
  for (p in prevalences) {
    n_neg_req <- round(n_pos_fixed * (1 - p) / p)
    message(sprintf(">>> Prevalence = %.3f → pos=%d / neg=%d", p, n_pos_fixed, n_neg_req))
    for (i in seq_len(n_iter)) {
      message(sprintf("    Iter %3d/%3d (prevalence=%.3f)", i, n_iter, p))
      samp_pres <- sample(idx_pres_pool, n_pos_fixed, replace = TRUE)
      samp_abs  <- sample(idx_abs_pool,  n_neg_req,   replace = TRUE)
      train_df  <- pool_df[c(samp_pres, samp_abs), ]
      
      mod <- randomForest(
        pres ~ x + y + time + temp + bathy,
        data       = train_df,
        ntree      = ntree,
        mtry       = mtry,
        replace    = TRUE,   
        importance = FALSE
      )
      probs_test <- predict(mod, newdata = test_df, type = "prob")[,2]
      
      roc_obj <- roc(test_df$pres, probs_test, quiet = TRUE)
      best    <- coords(roc_obj, x = "best", best.method = "youden",
                        ret = c("threshold","sensitivity","specificity"))
      thr  <- best["threshold"] 
      sens <- best["sensitivity"]
      spec <- best["specificity"]
      tss  <- sens + spec - 1
      aucv <- as.numeric(auc(roc_obj))
      brier <- mean((probs_test - test_df$pres_num)^2)  
      
      res[[length(res) + 1]] <- data.frame(
        Prevalence  = p,
        AUC         = aucv,
        Threshold   = thr,
        Sensitivity = sens,
        Specificity = spec,
        TSS         = tss,
        BrierScore  = brier    
      )
    }
  }
})

df_results_prevalences<- bind_rows(res)

names(df_results_prevalences)[names(df_results_prevalences) == "sensitivity.1"] <- "TSS"

df_results_prevalences_fix$Prevalence <- as.factor(df_results_prevalences$Prevalence)

summary(df_results_prevalences)

plot <- function(df, metric) {
  df <- df %>%
    
    mutate(Prevalence = fct_reorder(as.factor(Prevalence),
                                    as.numeric(as.character(Prevalence)),
                                    .desc = TRUE))
  
  ggplot(df,aes(x = Prevalence, y = .data[[metric]], fill = Prevalence)) +
    geom_boxplot(
      width         = 0.1,
      position      = position_dodge(width = 0.3),
      outlier.shape = NA,
      alpha         = 0.8
    ) +
    scale_fill_viridis(discrete = TRUE, option = "C") +
    labs(
      title = metric,
      x     = "Prevalence",
      y     = NULL
    ) +
    theme_classic() +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      strip.background = element_blank()
    )
}

# Llamadas a la función
p_auc  <- plot(df_results_prevalences_fix , "AUC")
p_sens <- plot(df_results_prevalences_fix , "sensitivity")
p_spec <- plot(df_results_prevalences_fix , "specificity")
p_tss  <- plot(df_results_prevalences_fix , "TSS")
p_thr  <- plot(df_results_prevalences_fix , "threshold")
p_log <- plot(df_results_prevalences_fix, "BrierScore")

windows()
(p_auc | p_sens | p_thr) / (p_spec | p_tss | p_log) 


# MODELS AND PREDICTION --------------------------------

df_simulacion_fit$pres      <- factor(df_simulacion_fit$pres,  levels = c(0, 1))

prevalences <- c(0.5,  0.20, 0.10, 0.01)

n_pos_fixed <- 642      
ntree       <- 200
mtry        <- sqrt(ncol(df_simulacion_fit) - 1)

idx_pres <- which(df_simulacion_fit$pres == 1)
idx_abs  <- which(df_simulacion_fit$pres == 0)


rf_models <- list()

with_seed(123456789, { 
  
  for (p in prevalences) {
    # ceros necesarios para cumplir los ratios
    n_neg_req <- round(n_pos_fixed * (1 - p) / p)
    # indices de las filas (ID unicos)
    samp_pres <- sample(idx_pres, n_pos_fixed, replace = TRUE)
    samp_abs  <- sample(idx_abs,  n_neg_req,    replace = TRUE)
    samp_idx  <- c(samp_pres, samp_abs)
    
    #modelos
    mod <- randomForest(
      pres ~ x + y + temp + time + bathy,
      data       = df_simulacion_fit[samp_idx, ],
      ntree      = ntree,
      mtry       = mtry,
      replace    = FALSE,   
      importance = TRUE
    )
    
    name <- paste0("rf_", sub("^0\\.", "", sprintf("%0.3f", p)))
    rf_models[[name]] <- mod
  }
  
})


#Prediction 

df_simulacion_predict <- df_simulacion_predict %>%
  dplyr::select(x, y, temp, time, bathy)

with_seed(123456789, {
  
  for (nm in names(rf_models)) {
    mod   <- rf_models[[nm]]
    
    preds <- predict(mod, df_simulacion_predict, type = "prob")[, 2]
    
    pred_name <- paste0(nm, "_pred")
    
    assign(pred_name, preds, envir = .GlobalEnv)
  }
  
  prevalences <- c(0.5,  0.20, 0.10, 0.01)
  
  
  df_preds <- do.call(rbind, lapply(prevalences, function(p) {
    
    nm <- paste0("rf_", sub("^0\\.", "", sprintf("%0.3f", p)), "_pred")
    preds <- get(nm)  
    data.frame(
      prevalences = factor(p, levels = prevalences),
      Prediction = preds
    )
  }))
  
})

head(df_preds)
library(viridis)


# Plot the maps

df_pred_long <- df_simulacion_predict %>%
  mutate(
    `0.5`   = rf_500_pred,
    `0.20`  = rf_200_pred,
    `0.10`  = rf_100_pred, 
    `0.01`  = rf_010_pred
  ) %>%
  pivot_longer(
    cols = c(`0.5`,`0.20`,
             `0.10`,`0.01`),
    names_to  = "Prevalence",
    values_to = "Probability"
  ) %>%
  mutate(ratio = factor(Prevalence,
                        levels = c("0.5","0.20", 
                                  "0.10","0.01")))

windows();ggplot(df_pred_long, aes(x =x, y = y, fill = Probability)) +
  geom_tile() +
  scale_fill_gradientn(
    name   = "Probabilidad (raw)",
    colours = rainbow(11),
    values = scales::rescale(seq(0, 1, by = 0.1)),
    guide  = "colourbar"
  ) +
  coord_equal()+
  labs(
    x     = "x",
    y     = "y"
  ) +
  facet_wrap(~ ratio, ncol = 5) +
  theme_classic()






# REALIBILITY DIAGRAM  ----------------------------------------------------

#Without calibration 

cal_plot_df <- imap_dfr(rf_models, function(rf_mod, nm) {
  suf <- sub("^rf_", "", nm)
  r   <- as.numeric(suf) / 1000
  
  prob_raw <- predict(rf_mod, df_simulacion_fit, type = "prob")[,2]
  # Obs bien en 0/1
  obs_bin  <- if_else(df_simulacion_fit$pres == 1, 1, 0)
  
  tibble(
    ratio = factor(
      r,
      levels = sort(unique(r)),
      labels = as.character(sort(unique(r)))
    ),
    Obs        = obs_bin,
    Raw        = prob_raw
  )
})

raw_bins <- cal_plot_df %>%
  mutate(
    bin = cut(Raw,
              breaks = seq(0, 1, length.out = 11),
              include.lowest = TRUE)
  ) %>%
  group_by(ratio, bin) %>%
  summarise(
    mean_prob = mean(Raw, na.rm = TRUE),
    mean_obs  = mean(Obs, na.rm = TRUE),
    .groups   = "drop"
  )

ggplot(raw_bins, aes(x = mean_prob, y = mean_obs)) +
  geom_line(color = "steelblue", size = 0.7) +
  geom_point(color = "steelblue", size = 1.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray40") +
  facet_wrap(~ ratio, ncol = 5, scales = "free_x") +
  expand_limits(x = 0, y = 0) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Reliability Diagram – Raw probabilities",
    x     = "Mean predicted probability",
    y     = "Observed frequency"
  ) +
  # theme_ipsum(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5)
  ) +
  theme_classic() 

# CALIBRATION WITH BOOTSTRAP --------------------------------------------

eps       <- 1e-6
prev_real <- 0.018

df_simulacion_fit$pres_num <- as.numeric(as.character(df_simulacion_fit$pres))

boot_ctrl <- trainControl(
  method           = "boot",
  number           = 100,
  classProbs       = TRUE,
  summaryFunction  = twoClassSummary,
  savePredictions  = "final"
)


calib_models <- imap(rf_models, function(rf_mod, nm) {
  
  
  prob_raw  <- predict(rf_mod, df_simulacion_fit, type = "prob")[,2]
  logit_raw <- qlogis(pmin(pmax(prob_raw, eps), 1 - eps))
  
 
  prev_train <- mean(df_simulacion_fit$pres_num)
  w <- ifelse(
    df_simulacion_fit$pres_num == 1,
    prev_real      / prev_train,
    (1 - prev_real) / (1 - prev_train)
  )
  
 
  df_cal <- df_simulacion_fit %>%
    transmute(
      pres = factor(if_else(pres_num == 1, "yes", "no")),
      logit_raw
    )
  
  #  GLM + bootstrap (Platt Scaling)
  train(
    pres ~ logit_raw,
    data      = df_cal,
    method    = "glm",
    family    = binomial(link = "logit"),
    weights   = w,
    metric    = "ROC",
    trControl = boot_ctrl
  )
})

names(prevalences) <- names(rf_models)

df_cal_preds <- imap_dfr(rf_models, function(rf_mod, nm) {
  
  p <- prevalences[nm]  
  
  prob_raw  <- predict(rf_mod, df_simulacion_fit, type = "prob")[,2]
  logit_raw <- qlogis(pmin(pmax(prob_raw, eps), 1 - eps))
  
  cal_mod   <- calib_models[[nm]]
  prob_cal  <- predict(
    cal_mod,
    newdata = data.frame(logit_raw = logit_raw),
    type    = "prob"
  )[,"yes"]
  
  tibble(
    prevalences = factor(rep(p, 2 * nrow(df_simulacion_fit)), levels = prevalences),
    pres_num    = rep(df_simulacion_fit$pres_num, 2),
    Estado      = rep(c("Without calibration", "Calibrated"), each = nrow(df_simulacion_fit)),
    Probability = c(prob_raw, prob_cal)
  )
})

df_bins <- df_cal_preds %>%
  group_by(prevalences, Estado) %>%
  mutate(bin = cut(Probability, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
  group_by(prevalences, Estado, bin) %>%
  summarise(
    mean_prob = mean(Probability),
    obs_rate  = mean(pres_num),
    .groups   = "drop"
  )

#Realibility diagram for calibrated and raw probabilities
ggplot(df_bins, aes(x = mean_prob, y = obs_rate, color = Estado)) +
  geom_line() +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~prevalences) +
  labs(
    x = "Predicted probability",
    y = "Observed frequency"
  ) +
  theme_minimal()


##Calibrated predictions

df_pred_cal <- imap_dfr(rf_models, function(rf_mod, nm) {
  p <- prevalences[nm]  
  
  
  prob_raw  <- predict(rf_mod, df_simulacion_predict, type = "prob")[,2]
  logit_raw <- qlogis(pmin(pmax(prob_raw, eps), 1 - eps))
  
  
  cal_mod   <- calib_models[[nm]]
  prob_cal  <- predict(
    cal_mod,
    newdata = data.frame(logit_raw = logit_raw),
    type    = "prob"
  )[,"yes"]
  
  
  tibble(
    x            = df_simulacion_predict$x,
    y            = df_simulacion_predict$y,
    Probability  = prob_cal,
    ratio        = factor(p, levels = prevalences)  # mismo facet que usabas
  )
})


windows(); ggplot(df_pred_cal, aes(x = x, y = y, fill = Probability)) +
  geom_tile() +
  scale_fill_gradientn(
    name   = "Probabilidad",
    colors = c("#08589e","#2b8cbe","#74a9cf","#bfd3e6",
               "#fdbb84","#e34a33","#a50f15","#5b0100"),
    values = scales::rescale(c(0, 0.15, 0.30, 0.45, 0.55, 0.70, 0.85, 1)),
    guide  = "colourbar"
  ) +
  coord_equal() +
  labs(x = "Longitud", y = "Latitud",
       title = "Mapas de probabilidad CALIBRADA por ratio de prevalencia") +
  facet_wrap(~ ratio, ncol = 5) +
  theme_classic()


# PARTIAL PLOTS -----------------------------------------------------------

#UNCALIBRATED

mod_050 <- rf_models[["rf_500"]]  # p = 0.5
mod_001 <- rf_models[["rf_010"]]  # p = 0.01

vars <- c("temp", "time", "bathy")
grid.res <- 30   
cls      <- "1"  

make_pdp_plot <- function(rf_model, calib_model, varname, titulo_modelo) {
  
  
  pd <- pdp::partial(
    object          = rf_model,
    pred.var        = varname,
    train           = df_simulacion_fit,
    grid.resolution = grid.res,
    prob            = TRUE,
    which.class     = cls,
    progress        = "none"
  )
  
  
  eps <- 1e-6
  pd$prob_raw   <- pd$yhat
  pd$logit_raw  <- qlogis(pmin(pmax(pd$prob_raw, eps), 1 - eps))
  pd$prob_calib <- predict(calib_model, newdata = pd, type = "prob")[, "yes"]
  
  
  pd_long <- pd %>%
    dplyr::select(!!varname, prob_raw, prob_calib) %>%
    tidyr::pivot_longer(
      cols      = c(prob_raw, prob_calib),
      names_to  = "Calibracion",
      values_to = "Dependencia"
    ) %>%
    mutate(
      Calibracion = recode(Calibracion,
                           "prob_raw"   = "Raw",
                           "prob_calib" = "Calibrated")
    )
  
 
  ggplot(pd_long, aes_string(x = varname, y = "Dependencia", color = "Calibracion")) +
    geom_line(size = 1) +
    scale_color_manual(
      values = c("Raw" = "#E69F00", "Calibrated" = "#56B4E9"),
      name= NULL
    ) +
    labs(
      x     = varname,
      y     = "Probability"
    ) +
    theme_classic() +
    theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14),
          legend.position = "top",
          axis.line       = element_line(color = "black", linewidth = 0.4),
          axis.ticks      = element_line(color = "black", linewidth = 0.3),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid      = element_blank(),
          text = element_text(family = "Times New Roman")
    )
}


rf_050 <- rf_models[["rf_500"]]
rf_010 <- rf_models[["rf_010"]]

modcal_050 <- calib_models[["rf_500"]]
modcal_001 <- calib_models[["rf_010"]]


pdp_comp_050 <- lapply(vars, \(v) make_pdp_plot_comparativo(rf_050, modcal_050, v, "0.5"))
pdp_comp_001 <- lapply(vars, \(v) make_pdp_plot_comparativo(rf_010, modcal_001, v, "0.01"))


label_row_050 <- plot_spacer() + 
  plot_annotation(title = "Prevalence = 0.5") & 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

label_row_001 <- plot_spacer() + 
  plot_annotation(title = "Prevalence = 0.01") & 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))


fila_050 <- pdp_comp_050[[1]] /
  pdp_comp_050[[2]] /
  pdp_comp_050[[3]]

fila_001 <- pdp_comp_001[[1]] /
  pdp_comp_001[[2]] /
  pdp_comp_001[[3]]


titulo_columna <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "High", size = 6, fontface = "bold") +
  theme_void()

titulo_columna2 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Extremely low", size = 6, fontface = "bold") +
  theme_void()

p_final <-  (titulo_columna | titulo_columna2) / (fila_050 | fila_001) +
  plot_layout(guides = "collect", heights = c(0.06, 1, 0.06, 1)) &
  theme(legend.position = "bottom")


p_final

# BOOSTRAP WITH CALIBRATION FOR REASSEMENT OF DISCRIMINATION  -----------------------------------------------------------


df <- df_simulacion_fit
df$pres <- factor(df$pres, levels = c(0,1))
df$pres_num <- as.numeric(as.character(df$pres))

idx1 <- which(df$pres == 1)
idx0 <- which(df$pres == 0)
test1 <- sample(idx1, floor(0.3 * length(idx1)))
test0 <- sample(idx0, floor(0.3 * length(idx0)))
test_idx <- sort(c(test1, test0))

test_df  <- df[test_idx, ]
pool_df  <- df[-test_idx, ]  

idx_pres_pool <- which(pool_df$pres == 1)
idx_abs_pool  <- which(pool_df$pres == 0)

prevalences <- c(0.5, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.025, 0.01)
n_pos_fixed <- 642
n_iter      <- 100
ntree       <- 200
mtry        <- floor(sqrt(ncol(pool_df) - 1))

results_aftercal <- list()

with_seed(123456789, {
  for (p in prevalences) {
    n_neg_req <- round(n_pos_fixed * (1 - p) / p)
    message(sprintf(">>> Prevalence = %.3f → pos=%d / neg=%d", p, n_pos_fixed, n_neg_req))
    for (i in seq_len(n_iter)) {
      message(sprintf("    Iter %3d/%3d (prevalence=%.3f)", i, n_iter, p))
      samp_pres <- sample(idx_pres_pool, n_pos_fixed, replace = TRUE)
      samp_abs  <- sample(idx_abs_pool,  n_neg_req,   replace = TRUE)
      train_df  <- pool_df[c(samp_pres, samp_abs), ]
      
      mod <- randomForest(
        pres ~ x + y + time + temp + bathy,
        data       = train_df,
        ntree      = ntree,
        mtry       = mtry,
        replace    = TRUE,   
        importance = FALSE
      )
      
      
      probs_train <- predict(mod, newdata = train_df, type = "prob")[,2]
      probs_test_raw <- predict(mod, newdata = test_df, type = "prob")[,2]
      
      
      obs_train <- as.numeric(as.character(train_df$pres))
      df_cal  <- data.frame(prob_raw = probs_train, obs = obs_train)
      cal_mod <- glm(obs ~ prob_raw, family = binomial, data = df_cal)
      
      
      probs_cal <- predict(
        cal_mod,
        newdata = data.frame(prob_raw = probs_test_raw),
        type = "response"
      )
      
      roc_obj <- roc(test_df$pres, probs_cal, quiet = TRUE)
      best    <- coords(roc_obj, x = "best", best.method = "youden",
                        ret = c("threshold","sensitivity","specificity"))
      thr  <- best["threshold"]; sens <- best["sensitivity"]; spec <- best["specificity"]
      tss  <- sens + spec - 1
      aucv <- as.numeric(auc(roc_obj))
      
      
      results_aftercal[[length(results_aftercal) + 1]] <- data.frame(
        Prevalence  = p,
        AUC         = aucv,
        Threshold   = thr,
        Sensitivity = sens,
        Specificity = spec,
        TSS         = tss,
        
      )
    }
  }
})

df_cal_prevalencesRF <- bind_rows(results_aftercal)

names(df_cal_prevalencesRF)[names(df_cal_prevalencesRF) == "sensitivity.1"] <- "TSS"

df_cal_prevalencesRF$Prevalence <- as.factor(df_cal_prevalencesRF$Prevalence)

summary(df_cal_prevalencesRF)
