#Marine BEACON Project (PhD Chapter One)
#Authors: Lola Riesgo 

#General objetive: evaluacion del poder predictivo de los modelos en el cambio de las
#Prevalencias con RANDOM FOREST

#1. Boostraping para métricas de calibración (LogLoss + Brier Score) y métricas de discriminación
#(AUC y TSS)
#2. Modelos individuales para cada nivel de la prevalencia 
#3. Predicciones 
#4. Realibility diagram 

#R version 4.4.2

#If necessary 
#rm(list=ls(all=TRUE)) 
set.seed(123456789)

#libraries 
#libraries 
library(randomForest)
library(gghalves)
library(caret)
library(dplyr)
library(pROC)
library(rworldmap)
require(rworldxtra)
library(openxlsx)
library(reshape)
library(patchwork)
library(mapdata)
library(dplyr)
library(viridis)
library(purrr)
library(pdp)
library(patchwork)
library(hrbrthemes)
library(RColorBrewer)
#library(maptools) #Not in the CRAN Repository (need to find a susbtitute)
world <- map_data("world")
class(world)

#Set working directory 
setwd("D:/MarineBeacon/CopernicusData/data_final")

#Upload data
data <- readRDS("data_complete_v2.rds")
data <-  as.data.frame(data)
class(data)


# Prepare the fitted data --------------------------------------------------------

data <- data %>%
  dplyr::select(set_id, program, observation_date, latitude, longitude, RMM_number, mnkc_epi, 
                chl_integrated, termocline_intensity, sst_front, year, month)

#Transform RMM column into presence/absence

data$RMM_pres <- ifelse(data$RMM_number > 0, 1, 0)

glimpse(data)

#Filter the public data 

data <- data %>%
  filter(program %in% c("DCF Senne (AZTI)", "DCF Senne (IEO)", "DCF Senne (IRD)", "Moratoire ICCAT 2013-? (IRD)",
                        "OCUP - OB (ORTHONGEL)")) %>%
  droplevels()

data$termocline_intensity <- data$termocline_intensity * -1
glimpse(data)

data_fitted <- data


# Prepare the prediction data ---------------------------------------------

load("D:/MarineBeacon/Rcodes/Exploration_analysis/RData_Predictions/df_monthly_1deg.RData")

df_monthly_1deg <- df_monthly_1deg %>%
  dplyr::select(latitude1, longitude1, mnkc_epi, 
                chl_integrated, termocline_intensity, sst_front, year, month)


data_prediction <- df_monthly_1deg

data_prediction <- data_prediction %>%
  dplyr::rename(
    latitude  = latitude1,
    longitude = longitude1
  )

glimpse(data_prediction)

# BOOTSTRAP FOR THE METRICS -----------------------------------------------


df0 <- data_fitted

predictors <- c(
  "chl_integrated",
  "termocline_intensity",
  "sst_front",
  "mnkc_epi",
  "latitude",
  "longitude",
  "month",
  "year"
)

df <- df0 %>%
  dplyr::select(RMM_pres, dplyr::all_of(predictors)) %>%
  mutate(
    # Coacción segura de predictores a numérico si vinieran como factor/char
    across(all_of(predictors), ~ suppressWarnings(as.numeric(as.character(.))))
  ) %>%
  tidyr::drop_na(RMM_pres, dplyr::all_of(predictors)) %>%
  mutate(
    # Garantiza factor binario con niveles "0","1"
    RMM_pres = factor(as.character(RMM_pres), levels = c("0","1")),
    pres_num = as.integer(as.character(RMM_pres))
  )

message("Tamaño tras limpieza: ", nrow(df))
message("Tabla de la respuesta:\n"); print(table(df$RMM_pres))

idx1 <- which(df$RMM_pres == "1")
idx0 <- which(df$RMM_pres == "0")

n1_test <- max(1, floor(0.30 * length(idx1)))
n0_test <- max(1, floor(0.30 * length(idx0)))

test1    <- sample(idx1, n1_test, replace = FALSE)
test0    <- sample(idx0, n0_test, replace = FALSE)
test_idx <- sort(c(test1, test0))

test_df <- df[test_idx, , drop = FALSE]
pool_df <- df[-test_idx, , drop = FALSE]


idx_pres_pool <- which(pool_df$RMM_pres == "1")
idx_abs_pool  <- which(pool_df$RMM_pres == "0")

prevalences <- c(0.50, 0.20, 0.10, 0.01)
n_pos_fixed_target <- 420       # objetivo de presencias fijas
n_iter              <- 100
ntree               <- 535
mtry                <- 2


res <- vector("list", length(prevalences) * n_iter)
k   <- 0L

with_seed(123456789, {
  for (p in prevalences) {
    n_neg_req <- round(n_pos_fixed_target * (1 - p) / p)
    message(sprintf("\n>>> Prevalence = %.3f → pos=%d / neg=%d",
                    p, n_pos_fixed_target, n_neg_req))
    for (i in seq_len(n_iter)) {
      # Muestreo ESTRATIFICADO desde el pool, con reemplazo para garantizar tamaños
      samp_pres <- sample(idx_pres_pool, n_pos_fixed_target, replace = TRUE)
      samp_abs  <- sample(idx_abs_pool,  n_neg_req,         replace = TRUE)
      
      train_df <- pool_df[c(samp_pres, samp_abs), , drop = FALSE]
      
      # Seguridad extra: garantizamos nuevamente no NAs (no debería haber)
      train_df <- train_df %>% tidyr::drop_na(RMM_pres, dplyr::all_of(predictors))
      
      # Ajuste RF (na.action=na.omit por si algo extraño se cuela)
      mod <- randomForest(
        formula    = RMM_pres ~ chl_integrated + termocline_intensity + sst_front +
          mnkc_epi + latitude + longitude + month + year,
        data       = train_df,
        ntree      = ntree,
        mtry       = mtry,
        replace    = FALSE,        # ya controlamos el muestreo fuera
        importance = TRUE,
        na.action  = na.omit
      )
      
      # Predicción en test fijo
      probs_test <- predict(mod, newdata = test_df, type = "prob")[, 2]
      
      # AUC / TSS (pROC)
      roc_obj <- roc(response = test_df$RMM_pres, predictor = probs_test, quiet = TRUE)
      best    <- coords(roc_obj, x = "best", best.method = "youden",
                        ret = c("threshold", "sensitivity", "specificity"))
      thr  <- as.numeric(best["threshold"])
      sens <- as.numeric(best["sensitivity"])
      spec <- as.numeric(best["specificity"])
      tss  <- sens + spec - 1
      aucv <- as.numeric(auc(roc_obj))
      
      # LogLoss (con clipping para evitar -Inf)
      eps        <- .Machine$double.eps
      probs_clip <- pmin(pmax(probs_test, eps), 1 - eps)
      logloss    <- -mean(test_df$pres_num * log(probs_clip) +
                            (1 - test_df$pres_num) * log(1 - probs_clip))
      
      # Brier
      brier <- mean((probs_test - test_df$pres_num)^2)
      
      # Guardar resultados
      k <- k + 1L
      res[[k]] <- data.frame(
        Prevalence  = p,
        Iteration   = i,
        AUC         = aucv,
        Threshold   = thr,
        Sensitivity = sens,
        Specificity = spec,
        TSS         = tss,
        LogLoss     = logloss,
        BrierScore  = brier,
        stringsAsFactors = FALSE
      )
      if (i %% 10 == 0) {
        message(sprintf("Iteración %d/%d", i, n_iter))
      }
    }
  }
})

df_RF_mobmobular <- bind_rows(res)

summary(df_RF_mobmobular)

save(df_RF_mobmobular, 
     file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/RF_Simulated/df_RF_mobmobular.RData")


# ONE MODEL PER PREVALENCE ------------------------------------------------

data_fitted$RMM_pres <- as.factor(data_fitted$RMM_pres)

prevalences      <- c(0.5,  0.20,
                 0.10,  0.01)
n_pos_fixed <- 420    # siempre 420 presencias
n_iter      <- 100    # iteraciones por ratio
ntree       <- 535
mtry        <- 2
#Cogemos los ceros de la primera iteraccion
data_fitted <- na.omit(data_fitted)
idx_pres <- which(data_fitted$RMM_pres == 1)
idx_abs  <- which(data_fitted$RMM_pres == 0)

# 3) Lista para guardar modelos
rf_models <- list()

for (p in prevalences) {
  # 4.1) calcular cuántos ceros hacen falta
  n_neg_req <- round(n_pos_fixed * (1 - p) / p)
  
  # 4.2) muestreo de índices (una sola vez)
  samp_pres <- sample(idx_pres, n_pos_fixed, replace = TRUE)
  samp_abs  <- sample(idx_abs,  n_neg_req,    replace = TRUE)
  samp_idx  <- c(samp_pres, samp_abs)
  
  # 4.3) entrenar el RF sobre ese subconjunto
  mod <- randomForest(
    RMM_pres ~ chl_integrated + termocline_intensity + sst_front +
      mnkc_epi + latitude + longitude + month + year,
    data       = data_fitted[samp_idx, ],
    ntree      = ntree,
    mtry       = mtry,
    replace    = FALSE,   # ya muestreamos afuera
    importance = TRUE
  )
  
  # 4.4) guardar en la lista con un nombre legible
  #    p.ej. rf_05, rf_30, rf_25, ...
  name <- paste0("rf_", sub("^0\\.", "", sprintf("%0.3f", p)))
  rf_models[[name]] <- mod
}


#Para cada modelo sacamos la prediccion 

for (nm in names(rf_models)) {
  mod   <- rf_models[[nm]]
  # Predecimos probabilidades de la clase "1" (presencia)
  preds <- predict(mod, data_prediction, type = "prob")[, 2]
  # Nombre de la variable de salida
  pred_name <- paste0(nm, "_pred")
  # Asignamos en el entorno global
  assign(pred_name, preds, envir = .GlobalEnv)
}

prevalences <- c(0.5, 0.20, 
            0.10,  0.01)

# Data.frame uniendo todas las predicciones
df_preds <- do.call(rbind, lapply(prevalences, function(p) {
  # nombre de la variable con las preds de ese modelo
  nm <- paste0("rf_", sub("^0\\.", "", sprintf("%0.3f", p)), "_pred")
  preds <- get(nm)  
  data.frame(
    Prevalence = factor(p, levels = prevalences),
    Prediction = preds
  )
}))


df_map_long_RFmobmobular <- df_monthly_1deg %>%
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
  mutate(Prevalence = factor(Prevalence,
                             levels = c("0.5","0.20",
                                        "0.10","0.01")))

# save(df_map_long_RFmobmobular, 
#      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/RF_Simulated/df_map_long_RFmobmobular.Rdata")



# CALIBRATION  -----------------------------------------------

# a) Calibracion con bootstrap

eps       <- 1e-6
prev_real <- 0.018

data_fitted$pres_num <- as.numeric(as.character(data_fitted$RMM_pres))

boot_ctrl <- trainControl(
  method           = "boot",
  number           = 50,
  classProbs       = TRUE,
  summaryFunction  = twoClassSummary,
  savePredictions  = "final"
)

# Entrenar calibradores con bootstrap ---
calib_models <- imap(rf_models, function(rf_mod, nm) {
  
  # Probabilidades crudas del RF en el fit dataset
  prob_raw  <- predict(rf_mod, data_fitted, type = "prob")[,2]
  logit_raw <- qlogis(pmin(pmax(prob_raw, eps), 1 - eps))
  
  # Prevalencia en train y pesos
  prev_train <- mean(data_fitted$pres_num)
  w <- ifelse(
    data_fitted$pres_num == 1,
    prev_real      / prev_train,
    (1 - prev_real) / (1 - prev_train)
  )
  
  # Data frame para caret
  df_cal <-data_fitted %>%
    transmute(
      pres = factor(if_else(pres_num == 1, "yes", "no")),
      logit_raw
    ) 
  
  # Entrenar calibrador GLM con bootstrap
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
  
  prob_raw  <- predict(rf_mod, data_fitted, type = "prob")[,2]
  logit_raw <- qlogis(pmin(pmax(prob_raw, eps), 1 - eps))
  
  cal_mod   <- calib_models[[nm]]
  prob_cal  <- predict(
    cal_mod,
    newdata = data.frame(logit_raw = logit_raw),
    type    = "prob"
  )[,"yes"]
  
  tibble(
    prevalences = factor(rep(p, 2 * nrow(data_fitted)), levels = prevalences),
    pres_num    = rep(data_fitted$pres_num, 2),
    Estado      = rep(c("Sin calibrar", "Calibrado"), each = nrow(data_fitted)),
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

df_bins_calibrationRF_Mobmobular <- df_bins

# save(df_bins_calibrationRF_Mobmobular, 
#      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/RF_Simulated/df_bins_calibrationRF_Mobmobular.RData")


ggplot(df_bins_calibrationRF_Mobmobular, aes(x = mean_prob, y = obs_rate, color = Estado)) +
  geom_line() +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(   # ← AQUÍ está el cambio
    values = c("Sin Calibrar" = "#E69F00", "Calibrado" = "#56B4E9"),
    name = NULL
  ) +
  facet_wrap(~prevalences) +
  labs(
    x = "Predicted probability",
    y = "Observed frequency"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.line       = element_line(color = "black", linewidth = 0.4),
    axis.ticks      = element_line(color = "black", linewidth = 0.3),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid      = element_blank()
  )


# Mapas de calibracion 

df_pred_cal <- imap_dfr(rf_models, function(rf_mod, nm) {
  p <- prevalences[nm]  # devuelve, p.ej., 0.5 para "rf_500"
  
  # a) Probabilidades crudas del RF en el grid de predicción
  prob_raw  <- predict(rf_mod, data_prediction, type = "prob")[,2]
  logit_raw <- qlogis(pmin(pmax(prob_raw, eps), 1 - eps))
  
  # b) Pasar por el calibrador correspondiente
  cal_mod   <- calib_models[[nm]]
  prob_cal  <- predict(
    cal_mod,
    newdata = data.frame(logit_raw = logit_raw),
    type    = "prob"
  )[,"yes"]
  
  # c) Devolver un tibble con coordenadas y ratio
  tibble(
    x            = data_prediction$longitude,
    y            = data_prediction$latitude,
    Probability  = prob_cal,
    prevalences        = factor(p, levels = prevalences)  # mismo facet que usabas
  )
})

df_pred_calRF_mobMobular <- df_pred_cal

# save(df_pred_calRF_mobMobular, 
#       file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/RF_Simulated/df_pred_calRF_mobMobular.RData")


pal <- wes_palette("Zissou1", 100, type = "continuous")

p <- ggplot(df_pred_calRF_mobMobular, aes(x = x, y = y, fill = Probability)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal, limits = c(0,0.75), name  = "Calibrated Probability") + 
  coord_equal() +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  facet_wrap(~ prevalences) +
  # geom_tile(
  #   data = df_presencias,
  #   aes(x = longitude1, y = latitude1),
  #   inherit.aes = FALSE,
  #   width  = 1,
  #   height = 1,
  #   fill   = NA,
  #   color  = "black",
  #   linewidth = 0.35
  # ) +
  geom_map(data=world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = "grey") + 
  coord_fixed(xlim = c(-60, 20), ylim = c(-30, 30)) +
  theme(
    legend.position   = "bottom",
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.grid        = element_blank(),
    strip.background  = element_blank(),                        # ← quita el fondo de la etiqueta
    strip.text        = element_text(face = "bold", size = 12)  # ← estilo de letra para A/B/C/D
  )

p

df_diff <- imap_dfr(rf_models, function(rf_mod, nm) {
  p <- prevalences[nm]
  
  # Probabilidad sin calibrar
  prob_raw  <- predict(rf_mod, data_prediction, type = "prob")[,2]
  logit_raw <- qlogis(pmin(pmax(prob_raw, eps), 1 - eps))
  
  # Probabilidad calibrada
  prob_cal  <- predict(
    calib_models[[nm]],
    newdata = data.frame(logit_raw = logit_raw),
    type    = "prob"
  )[,"yes"]
  
  tibble(
    x          =  data_prediction$longitude,
    y          =  data_prediction$latitude,
    prevalences     = factor(p, levels = prevalences),
    Diferencia =  prob_raw - prob_cal
  )
})

df_diffRFMobmobular <- df_diff


library(paletteer)

min_val <- min(df_diffRFMobmobular$Diferencia, na.rm = TRUE)
max_val <- max(df_diffRFMobmobular$Diferencia, na.rm = TRUE)
lim <- max(abs(min_val), abs(max_val))  # límites simétricos alrededor de 0

lim

save(df_diffRFMobmobular, 
     file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/RF_Simulated/df_diffRFMobmobular.RData")

ggplot(df_diffRFMobmobular, aes(x = x, y = y, fill = Diferencia)) +
  geom_tile(width = 1, height = 1) +  # opcional: fija tamaño si tu grid es de 1x1
  scale_fill_paletteer_c(
    "grDevices::Blue-Red 3",
    name   = "Differences",
    limits = c(-0.68, 0.9)
  ) + labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  geom_map(data=world, map = world, aes(long, lat, map_id = region),
           color = "black", fill = "grey") + 
  coord_fixed(xlim = c(-60, 20), ylim = c(-30, 30)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ prevalences) +
  theme_classic()+
  theme(
    legend.position  = "bottom",
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.3),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid       = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 12)
  )

# Variacion de las metricas una vez calibradas las probabilidades

# VARIACIÓN DE LAS MÉTRICAS  -----------------------------------------------------------


df <- data_fitted
df <- na.omit(df)
df$RMM_pres <- factor(df$RMM_pres, levels = c(0,1))
df$pres_num <- as.numeric(as.character(df$RMM_pres))

idx1 <- which(df$RMM_pres == 1)
idx0 <- which(df$RMM_pres == 0)
test1 <- sample(idx1, floor(0.3 * length(idx1)))
test0 <- sample(idx0, floor(0.3 * length(idx0)))
test_idx <- sort(c(test1, test0))

test_df  <- df[test_idx, ]
pool_df  <- df[-test_idx, ]  # solo de aquí tomamos train en cada prevalencia

# Índices en el pool para muestrear
idx_pres_pool <- which(pool_df$RMM_pres == 1)
idx_abs_pool  <- which(pool_df$RMM_pres == 0)

prevalences <- c(0.5,  0.20, 0.10,  0.01)

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
        formula    = RMM_pres ~ chl_integrated + termocline_intensity + sst_front +
          mnkc_epi + latitude + longitude + month + year,
        data       = train_df,
        ntree      = ntree,
        mtry       = mtry,
        replace    = T,        # ya controlamos el muestreo fuera
        importance = F,
        na.action  = na.omit
      )
      
      # --- EVALUAR SIEMPRE EN TEST FIJO ---
      probs_train <- predict(mod, newdata = train_df, type = "prob")[,2]
      probs_test_raw <- predict(mod, newdata = test_df, type = "prob")[,2]
      
      # Calibrador (Platt Scaling con regresión logística)
      obs_train <- as.numeric(as.character(train_df$RMM_pres))
      df_cal  <- data.frame(prob_raw = probs_train, obs = obs_train)
      cal_mod <- glm(obs ~ prob_raw, family = binomial, data = df_cal)
      
      # aplicar calibrador a las probabilidades del TEST
      probs_cal <- predict(
        cal_mod,
        newdata = data.frame(prob_raw = probs_test_raw),
        type = "response"
      )
      
      roc_obj <- roc(test_df$RMM_pres, probs_cal, quiet = TRUE)
      best    <- coords(roc_obj, x = "best", best.method = "youden",
                        ret = c("threshold","sensitivity","specificity"))
      thr  <- best["threshold"]; sens <- best["sensitivity"]; spec <- best["specificity"]
      tss  <- sens + spec - 1
      aucv <- as.numeric(auc(roc_obj))
      
      eps <- .Machine$double.eps
      probs_clip <- pmin(pmax(probs_cal, eps), 1 - eps)
      logloss <- -mean(test_df$pres_num * log(probs_clip) +
                         (1 - test_df$pres_num) * log(1 - probs_clip))
      
      brier <- mean((probs_cal - test_df$pres_num)^2)
      
      results_aftercal[[length(results_aftercal) + 1]] <- data.frame(
        Prevalence  = p,
        AUC         = aucv,
        Threshold   = thr,
        Sensitivity = sens,
        Specificity = spec,
        TSS         = tss,
        LogLoss     = logloss,
        BrierScore  = brier
      )
    }
  }
})

df_cal_prevalencesRF_Mobmobular <- bind_rows(results_aftercal)

names(df_cal_prevalencesRF_Mobmobular)[names(df_cal_prevalencesRF_Mobmobular) == "sensitivity.1"] <- "TSS"
summary(df_cal_prevalencesRF_Mobmobular)

df_cal_prevalencesRF_Mobmobular$Prevalence <- as.factor(df_cal_prevalencesRF_Mobmobular$Prevalence)

load("D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/RF_Simulated/df_RF_mobmobular.RData")

df_RF_mobmobular <-df_RF_mobmobular  %>% 
  mutate(Type = "Uncalibrated")
df_cal_prevalencesRF_Mobmobular  <- df_cal_prevalencesRF_Mobmobular %>% 
  mutate(Type = "Calibrated")

df_combined_RMM <- bind_rows(df_RF_mobmobular, df_cal_prevalencesRF_Mobmobular) %>%
  pivot_longer(
    cols      = c(AUC, sensitivity, specificity, TSS, LogLoss, BrierScore),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence)))



df_combined_RMM <- df_combined_RMM %>%
  mutate(Type = recode(Type, "Uncalibrated" = "Raw"))

##AUC

df_auc <- df_combined_RMM %>%
  filter(Metric %in% c("AUC"))  

p_AUC <- ggplot(df_auc, aes(x = Prevalence, y = Value, fill = Type)) +
  geom_violin(
    aes(group = interaction(Prevalence, Type)),
    position = position_dodge(width = 0.8),
    width = 0.7, trim = FALSE,
    color = "black", size = 0.4, alpha = 0.65
  ) +
  geom_boxplot(
    aes(group = interaction(Prevalence, Type)),
    position = position_dodge(width = 0.8),
    width = 0.12,
    color = "black", fill = "white",
    size = 0.4, outlier.shape = NA
  ) +
  labs(x = "Prevalence", y = "AUC", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.line       = element_line(color = "black", linewidth = 0.4),
    axis.ticks      = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid      = element_blank()
  )

font_add(family = "Times New Roman", regular = "C:/Windows/Fonts/times.ttf")
p_AUC + theme(text = element_text(family = "Times New Roman"))

##TSS

df_TSS <- df_combined_RMM %>%
  filter(Metric %in% c("TSS"))  

p_TSS <- ggplot(df_TSS, aes(x = Prevalence, y = Value, fill = Type)) +
  geom_violin(
    aes(group = interaction(Prevalence, Type)),
    position = position_dodge(width = 0.8),
    width = 0.7, trim = FALSE,
    color = "black", size = 0.4, alpha = 0.65
  ) +
  geom_boxplot(
    aes(group = interaction(Prevalence, Type)),
    position = position_dodge(width = 0.8),
    width = 0.12,
    color = "black", fill = "white",
    size = 0.4, outlier.shape = NA
  ) +
  labs(x = "Prevalence", y = "TSS", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.line       = element_line(color = "black", linewidth = 0.4),
    axis.ticks      = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid      = element_blank()
  )

font_add(family = "Times New Roman", regular = "C:/Windows/Fonts/times.ttf")
p_TSS + theme(text = element_text(family = "Times New Roman"))

#LOGLOSS

df_LogLoss <- df_combined_RMM %>%
  filter(Metric %in% c("LogLoss"))  

p_LogLoss <- ggplot(df_LogLoss, aes(x = Prevalence, y = log(Value), fill = Type)) +
  geom_violin(
    aes(group = interaction(Prevalence, Type)),
    position = position_dodge(width = 0.8),
    width = 0.7, trim = FALSE,
    color = "black", size = 0.4, alpha = 0.65
  ) +
  geom_boxplot(
    aes(group = interaction(Prevalence, Type)),
    position = position_dodge(width = 0.8),
    width = 0.12,
    color = "black", fill = "white",
    size = 0.4, outlier.shape = NA
  ) +
  labs(x = "Prevalence", y = "log(Log Loss)", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.line       = element_line(color = "black", linewidth = 0.4),
    axis.ticks      = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid      = element_blank()
  )

font_add(family = "Times New Roman", regular = "C:/Windows/Fonts/times.ttf")
p_LogLoss + theme(text = element_text(family = "Times New Roman"))

##Brier Score

df_BrierScore <- df_combined_RMM %>%
  filter(Metric %in% c("BrierScore"))  

p_BrierScore <- ggplot(df_BrierScore, aes(x = Prevalence, y = log(Value), fill = Type)) +
  geom_violin(
    aes(group = interaction(Prevalence, Type)),
    position = position_dodge(width = 0.8),
    width = 0.7, trim = FALSE,
    color = "black", size = 0.4, alpha = 0.65
  ) +
  geom_boxplot(
    aes(group = interaction(Prevalence, Type)),
    position = position_dodge(width = 0.8),
    width = 0.12,
    color = "black", fill = "white",
    size = 0.4, outlier.shape = NA
  ) + 
  labs(x = "Prevalence", y = "log(Brier Score)", fill = NULL) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.line       = element_line(color = "black", linewidth = 0.4),
    axis.ticks      = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid      = element_blank()
  )

font_add(family = "Times New Roman", regular = "C:/Windows/Fonts/times.ttf")
p_BrierScore + theme(text = element_text(family = "Times New Roman"))

common_fill <- scale_fill_manual(
  values = c("Raw" = "#E69F00", "Calibrated" = "#56B4E9"),
  name = NULL
)

final_plot <- (p_AUC | p_TSS) / (p_LogLoss | p_BrierScore) +
  plot_layout(guides = "collect") &
  common_fill &
  theme(legend.position = "bottom") 

final_plot






# PARTIAL DEPENDENCE PLOT -------------------------------------------------

vars <- c("month", "mnkc_epi", "chl_integrated", "termocline_intensity", "sst_front")
grid.res <- 30   # resolución de la rejilla
cls      <- "1"  # prob. de la clase 1 (presencia)

library(ggplot2)
library(pdp)

make_pdp_plot_comparativo <- function(rf_model, calib_model, varname, titulo_modelo) {
  
  # PDP sin calibrar
  pd <- pdp::partial(
    object          = rf_model,
    pred.var        = varname,
    train           = data_fitted,
    grid.resolution = grid.res,
    prob            = TRUE,
    which.class     = cls,
    progress        = "none"
  )
  
  # Calibrar probabilidades
  eps <- 1e-6
  pd$prob_raw   <- pd$yhat
  pd$logit_raw  <- qlogis(pmin(pmax(pd$prob_raw, eps), 1 - eps))
  pd$prob_calib <- predict(calib_model, newdata = pd, type = "prob")[, "yes"]
  
  # Convertir a formato largo para ggplot
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
  
  # Plot con ggplot2
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

# Modelos RF
rf_050 <- rf_models[["rf_500"]]
rf_010 <- rf_models[["rf_010"]]

# Modelos calibradores
modcal_050 <- calib_models[["rf_500"]]
modcal_001 <- calib_models[["rf_010"]]

# PDPs comparativos
pdp_comp_050 <- lapply(vars, \(v) make_pdp_plot_comparativo(rf_050, modcal_050, v, "0.5"))
pdp_comp_001 <- lapply(vars, \(v) make_pdp_plot_comparativo(rf_010, modcal_001, v, "0.01"))

# Fila superior y fila inferior
fila_050 <- pdp_comp_050[[1]] /
  pdp_comp_050[[2]] /
  pdp_comp_050[[3]] /
  pdp_comp_050[[4]] /
  pdp_comp_050[[5]]

fila_001 <- pdp_comp_001[[1]] /
  pdp_comp_001[[2]] /
  pdp_comp_001[[3]] /
  pdp_comp_001[[4]] /
  pdp_comp_001[[5]]

# Crear títulos para las columnas
titulo_columna <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "High", size = 6, fontface = "bold") +
  theme_void()

titulo_columna2 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Extremely low", size = 6, fontface = "bold") +
  theme_void()


p_final <-  (titulo_columna | titulo_columna2) / (fila_050 | fila_001) +
  plot_layout(guides = "collect",heights = c(0.05, 1)) &
  theme(legend.position = "bottom")

p_final








