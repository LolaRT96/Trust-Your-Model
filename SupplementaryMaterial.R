#Authors: Maria Dolores Riesgo (IEO-CSIC)
#Research paper: "Trust your model"
#General objective: sensitivity analysis of reliability diagrams

#R version 4.4.2

#If necessary 

#R version 4.4.2

#If necessary 
#R version 4.4.2
library(randomForest)
library(dplyr)
library(pROC)

set.seed(123456789)
print(.Random.seed[1:4])

# load data simulation -----------------------------------------

# load("D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/df_simulacion_fit.RData")
# 
# load("D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/datos_simulacion_predict.RData")

head(df_simulacion_fit)
head(df_simulacion_predict)

table(df_simulacion_fit$pres)

mean(df_simulacion_fit$pres)

df_simulacion_fit %>% 
  group_by(time) %>% 
  summarise(prev = mean(pres))

df_simulacion_fit %>%
  group_by(time) %>%
  summarise(
    n_pres = sum(pres == 1),
    n_abs  = sum(pres == 0),
    ratio  = n_pres / n_abs
  )

# RF -------------------------------------------------------------

df <- df_simulacion_fit
df$pres <- factor(df$pres, levels = c(0,1))
df$pres_num <- as.numeric(as.character(df$pres))


idx1 <- which(df$pres == 1)
idx0 <- which(df$pres == 0)
test1 <- sample(idx1, floor(0.3 * length(idx1)))
test0 <- sample(idx0, floor(0.3 * length(idx0)))
test_idx <- sort(c(test1, test0))

test_df <- df[test_idx, ]
pool_df <- df[-test_idx, ]   

# Índices en el pool
idx_pres_pool <- which(pool_df$pres == 1)
idx_abs_pool  <- which(pool_df$pres == 0)


prevalences <- c(0.50, 0.20, 0.10, 0.01)
n_pos_fixed <- 642
n_iter      <- 100
ntree       <- 200
mtry        <- floor(sqrt(ncol(pool_df) - 1))
prev_real   <- 0.01


res_rel     <- list()     # reliability bins


with_seed(123456789, {
  
  for (p in prevalences) {
    
    n_neg_req <- round(n_pos_fixed * (1 - p) / p)
    message(sprintf(">>> Prevalence = %.3f → pos=%d / neg=%d",
                    p, n_pos_fixed, n_neg_req))
    
    for (i in seq_len(n_iter)) {
      
      message(sprintf("    Iter %3d/%3d (prevalence=%.3f)", i, n_iter, p))
      
      samp_pres <- sample(idx_pres_pool, n_pos_fixed, replace = TRUE)
      samp_abs  <- sample(idx_abs_pool,  n_neg_req,   replace = TRUE)
      train_df  <- pool_df[c(samp_pres, samp_abs), ]
      
     
      # Train RF
     
      mod <- randomForest(
        pres ~ x + y + time + temp + bathy,
        data       = train_df,
        ntree      = ntree,
        mtry       = mtry,
        replace    = TRUE,
        importance = FALSE
      )
      
      
      prob_raw <- predict(mod, newdata = test_df, type = "prob")[,2]
      
      # Platt scaling (GLM) using TRAIN data
      
      prob_raw_train <- predict(mod, newdata = train_df, type = "prob")[,2]
      logit_raw_train <- qlogis(pmin(pmax(prob_raw_train, eps), 1 - eps))
      
      prev_train <- mean(train_df$pres_num)
      
      w <- ifelse(
        train_df$pres_num == 1,
        prev_real      / prev_train,
        (1 - prev_real) / (1 - prev_train)
      )
      
      df_cal <- data.frame(
        pres = train_df$pres_num,
        logit_raw = logit_raw_train
      )
      
      mod_cal <- glm(
        pres ~ logit_raw,
        data    = df_cal,
        family  = binomial(link = "logit"),
        weights = w
      )
      
     
      logit_test <- qlogis(pmin(pmax(prob_raw, eps), 1 - eps))
      logit_cal  <- predict(mod_cal, newdata = data.frame(logit_raw = logit_test))
      prob_cal   <- plogis(logit_cal)
      
      
      # reliability bins raw + calibrated
      
      df_rel <- tibble(
        Prevalence = p,
        Iter       = i,
        pres_num   = rep(test_df$pres_num, 2),
        Estado     = rep(c("raw","calibrated"), each = nrow(test_df)),
        Probability = c(prob_raw, prob_cal)
      ) %>%
        mutate(
          bin = cut(Probability, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
        ) %>%
        group_by(Prevalence, Iter, Estado, bin) %>%
        summarise(
          mean_prob = mean(Probability),
          obs_rate  = mean(pres_num),
          .groups   = "drop"
        )
      
      res_rel[[length(res_rel) + 1]] <- df_rel
      
    }
  }
})

df_reliability_bins        <- bind_rows(res_rel)


library(dplyr)
library(ggplot2)

# --- average 

df_mean <- df_reliability_bins %>%
  group_by(Prevalence, bin) %>%
  summarise(
    mean_prob = mean(mean_prob),
    mean_obs  = mean(obs_rate),
    .groups = "drop"
  )

ggplot() +
  geom_line(
    data = df_mean,
    aes(x = mean_prob, y = mean_obs),
    color = "blue", linewidth = 1.2
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ Prevalence) +
  coord_equal() +
  theme_bw() +
  labs(
    x = "Predicted probability",
    y = "Observed frequency",
    title = "Reliability diagram bootstrapped"
  )


df_last <- df_reliability_bins %>%
  filter(Iter == max(Iter), Estado == "raw") %>%
  arrange(Prevalence, mean_prob)  

ggplot(df_last, aes(x = mean_prob, y = obs_rate)) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_point(color = "blue", size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ Prevalence) +
  coord_equal() +
  theme_bw() +
  labs(
    x = "Predicted probability",
    y = "Observed frequency",
    title = "Reliability diagram – Último bootstrap (raw)"
  )

# GAM - INLA    ----------------------------------------

lista_fit <- list()

prevalences <- c(0.5, 0.2, 0.1, 0.01)

n_iter <- 100

for (p in prevalences) {
  message("Procesando prevalencia = ", p)
  
  lista_iter <- vector("list", n_iter)
  
  for (i in 1:n_iter) {
    
    data_pres <- df_simulacion_fit %>% filter(pres == 1)
    data_abs  <- df_simulacion_fit %>% filter(pres == 0)
    n_pres    <- nrow(data_pres)
    n_abs     <- round(n_pres * (1 - p) / p)
    
    sampled_abs <- data_abs %>% slice_sample(n = n_abs)
    
    data_ratio <- bind_rows(data_pres, sampled_abs) %>%
      drop_na(all_of(vars_keep))
    
    
    loc_ratio <- as.matrix(data_ratio[, c("x","y")])
    A.inf     <- inla.spde.make.A(mesh, loc = loc_ratio)
    
   
    stack.fit <- inla.stack(
      data   = list(y = data_ratio$pres),
      A      = list(A.inf, 1),
      effects = list(
        s.index,  
        data_ratio %>% transmute(
          intercept = 1,
          temp, time, bathy
        )
      ),
      tag = "fit"
    )
    
    stack.pred <- inla.stack(
      data    = list(y = rep(NA, n_pred)),
      A       = list(A.pred, 1),
      effects = list(
        s.index,
        df_simulacion_predict %>% transmute(
          intercept = 1,
          temp, time, bathy
        )
      ),
      tag = "pred"
    )
    
    stack.full <- inla.stack(stack.fit, stack.pred)
    
  
    mod <- inla(
      f,
      data              = inla.stack.data(stack.full),
      family            = "binomial",
      control.predictor = list(
        A       = inla.stack.A(stack.full),
        compute = TRUE
      ),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
      verbose           = FALSE
    )
    
    idx_fit <- inla.stack.index(stack.full, tag = "fit")$data
    res_fit <- data.frame(
      prob_mean = plogis(mod$summary.fitted.values[idx_fit, "mean"]),
      pres       = data_ratio$pres,
      ratio      = as.character(p),
      iter       = i
    )
    
    lista_iter[[i]] <- res_fit
  }
  
  # Combinar resultados de las 100 iteraciones
  lista_fit[[as.character(p)]] <- bind_rows(lista_iter)
}


df_fit_all <- bind_rows(lista_fit) %>%
    mutate(prevalence = factor(ratio, 
                                levels = sort(unique(ratio), decreasing = TRUE)))


df_bins_average <- df_fit_all %>%
  group_by(ratio) %>%
  mutate(bin = ntile(prob_mean, 10)) %>%
  group_by(ratio, bin) %>%
  summarise(
    mean_pred = mean(prob_mean, na.rm = TRUE),
    obs_freq  = mean(pres, na.rm = TRUE),
    .groups   = "drop" )

windows();ggplot(df_bins_average, aes(x = mean_pred, y = obs_freq)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_line(color = "black") +
  geom_point(size = 2, shape = 21, fill = "gray70", color = "black") +
  facet_wrap(~ratio, ncol = 5, scale = "free_x") +
  theme_classic() 


df_last <- df_fit_all %>%
  filter(iter == min(iter)) 

df_bins_last <- df_last %>%
  group_by(ratio) %>%
  mutate(bin = ntile(prob_mean, 10)) %>%
  group_by(ratio, bin) %>%
  summarise(
    mean_pred = mean(prob_mean, na.rm = TRUE),
    obs_freq  = mean(pres, na.rm = TRUE),
    .groups   = "drop" )

windows();ggplot(df_bins_last, aes(x = mean_pred, y = obs_freq)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_line(color = "black") +
  geom_point(size = 2, shape = 21, fill = "gray70", color = "black") +
  facet_wrap(~ratio, ncol = 5, scale = "free_x") +
  theme_classic() 

