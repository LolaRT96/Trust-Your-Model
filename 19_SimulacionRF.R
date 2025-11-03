
#Marine BEACON Project (PhD Chapter One)
#Authors: Lola Riesgo 
#General objetive: evaluacion del poder predictivo de los modelos en el cambio de los ratios con datos 
#SIMULADOS 

#1. Generar los datos simulados (mismos que en GAM-INLA scripts)
#Datos prediccion 
#Datos fitted models 
# 
#2. Modelos RF optimización de hiperparametros

#3. Modelos con diferentes PREVALENCIAS (nº presenicas/total)

#4. Evaluacion de los modelos
#- Métricas AUC, Sensibilidad, Especificidad, TSS 
#- Density plots
#- Mapa de prediccion
#- Realibility diagram + Brier Score
#- Calibrado platt scaling 
#- Resto de mapas + diagramas + etc

#R version 4.4.2

#If necessary 

rm(list=ls(all=TRUE)) 

library(sp) # para objetos espaciales
#library(rgeos)
library(geoR)
library(dismo)
library(hSDM)
#library(rgdal)
library(spdep)
library(spData)
library(fields)
library(raster)
#library(maptools)
library(gridExtra)
library(ggplot2)
library(rworldmap)
library(rworldxtra)
library(openxlsx)
library(INLA)
library(GGally)
library(ggplot2)
library(ncdf4)
library(corrplot)
library(rasterVis)
library(rnaturalearth)
library(beepr)
library(tidyverse)
library(gstat)     # para simulación de campos gaussianos
library(randomForest)
library(caret)
library(pROC)
library(ranger)     # Random Forest rápido y multihilo
library(foreach)    # Bucles paralelos "foreach"
library(doParallel) # Backend multiproceso para foreach
library(lobstr)  
library(fmesher)
library(withr)

library(viridis)
library(gghalves)
library(patchwork)

#Set the seed

set.seed(123456789)
print(.Random.seed[1:4])

# CARGAMOS DATOS DE LA SIMULACION -----------------------------------------

load("D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/df_simulacion_fit.RData")

load("D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/datos_simulacion_predict.RData")

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



# RANDOM FOREST APLICACION  ------------------------------------------------

df_simulacion_fit <- df_simulacion_fit %>%
  dplyr::select(x, y, temp, time, bathy, pres)

glimpse(df_simulacion_fit)

# 1) Seleccion de hiperparametros  ----------------------------------------

#Transformamos a factor
df_simulacion_fit$pres <- as.factor(df_simulacion_fit$pres)
df_simulacion_fit <- na.omit(df_simulacion_fit) #RF no acepta valores nulos 
minority_size <- min(table(df_simulacion_fit$pres))

# Entrenamos al modelo
rf_model1 <- randomForest(pres ~ ., 
                          data = df_simulacion_fit, 
                          importance = TRUE, 
                          ntree = 250, 
                          sampsize = c(minority_size, minority_size),
                          mtry        <- floor(sqrt(ncol(df_simulacion_fit) - 1)))


print(rf_model1)

#Control del numero de arboles por errores OOB (mayor ntree mas estabilidad, añadir mas arboles no mejora el rendimiento)
#OOB, cada arbol se entrena con bootstrap (muestras aleatorias con reemplazo)
#los datos NO usados son los OOB, y se usan para hacer un test interno en cada arbol 
#El error OOB es la tasa de error promedio en los predichos OOB de todos los arboles 
#El número de arboles determina cuanto crece el modelo

error_oob <- rf_model1$err.rate[, 1]

# Graficar el error OOB
plot(error_oob, type = "l", 
     main = "Error OOB vs. Número de Árboles",
     xlab = "Número de Árboles", ylab = "Error OOB")
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
plot(trees, errors, type = "b", xlab = "Número de árboles", ylab = "Error OOB")
# Ver el mínimo error OOB y en qué número de árbol ocurre
min_error <- min(errors)
best_n_tree <- trees[which.min(errors)]
cat("El menor error OOB es", min_error, "y ocurre con", best_n_tree, "árboles.\n")

##MEJOR NUMERO DE ARBOLES 200

# 3) Juego de prevalencias (CORRECTO) ----------------------------------------


df <- df_simulacion_fit
df$pres <- factor(df$pres, levels = c(0,1))
df$pres_num <- as.numeric(as.character(df$pres))

idx1 <- which(df$pres == 1)
idx0 <- which(df$pres == 0)
test1 <- sample(idx1, floor(0.3 * length(idx1)))
test0 <- sample(idx0, floor(0.3 * length(idx0)))
test_idx <- sort(c(test1, test0))

test_df  <- df[test_idx, ]
pool_df  <- df[-test_idx, ]  # solo de aquí tomamos train en cada prevalencia

# Índices en el pool para muestrear
idx_pres_pool <- which(pool_df$pres == 1)
idx_abs_pool  <- which(pool_df$pres == 0)

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
      
      # --- EVALUAR SIEMPRE EN TEST FIJO ---
      probs_test <- predict(mod, newdata = test_df, type = "prob")[,2]
      
      roc_obj <- roc(test_df$pres, probs_test, quiet = TRUE)
      best    <- coords(roc_obj, x = "best", best.method = "youden",
                        ret = c("threshold","sensitivity","specificity"))
      thr  <- best["threshold"]; sens <- best["sensitivity"]; spec <- best["specificity"]
      tss  <- sens + spec - 1
      aucv <- as.numeric(auc(roc_obj))
      
      eps <- .Machine$double.eps
      probs_clip <- pmin(pmax(probs_test, eps), 1 - eps)
      logloss <- -mean(test_df$pres_num * log(probs_clip) +
                         (1 - test_df$pres_num) * log(1 - probs_clip))
      
      brier <- mean((probs_test - test_df$pres_num)^2)  
      
      res[[length(res) + 1]] <- data.frame(
        Prevalence  = p,
        AUC         = aucv,
        Threshold   = thr,
        Sensitivity = sens,
        Specificity = spec,
        TSS         = tss,
        LogLoss     = logloss,
        BrierScore  = brier     # <- Aquí también
      )
    }
  }
})

df_results_prevalences_fix <- bind_rows(res)

names(df_results_prevalences_fix)[names(df_results_prevalences_fix) == "sensitivity.1"] <- "TSS"
summary(df_results_prevalences_fix)

df_results_prevalences_fix$Prevalence <- as.factor(df_results_prevalences_fix$Prevalence)
summary(df_results_prevalences_fix)


plot_raincloud <- function(df, metric) {
  df <- df %>%
    # reordenamos Prevalence de mayor a menor
    mutate(Prevalence = fct_reorder(as.factor(Prevalence),
                                    as.numeric(as.character(Prevalence)),
                                    .desc = TRUE))
  
  ggplot(df, aes(x = Prevalence, y = .data[[metric]], fill = Prevalence)) +
    # Medio violín (lado derecho)
    geom_half_violin(
      position = position_nudge(x =  +0.15),
      side     = "r",
      alpha    = 0.6,
      width    = 0.8
    ) +
    # Boxplot centrado
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
p_auc  <- plot_raincloud(df_results_prevalences_fix , "AUC")
p_sens <- plot_raincloud(df_results_prevalences_fix , "sensitivity")
p_spec <- plot_raincloud(df_results_prevalences_fix , "specificity")
p_tss  <- plot_raincloud(df_results_prevalences_fix , "TSS")
p_thr  <- plot_raincloud(df_results_prevalences_fix , "threshold")
p_log <- plot_raincloud(df_results_prevalences_fix, "LogLoss")

windows()
(p_auc | p_sens | p_thr) / (p_spec | p_tss | p_log) 

df_metrics_RFSimulation <- df_results_prevalences_fix

save(df_metrics_RFSimulation, 
     file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/RF_Simulated/df_metrics_RFSimulation.RData")

load("D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/df_results_prevalences_fix.RData")

# 4) Un modelo por ratio  --------------------------------------

#a) Creación de los modelos 
df_simulacion_fit$pres      <- factor(df_simulacion_fit$pres,  levels = c(0, 1))

prevalences <- c(0.5,  0.20, 0.10, 0.01)

n_pos_fixed <- 642      
ntree       <- 200
mtry        <- sqrt(ncol(df_simulacion_fit) - 1)

#Cogemos los ceros de la primera iteraccion

idx_pres <- which(df_simulacion_fit$pres == 1)
idx_abs  <- which(df_simulacion_fit$pres == 0)

# 3) Lista para guardar modelos
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
      replace    = FALSE,   # ya muestreamos afuera
      importance = TRUE
    )
    
    #guardamos los modelos en la lista de los modelos 
    name <- paste0("rf_", sub("^0\\.", "", sprintf("%0.3f", p)))
    rf_models[[name]] <- mod
  }
  
})

#5) Prediccion de los modelos  --------------------------------------------


df_simulacion_predict <- df_simulacion_predict %>%
  dplyr::select(x, y, temp, time, bathy)

with_seed(123456789, {
  
  for (nm in names(rf_models)) {
    mod   <- rf_models[[nm]]
    # Predecimos probabilidades de la clase "1" (presencia)
    preds <- predict(mod, df_simulacion_predict, type = "prob")[, 2]
    # Nombre de la variable de salida
    pred_name <- paste0(nm, "_pred")
    # Asignamos en el entorno global
    assign(pred_name, preds, envir = .GlobalEnv)
  }
  
  prevalences <- c(0.5,  0.20, 0.10, 0.01)
  
  # Data.frame uniendo todas las predicciones
  df_preds <- do.call(rbind, lapply(prevalences, function(p) {
    # nombre de la variable con las preds de ese modelo
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

# Density plot prediccion 

windows();ggplot(df_preds, aes(x = Prediction, 
                               color = prevalences, 
                               fill  = prevalences)) +
  geom_density(alpha = 0.4, size = 0.3) +
  scale_color_viridis(discrete = TRUE, option = "C") +
  scale_fill_viridis(discrete = TRUE, option = "C") +
  facet_wrap(~ prevalences, ncol = 5, scales = "free_y") +  # eje y libre por panel
  labs(
    title = "Distribución de probabilidades predichas por ratio de prevalencia",
    x     = "Probabilidad predicha de presencia",
    y     = "Densidad"
  ) +
  # theme_ipsum(base_size = 12) +
  theme(
    legend.position  = "none",              
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5)
  ) +
  theme_classic()

df_prediccion_RF_densityPlot <- df_preds

save(df_prediccion_RF_densityPlot, 
     file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/df_prediccion_RF_densityPlot.RData")



# Mapa de la prediccion

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
    x     = "Longitud",
    y     = "Latitud",
    title = "Mapas de probabilidad por ratio de prevalencia"
  ) +
  facet_wrap(~ ratio, ncol = 5) +
  theme_classic()

df_prediccion_RF_Simulated <- df_pred_long

save(df_prediccion_RF_Simulated,
      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/RF_Simulated/df_prediccion_RF_Simulated.RData")


# Métricas de calibración  ------------------------------------------------

# a) Realiability diagram 

#Para las predicciones sin calibrar

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
    title = "Reliability Diagram – Probabilidades RAW por ratio",
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


# b) CALIBRACIÓN SIN BOOTSRAP ---------------------------------------------

calib_models <- imap(rf_models, function(rf_mod, nm) {
  # 1) Sacar la probabilidad "raw" sobre df_final
  prob_raw_train <- predict(rf_mod, df_simulacion_fit, type = "prob")[,2]
  
  # 2) Data frame para el glm de calibración
  df_cal  <- data.frame(
    prob_raw   = prob_raw_train,
    pres   = df_simulacion_fit$pres
  )
  
  # 3) Ajustar un GLM binomial (Platt scaling)
  glm(pres ~ prob_raw,
      data   = df_cal,
      family = binomial(link = "logit"))
})


ratios   <- c(0.5, 0.20,0.10, 0.01)
suffixes <- c("500","200","100","010")
names(ratios) <- suffixes

# calib_models: lista de glms de Platt‐Scaling, con mismos nombres que rf_models

df_cal_preds <- imap_dfr(rf_models, function(rf_mod, nm) {
  # 1) extraigo el sufijo quitando "rf_"
  suf <- sub("^rf_", "", nm)
  r   <- ratios[suf]           # ahora sí encuentra 0.5, 0.30, …
  
  # 2) predicción cruda y calibrada
  prob_raw <- predict(rf_mod, df_simulacion_predict, type = "prob")[,2]
  cal_mod  <- calib_models[[nm]]
  prob_cal <- predict(
    cal_mod,
    newdata = data.frame(prob_raw = prob_raw),
    type    = "response"
  )
  
  # 3) devuelvo un tibble con PREVALENCE y PROB_CAL
  tibble(
    ratio = factor(r, levels = as.character(ratios)),
    Probability = prob_cal
  )
})


# Density‐plot faceteado en 5×2, eje y libre
ggplot(df_cal_preds, aes(x = Probability, fill = ratio, color = ratio)) +
  geom_density(alpha = 0.4, size = 0.3) +
  scale_color_viridis(discrete = TRUE, option = "C") +
  scale_fill_viridis(discrete = TRUE, option = "C") +
  facet_wrap(~ ratio, ncol = 5, scales = "free_y") +
  labs(
    title = "Distribución de probabilidades calibradas (Platt Scaling)",
    x     = "Probabilidad calibrada de presencia",
    y     = "Densidad"
  ) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )

ggplot(df_cal_preds, aes(x = Probability, fill = ratio, color = ratio)) +
  geom_density(alpha = 0.4, size = 0.3) +
  scale_color_viridis(discrete = TRUE, option = "C") +
  scale_fill_viridis(discrete = TRUE, option = "C") +
  facet_wrap(~ ratio, ncol = 5, scales = "free_y") +
  coord_cartesian(xlim = c(0, 0.2)) +   # <-- límite en 0–0.5
  labs(
    title = "Distribución de probabilidades calibradas (Platt Scaling)",
    x     = "Probabilidad calibrada de presencia",
    y     = "Densidad"
  ) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )



cal_plot_df_cal <- imap_dfr(rf_models, function(rf_mod, nm) {
  suf <- sub("^rf_", "", nm)
  r   <- as.numeric(suf) / 1000
  
  # Predicciones raw
  prob_raw <- predict(rf_mod, df_simulacion_fit, type = "prob")[,2]
  # Calibradas
  prob_cal <- predict(
    calib_models[[nm]],
    newdata = data.frame(prob_raw = prob_raw),
    type    = "response"
  )
  # Observaciones 0/1
  obs_bin  <- if_else( df_simulacion_fit$pres == 1, 1, 0)
  
  tibble(
    ratio = factor(
      r,
      levels = sort(unique(r)),
      labels = as.character(sort(unique(r)))
    ),
    Obs      = obs_bin,
    Calibrated = prob_cal
  )
})

# 2) Agrupamos en bins y calculamos medias para Calibrated
cal_bins <- cal_plot_df_cal %>%
  mutate(
    bin = cut(Calibrated,
              breaks = seq(0, 1, length.out = 11),
              include.lowest = TRUE)
  ) %>%
  group_by(ratio, bin) %>%
  summarise(
    mean_prob = mean(Calibrated, na.rm = TRUE),
    mean_obs  = mean(Obs,       na.rm = TRUE),
    .groups   = "drop"
  )

# 3) Dibujamos el Reliability Diagram faceteado con eje X libre y Y entre 0 y 1
ggplot(cal_bins, aes(x = mean_prob, y = mean_obs)) +
  geom_line(color = "firebrick", size = 0.7) +
  geom_point(color = "firebrick", size = 1.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray40") +
  facet_wrap(~ ratio, ncol = 5, scales = "free_x") +
  expand_limits(x = 0, y = 0) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Reliability Diagram – Probabilidades CALIBRADAS por ratio",
    x     = "Mean calibrated probability",
    y     = "Observed frequency"
  ) +
  theme_ipsum(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )



# c) CALIBRACIÓN CON BOOTSRAPING --------------------------------------------

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

# Entrenar calibradores con bootstrap ---
calib_models <- imap(rf_models, function(rf_mod, nm) {
  
  # Probabilidades crudas del RF en el fit dataset
  prob_raw  <- predict(rf_mod, df_simulacion_fit, type = "prob")[,2]
  logit_raw <- qlogis(pmin(pmax(prob_raw, eps), 1 - eps))
  
  # Prevalencia en train y pesos
  prev_train <- mean(df_simulacion_fit$pres_num)
  w <- ifelse(
    df_simulacion_fit$pres_num == 1,
    prev_real      / prev_train,
    (1 - prev_real) / (1 - prev_train)
  )
  
  # Data frame para caret
  df_cal <- df_simulacion_fit %>%
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
    Estado      = rep(c("Sin calibrar", "Calibrado"), each = nrow(df_simulacion_fit)),
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

df_bins_calibrationRF <- df_bins

# save(df_bins_calibrationRF, 
#      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/df_bins_calibrationRF.RData")


ggplot(df_bins, aes(x = mean_prob, y = obs_rate, color = Estado)) +
  geom_line() +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~prevalences) +
  labs(
    x = "Predicted probability",
    y = "Observed frequency",
    title = "Reliability diagram antes y después de calibración"
  ) +
  theme_minimal()



df_pred_cal <- imap_dfr(rf_models, function(rf_mod, nm) {
  p <- prevalences[nm]  
  
  # a) Probabilidades crudas del RF en el grid de predicción
  prob_raw  <- predict(rf_mod, df_simulacion_predict, type = "prob")[,2]
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
    x            = df_simulacion_predict$x,
    y            = df_simulacion_predict$y,
    Probability  = prob_cal,
    ratio        = factor(p, levels = prevalences)  # mismo facet que usabas
  )
})

df_pred_calRFSimulated <- df_pred_cal

# save(df_pred_calRFSimulated, 
#      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/RF_Simulated/df_pred_calRFSimulated.RData")

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


df_diff <- imap_dfr(rf_models, function(rf_mod, nm) {
  p <- prevalences[nm]
  
  # Probabilidad sin calibrar
  prob_raw  <- predict(rf_mod, df_simulacion_predict, type = "prob")[,2]
  logit_raw <- qlogis(pmin(pmax(prob_raw, eps), 1 - eps))
  
  # Probabilidad calibrada
  prob_cal  <- predict(
    calib_models[[nm]],
    newdata = data.frame(logit_raw = logit_raw),
    type    = "prob"
  )[,"yes"]
  
  tibble(
    x          = df_simulacion_predict$x,
    y          = df_simulacion_predict$y,
    ratio      = factor(p, levels = prevalences),
    Diferencia =  prob_raw - prob_cal
  )
})

df_diffRFSimulated <- df_diff

save(df_diffRFSimulated, 
     file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/RF_Simulated/df_diffRFSimulated.RData")


# 2) Graficar
windows(); ggplot(df_diff, aes(x = x, y = y, fill = Diferencia)) +
  geom_tile() +
  scale_fill_gradient2(
    name   = "Calibrado - Sin calibrar",
    low    = "blue",
    mid    = "white",
    high   = "red",
    midpoint = 0
  ) +
  coord_equal() +
  labs(
    x = "Longitud",
    y = "Latitud",
    title = "Diferencia de probabilidad (calibrado - sin calibrar)"
  ) +
  facet_wrap(~ ratio, ncol = 5) +
  theme_classic()


# PARTIAL PLOTS -----------------------------------------------------------

#UNCALIBRATED

mod_050 <- rf_models[["rf_500"]]  # p = 0.5
mod_001 <- rf_models[["rf_010"]]  # p = 0.01

vars <- c("temp", "time", "bathy")
grid.res <- 30   # resolución de la rejilla
cls      <- "1"  # prob. de la clase 1 (presencia)


make_pdp_plot <- function(modelo, varname, titulo_modelo) {
  pd <- pdp::partial(
    object          = modelo,
    pred.var        = varname,
    train           = df_simulacion_fit,
    grid.resolution = grid.res,
    prob            = TRUE,         # para clasificación
    which.class     = cls
  )
  autoplot(pd) +
    labs(
      title = sprintf("%s — %s", titulo_modelo, varname),
      x     = varname,
      y     = "Partial effect"
    ) +
    theme_minimal() +
    theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
           legend.position = "top",
           axis.line       = element_line(color = "black", linewidth = 0.4),
           axis.ticks      = element_line(color = "black", linewidth = 0.3),
           panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
           panel.grid      = element_blank(),
           text = element_text(family = "Times New Roman")
    )
}


plots_050 <- lapply(vars, \(v) make_pdp_plot(mod_050, v, "RF p=0.5 (rf_500)"))
plots_001 <- lapply(vars, \(v) make_pdp_plot(mod_001, v, "RF p=0.01 (rf_010)"))


p_final <- (plots_050[[1]] | plots_050[[2]] | plots_050[[3]]) /
  (plots_001[[1]] | plots_001[[2]] | plots_001[[3]])

p_final


#CALIBRATED


# El calibrador no se puede usar directamente para PDP, porque su única función es transformar el logit crudo en una probabilidad calibrada.
# 
# Si quieres ver los PDP calibrados, lo que necesitas hacer es:
# 
# Calibrar luego las probabilidades que resultan del PDP usando el modelo modcal_050 o modcal_001

make_cal_pdp <- function(rf_model, calib_model, varname, titulo_modelo) {
  
  # PDP sin calibrar
  pd <- pdp::partial(
    object          = rf_model,
    pred.var        = varname,
    train           = df_simulacion_fit,
    grid.resolution = grid.res,
    prob            = TRUE,
    which.class     = cls,
    progress        = "none"
  )
  
  # Calibrar la predicción con el modelo de calibración
  eps <- 1e-6
  pd$prob_raw <- pd$yhat
  pd$logit_raw <- qlogis(pmin(pmax(pd$prob_raw, eps), 1 - eps))
  
  # predecir probabilidad calibrada
  pd$prob_cal <- predict(calib_model, newdata = pd, type = "prob")[, "yes"]
  
  # plot
  ggplot(pd, aes_string(x = varname, y = "prob_cal")) +
    geom_line(size = 1, color = "black") +
    labs(
      title = sprintf("Calibrated PDP — %s — %s", titulo_modelo, varname),
      y     = "Partial dependence (calibrated P[pres=1])"
    ) +
    theme_classic()
}

rf_050 <- rf_models[["rf_500"]]
rf_010 <- rf_models[["rf_010"]]

#PDPs calibrados
pdp_cal_050 <- lapply(vars, \(v) make_cal_pdp(rf_050, modcal_050, v, "RF p=0.5"))
pdp_cal_001 <- lapply(vars, \(v) make_cal_pdp(rf_010, modcal_001, v, "RF p=0.01"))

##COMBINACION CALIBRADO Y NO CALIBRADO 

library(ggplot2)
library(pdp)

make_pdp_plot_comparativo <- function(rf_model, calib_model, varname, titulo_modelo) {
  
  # PDP sin calibrar
  pd <- pdp::partial(
    object          = rf_model,
    pred.var        = varname,
    train           = df_simulacion_fit,
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

# Etiquetas visuales para filas
label_row_050 <- plot_spacer() + 
  plot_annotation(title = "Prevalence = 0.5") & 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

label_row_001 <- plot_spacer() + 
  plot_annotation(title = "Prevalence = 0.01") & 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))

# Fila superior y fila inferior
fila_050 <- pdp_comp_050[[1]] /
  pdp_comp_050[[2]] /
  pdp_comp_050[[3]]

fila_001 <- pdp_comp_001[[1]] /
  pdp_comp_001[[2]] /
  pdp_comp_001[[3]]


# Crear títulos para las columnas
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

# VARIACIÓN DE LAS MÉTRICAS  -----------------------------------------------------------


df <- df_simulacion_fit
df$pres <- factor(df$pres, levels = c(0,1))
df$pres_num <- as.numeric(as.character(df$pres))

idx1 <- which(df$pres == 1)
idx0 <- which(df$pres == 0)
test1 <- sample(idx1, floor(0.3 * length(idx1)))
test0 <- sample(idx0, floor(0.3 * length(idx0)))
test_idx <- sort(c(test1, test0))

test_df  <- df[test_idx, ]
pool_df  <- df[-test_idx, ]  # solo de aquí tomamos train en cada prevalencia

# Índices en el pool para muestrear
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
      
      # --- EVALUAR SIEMPRE EN TEST FIJO ---
      probs_train <- predict(mod, newdata = train_df, type = "prob")[,2]
      probs_test_raw <- predict(mod, newdata = test_df, type = "prob")[,2]
      
      # Calibrador (Platt Scaling con regresión logística)
      obs_train <- as.numeric(as.character(train_df$pres))
      df_cal  <- data.frame(prob_raw = probs_train, obs = obs_train)
      cal_mod <- glm(obs ~ prob_raw, family = binomial, data = df_cal)
      
      # aplicar calibrador a las probabilidades del TEST
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

df_cal_prevalencesRF <- bind_rows(results_aftercal)

names(df_cal_prevalencesRF)[names(df_cal_prevalencesRF) == "sensitivity.1"] <- "TSS"
summary(df_cal_prevalencesRF)

df_cal_prevalencesRF$Prevalence <- as.factor(df_cal_prevalencesRF$Prevalence)
summary(df_cal_prevalencesRF)

save(df_cal_prevalencesRF, 
     file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/RF_iterations_results/df_cal_prevalencesRF.RData")


# COMBINACIÓN SIN CALIBRAR + CALIBRADO
df_results_prevalences_fix <-df_results_prevalences_fix  %>% 
  mutate(Type = "Uncalibrated")
df_cal_prevalencesRF  <- df_cal_prevalencesRF %>% 
 mutate(Type = "Calibrated")
 
df_combined <- bind_rows(df_results_prevalences_fix, df_cal_prevalencesRF) %>%
  pivot_longer(
    cols      = c(AUC, sensitivity, specificity, TSS, LogLoss, BrierScore),
    names_to  = "Metric",
    values_to = "Value"
   ) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence)))

df_combinedRF <- df_combined %>%
  filter(Prevalence %in% c("0.5", "0.2", "0.1", "0.01"))  

save(df_combinedRF, 
     file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/df_combinedRF.RData")



















