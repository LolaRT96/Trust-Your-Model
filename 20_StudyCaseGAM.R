#Marine BEACON Project (PhD Chapter One)
#Authors: Lola Riesgo 

#General objetive: evaluacion del poder predictivo de los modelos en el cambio de las
#Prevalencias con GAM

#1. Boostraping para métricas de calibración (LogLoss + Brier Score) y métricas de discriminación
#(AUC y TSS)
#2. Modelos individuales para cada nivel de la prevalencia 
#3. Predicciones 
#4. Realibility diagram 

#R version 4.4.2

#If necessary 
#rm(list=ls(all=TRUE)) 
set.seed(123456789)

# Load necessary libraries
#Libraries
library(sp)
library(tidyr)
library(geoR)
library(INLA)
library(dismo)
library(hSDM)
library(spdep)
library(fields)
library(raster)
library(gridExtra)
library(ggplot2)
library(rworldmap)
require(rworldxtra)
library(openxlsx)
library(reshape)
library(patchwork)
library(viridis)
library(hrbrthemes)
library(INLA)
library(pROC)
library(dplyr)
library(caret)
library(mapdata)
world <- map_data("world")
class(world)

inla.setOption(num.threads = 4)
##READ THE DATA 

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

#Scale the data

vars_to_scale <- c("mnkc_epi", "chl_integrated", "termocline_intensity", "sst_front")

scaling_params_fitted_data <- data  %>% 
  summarise(across(
    all_of(vars_to_scale),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      sd   = ~sd(.x,   na.rm = TRUE)
    )
  ))


# Escalar solo las variables indicadas
data[vars_to_scale] <- scale(data[vars_to_scale])
summary(data)

#Treat the outliers 
data <- data %>%
  filter(if_all(all_of(vars_to_scale), ~ . >= -3.3 & . <= 3.3))
summary(data)

data_fitted <- data



# Prepare the prediction data ---------------------------------------------

load("D:/MarineBeacon/Rcodes/Exploration_analysis/RData_Predictions/df_monthly_1deg.RData")

vars_to_scale <- c("mnkc_epi", "chl_integrated", "termocline_intensity", "sst_front")

# Guardar medias y desviaciones estándar antes de escalar
scaling_params <- df_monthly_1deg %>% 
  summarise(across(
    all_of(vars_to_scale),
    list(
      mean = ~mean(.x, na.rm = TRUE),
      sd   = ~sd(.x,   na.rm = TRUE)
    )
  ))

# Escalar solo las variables indicadas
df_monthly_1deg[vars_to_scale] <- scale(df_monthly_1deg[vars_to_scale])
summary(df_monthly_1deg)

#Treat the outliers 
df_monthly_1deg_scale <- df_monthly_1deg %>%
  filter(if_all(all_of(vars_to_scale), ~ . >= -3.3 & . <= 3.3))
summary(df_monthly_1deg_scale)

data_prediction <- df_monthly_1deg_scale

data_prediction <- data_prediction %>%
  dplyr::rename(
    latitude  = latitude1,
    longitude = longitude1
  )

# BOOTSTRAP FOR THE METRICS -----------------------------------------------

#Mesh

mob_mobular <- data_fitted[, c("latitude", "longitude", "RMM_pres")]
loc <- cbind(mob_mobular$latitude, mob_mobular$longitude)
loc <- as.data.frame(loc)
colnames(loc) <- c("latitude", "longitude")
coordinates(loc) <- ~longitude + latitude
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(loc) <- proj

convhull <- inla.nonconvex.hull(loc) 
mesh <- INLA::inla.mesh.2d(boundary = convhull,
                           max.edge = c(3,6),
                           cutoff = 0.4)

mesh$n 

#Projector matrix

A <- inla.spde.make.A(mesh, loc = loc) 

#SPDE + Spatial Field 
spde <- inla.spde2.pcmatern(mesh, prior.range=c(2, 0.1), prior.sigma=c(1, 0.1))
s.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)


#Priors and Formula 

hyper_pc <- list(prec = list(prior = "pc.prec", param = c(3, 0.05)))

formula_generator <- function(spde) {
  y ~ -1 + intercept +
    f(inla.group(mnkc_epi, n = 20), model = "rw2", hyper = hyper_pc) +
    f(inla.group(sst_front, n = 20), model = "rw2", hyper = hyper_pc) +
    f(inla.group(chl_integrated, n = 20), model = "rw2", hyper = hyper_pc) +
    f(inla.group(termocline_intensity, n = 20), model = "rw2", hyper = hyper_pc) +
    f(month, model = "seasonal", season.length = 12) +
    f(spatial.field, model = spde)
}


#Bootstrapping de prevalencia 

sensitivity_prevalence_inla <- function(data, 
                                        y                 = "RMM_pres",
                                        coords            = c("longitude", "latitude"),
                                        formula_generator,
                                        hyper_pc,
                                        n_iter            = 100,
                                        prevalence        = 0.5,
                                        verbose           = TRUE) {
  
  #Separamos presencias y ausencias fijas
  data_pres <- filter(data, !!sym(y) == 1) #Subcojunto de presencias 
  data_abs  <- filter(data, !!sym(y) == 0) #Subcojunto de ausencias 
  n_pres    <- nrow(data_pres) #total de datos
  # para cada prevalencia deseada necesitamos cambiar los 0, los unos siempre estan fijos
  n_abs_req <- round(n_pres * (1 - prevalence) / prevalence)
  
  #dataframe para guardar las métricas 
  results <- data.frame(
    Iteration   = integer(),
    AUC         = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    TSS         = numeric(),
    LogLoss = numeric(),
    BrierScore = numeric()
  )
  
  for (i in seq_len(n_iter)) { #itera desde 1 =i hasta el final de la secuencia de las iteracciones en este caso 100 
    t0 <- Sys.time() #guardamos por que queremos ver cuanto tarda en ahcer cada iteraccion 
    if (verbose) cat("Prev =", prevalence, "| Iter", i, "de", n_iter, "...\n") #quiero saber por que número de la iteraccion va
    
    #Fijamos siempre la semilla como el numero de la iteracion para reproductibilidad 
    set.seed(i)
    #Aqui le decimos, de todas las observaciones de ausencias (data_abs) toma al azar n_abs_requeridas para alcanzar la prvalencia deseada
    
    sampled_abs <- slice_sample( #slice_sample es mejor sample()
      data_abs, 
      n       = n_abs_req,
      replace = (n_abs_req > nrow(data_abs)) #SI el numero de ausencias es mayor que el que hay en data_abs (NO VA A PASAR)
      #pero para que sea reproducible para otra gente, pues entonces el muestreo es con reemplazo para conseguir llegar el numero de ceros deseado
    )
    
    #Juntamos los ceros y los unos por que ese va a ser el data 
    dat_i <- bind_rows(data_pres, sampled_abs) %>% #presencias y ausencias 
      drop_na(
        !!sym(y), all_of(coords), #viene de lo que se defina en data original 
        mnkc_epi, sst_front, chl_integrated, termocline_intensity, month #viene de lo que se defina en data original
      )
    
    # Pasamos a todos los argumentos que necesitamos para inla
    #MESH Y SPDE + INDICE 
    coords_mat <- as.matrix(dat_i[, coords]) #matriz para construir el mesh 
    hull       <- inla.nonconvex.hull(coords_mat) 
    mesh_i     <- inla.mesh.2d(boundary = hull, max.edge = c(3, 6), 
                               cutoff = 0.4) #conservamos los argumentos
    spde_i     <- inla.spde2.pcmatern(mesh_i,
                                      prior.range = c(2,0.1),
                                      prior.sigma = c(1,0.1))
    A_i        <- inla.spde.make.A(mesh_i, loc = coords_mat)
    s.index    <- inla.spde.make.index("spatial.field", spde_i$n.spde)
    
    stack_i <- inla.stack(
      data    = list(y = dat_i[[y]]), #y=variable respuesta RMM_PRES
      A       = list(A_i, 1),
      effects = list(
        s.index,
        data.frame(
          intercept            = 1,
          mnkc_epi             = dat_i$mnkc_epi,
          sst_front            = dat_i$sst_front,
          chl_integrated       = dat_i$chl_integrated,
          termocline_intensity = dat_i$termocline_intensity,
          month                = dat_i$month
        )
      ),
      tag = "model" #fitted
    )
    
    formula_i <- formula_generator(spde_i)
    
    mod_i     <- inla(
      formula            = formula_i,
      data               = inla.stack.data(stack_i),
      family             = "binomial",
      control.predictor  = list(A = inla.stack.A(stack_i), compute = TRUE),
      control.compute    = list(dic = FALSE, waic = FALSE),
      verbose            = FALSE
    )
    
    idx_obs <- inla.stack.index(stack_i, tag = "model")$data #las observaciones
    fitted  <- mod_i$summary.fitted.values$mean[idx_obs] #las predicciones (probabilidades) sobre las observaciones
    true    <- dat_i[[y]] #Como INLA calcula tambien predicciones en los nodos que no hay datos directos pues necesitamos 
    #los datos que corresponden directamente a los puntos muestreados
    true    <- dat_i[[y]] #dat_i contiene para esa iteracción i todas las filas seleccionadas 
    
    roc_obj <- roc(response = true, predictor = fitted, quiet = TRUE)
    best    <- coords(roc_obj, x = "best", best.method = "youden",
                      ret = c("threshold","sensitivity","specificity"))
    thr  <- best["threshold"]; sens <- best["sensitivity"]; spec <- best["specificity"]
    tss  <- sens + spec - 1
    aucv <- as.numeric(auc(roc_obj))
    
    eps   <- .Machine$double.eps
    pclip <- pmin(pmax(fitted, eps), 1 - eps)
    logloss <- -mean(true * log(pclip) + (1 - true) * log(1 - pclip))
    
    brier <- mean((fitted - true)^2)
    
    results[i, c("AUC","Sensitivity","Specificity","TSS","Cutoff","LogLoss", "BrierScore", "Prevalence")] <-
      c(aucv, sens, spec, tss, thr, logloss, brier, prevalence)#almacenamos en results
    
    #El tiempo que tarda cada iteracción en cerrar el bucle desde t0 inicio 
    if (verbose) {
      dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
      cat("Duración iter", i, ":", dt, "minutos\n\n")
    }
  }
  
  # añadir columna de prevalencia
  results$Prevalence <- prevalence #añadimos en results la prevalencia que hemos usado para no perder track
  return(results) #que me aparezca el frame completo (se puede no poner)
}


prevalence <- c(0.50,  0.20,
                0.10,  0.01)

all_res <- lapply(prevalence, function(p) {
  sensitivity_prevalence_inla(
    data               = data_fitted,
    y                  = "RMM_pres",
    coords             = c("longitude","latitude"),
    formula_generator  = formula_generator,
    hyper_pc           = hyper_pc,
    n_iter             = 100,
    prevalence         = p,
    verbose            = TRUE
  )
})

mobula_GAM <- bind_rows(all_res)

save(mobula_GAM, 
     file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/GAM_Simulated/mobula_GAM.RData")

summary(mobula_GAM)

mobula_GAM$Prevalence <- as.factor(mobula_GAM$Prevalence)

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
p_auc  <- plot_raincloud(mobula_GAM , "AUC")
p_tss  <- plot_raincloud(mobula_GAM, "TSS")
p_brier <- plot_raincloud(mobula_GAM, "BrierScore")
p_log <- plot_raincloud(mobula_GAM, "LogLoss")

windows()
(p_auc | p_tss) / (p_brier | p_log) 


# Un modelo por Prevalencia  ----------------------------------------------
prevalences <- c(0.50,  0.20,
                0.10,  0.01)

vars_keep <- c("RMM_pres","longitude","latitude",
               "month","mnkc_epi","sst_front", "chl_integrated", "termocline_intensity")

loc_pred <- as.matrix(data_prediction[, c("longitude","latitude")])
A.pred   <- inla.spde.make.A(mesh, loc = loc_pred)
n_pred   <- nrow(loc_pred)

f <- y ~ -1 + intercept +
  f(inla.group(mnkc_epi, n = 20), model = "rw2", hyper = hyper_pc) +
  f(inla.group(sst_front, n = 20), model = "rw2", hyper = hyper_pc) +
  f(inla.group(chl_integrated, n = 20), model = "rw2", hyper = hyper_pc) +
  f(inla.group(termocline_intensity, n = 20), model = "rw2", hyper = hyper_pc) +
  f(month, model = "seasonal", season.length = 12) +
  f(spatial.field, model = spde)

data_prediction <- data_prediction %>% mutate(month = as.integer(month))

# Inicializar listas
# lista_probs         <- list()
# lista_fit           <- list()
# lista_spatial_field <- list()
lista_models <- list()

set.seed(123456789) # reproducibilidad del muestreo de ausencias

for (p in prevalences) {
  message("Procesando prevalencia = ", p)
  
  # 1. Submuestra
  data_pres <- data_fitted %>% filter(RMM_pres == 1)
  data_abs  <- data_fitted %>% filter(RMM_pres == 0)
  n_pres    <- nrow(data_pres)
  n_abs     <- round(n_pres * (1 - p) / p)
  sampled_abs <- data_abs %>% slice_sample(n = n_abs)
  
  data_ratio <- bind_rows(data_pres, sampled_abs) %>%
    drop_na(all_of(vars_keep)) %>%
    mutate(month = as.integer(month))
  
  # 2. Matriz A de inferencia (coincide con data_ratio)
  loc_ratio <- as.matrix(data_ratio[, c("longitude","latitude")])
  A.inf     <- inla.spde.make.A(mesh, loc = loc_ratio)
  
  # Comprobaciones rápidas
  stopifnot(nrow(A.inf)  == nrow(data_ratio))
  stopifnot(nrow(A.pred) == n_pred)
  
  # 3. Stacks
  #    OJO: el orden de A y de effects debe coincidir
  stack.fit <- inla.stack(
    data    = list(y = data_ratio$RMM_pres),
    A       = list(A.inf, 1),
    effects = list(
      s.index,
      data_ratio %>% transmute(
        intercept = 1,
        mnkc_epi, sst_front, chl_integrated, termocline_intensity, month
      )
    ),
    tag = "fit"
  )
  
  stack.pred <- inla.stack(
    data    = list(y = rep(NA, n_pred)),    # <- tantas filas como predicciones
    A       = list(A.pred, 1),               # <- A.pred hecho con data_prediction
    effects = list(
      s.index,
      data_prediction %>% transmute(
        intercept = 1,
        mnkc_epi, sst_front, chl_integrated, termocline_intensity, month
      )
    ),
    tag = "pred"
  )
  
  stack.full <- inla.stack(stack.fit, stack.pred)
  
  # 4. Ajuste del modelo
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

  # 5. Predicciones sobre la malla de predicción (stack.pred)
  idx_pred <- inla.stack.index(stack.full, tag = "pred")$data
  res_pred <- data.frame(
    x        = data_prediction$longitude,
    y        = data_prediction$latitude,
    link_mean= mod$summary.fitted.values[idx_pred, "mean"],
    link_sd  = mod$summary.fitted.values[idx_pred, "sd"],
    lower_q  = mod$summary.fitted.values[idx_pred, "0.025quant"],
    upper_q  = mod$summary.fitted.values[idx_pred, "0.975quant"]
  ) %>%
    mutate(
      prob_mean   = plogis(link_mean),
      prob_lower  = plogis(lower_q),
      prob_upper  = plogis(upper_q),
      prob_sd     = link_sd * prob_mean * (1 - prob_mean),
      odds        = prob_mean / (1 - prob_mean),
      ratio_glob  = nrow(data_pres) / nrow(sampled_abs), # opcional: usa p_eff si prefieres
      favorability= odds / (ratio_glob + odds),
      prevalence  = as.character(p)
    )

  # 6. Predicciones en datos de entrenamiento (stack.fit)
  idx_fit <- inla.stack.index(stack.full, tag = "fit")$data
  res_fit <- data.frame(
    prob_mean  = plogis(mod$summary.fitted.values[idx_fit, "mean"]),
    pres       = data_ratio$RMM_pres,      # <- nombre correcto
    prevalence = as.character(p)
  )

  # 7. Campo espacial en nodos de malla
  proj <- inla.mesh.projector(
    mesh,
    xlim = range(mesh$loc[,1]),
    ylim = range(mesh$loc[,2]),
    dims = c(300, 300)
  )
  field_mean <- inla.mesh.project(proj, field = mod$summary.random$spatial.field$mean)

  df_field <- reshape2::melt(field_mean)
  names(df_field) <- c("x_id", "y_id", "mean")
  df_field$x <- proj$x[df_field$x_id]
  df_field$y <- proj$y[df_field$y_id]
  df_field$prevalence <- as.character(p)

  # 8. Guardar
  lista_probs[[as.character(p)]]         <- res_pred
  lista_fit[[as.character(p)]]           <- res_fit
  lista_spatial_field[[as.character(p)]] <- df_field
  lista_models[[as.character(p)]]        <- mod 
}

# #Combinar todos los resultados
# df_prob_prevalence_GAM_mobular <- bind_rows(lista_probs) %>%
#   mutate(prevalence = factor(prevalence, 
#                         levels = sort(unique(prevalence),
#                                                   decreasing = TRUE)))
# 
# save(df_prob_prevalence_GAM_mobular, 
#      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/GAM_Simulated/df_prob_GAM_mobular.RData")
# 
# 
# df_fit_all_GAM_mobular <- bind_rows(lista_fit) %>%
#   mutate(prevalence = factor(prevalence, 
#                              levels = sort(unique(prevalence),
#                                            decreasing = TRUE)))
# 
# save(df_fit_all_GAM_mobular, 
#      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/GAM_Simulated/df_fit_all_GAM_mobular.RData")
# 
# 
# df_spatial_field_GAM_mobular <- bind_rows(lista_spatial_field)  %>%
#   mutate(prevalence = factor(prevalence, 
#                              levels = sort(unique(prevalence),
#                                            decreasing = TRUE)))
# 
# save(df_spatial_field_GAM_mobular, 
#      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/GAM_Simulated/df_spatial_field_GAM_mobular.RData")


# PARTIAL PLOTS -----------------------------------------------------------

mod_05 <- lista_models[["0.5"]]

##Modelo de 0.5


#Volvemos a los valores reales de las var ambientales antes de escalar
mu    <- scaling_params %>% select(ends_with("_mean"))
sdv   <- scaling_params %>% select(ends_with("_sd"))

#Intercepto
beta0 <- mod_05$summary.fixed["intercept", "mean"]

make_partial_plot <- function(var_short, mu, sdv, mod_obj, n_groups = 20, plot_on = c("prob","logit")) {
  plot_on <- match.arg(plot_on)
  
  #TABLA RESUMEN DEL EFECTO ALEATORIO 
  nm_rnd <- paste0("inla.group(", var_short, ", n = ", n_groups, ")") #nombre con el que lo guarda inla
  re_tbl <- mod_obj$summary.random[[nm_rnd]]
  
  # los valores escalados originales 
  x_scaled <- df_monthly_1deg_scale[[var_short]]
  
  #intervalores en los valores escalados para simular lo que ha hecho inla (20 intervalos continuos)
  brks <- seq(min(x_scaled, na.rm = TRUE),
              max(x_scaled, na.rm = TRUE),
              length.out = n_groups + 1)
  centers_scaled <- (brks[-1] + brks[-length(brks)]) / 2 #calculo de la media 
  
  #Devolvemos a la escala original 
  mu_i <- mu[[ paste0(var_short, "_mean") ]]
  sd_i <- sdv[[ paste0(var_short, "_sd") ]]
  centers_orig <- centers_scaled * sd_i + mu_i
  
  #datafrmae en escala original
  dfp <- data.frame(
    x_orig    = centers_orig, #centra cada uno de los 20 puntos desescalados 
    effect_m  = re_tbl$mean, #media del efecto aleatorio en ese bin
    lower_m   = re_tbl$`0.025quant`, #cuantil inferior
    upper_m   = re_tbl$`0.975quant`#cuantil superior
  )
  # lo que buscamos = intercepto + efecto
  dfp$logit_par <- beta0 + dfp$effect_m #logit parcial sumando la pendiente
  dfp$logit_lo  <- beta0 + dfp$lower_m
  dfp$logit_hi  <- beta0 + dfp$upper_m
  
  #trasformacion si hemos seleccionado la option "prob" por que es la porbabilidad de 0 a 1
  #plogis aplica la funcion logistica inversa
  if (plot_on == "prob") {
    dfp$y_mean <- plogis(dfp$logit_par)
    dfp$y_lo   <- plogis(dfp$logit_lo)
    dfp$y_hi   <- plogis(dfp$logit_hi)
    y_label    <- "Probabilidad parcial"
  } else {
    dfp$y_mean <- dfp$logit_par
    dfp$y_lo   <- dfp$logit_lo
    dfp$y_hi   <- dfp$logit_hi
    y_label    <- "Logit parcial"
  }
  
  # graficamos 
  magma_col <- magma(1)
  p <- ggplot(dfp, aes(x = x_orig, y = y_mean)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = y_lo, ymax = y_hi),
                alpha = 0.3, fill = magma_col) +
    labs(
      x = var_short,
      y = "Posterior mean",
    ) +
    theme_classic(base_size = 11) +
    theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
           legend.position = "top",
           axis.line       = element_line(color = "black", linewidth = 0.4),
           axis.ticks      = element_line(color = "black", linewidth = 0.3),
           panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
           panel.grid      = element_blank(),
           text = element_text(family = "Times New Roman")
    )
  
  return(p)
}

p1 <- make_partial_plot("mnkc_epi",          mu, sdv, mod_05, plot_on = "prob")
p2 <- make_partial_plot("sst_front",         mu, sdv, mod_05, plot_on = "prob")
p3 <- make_partial_plot("chl_integrated",    mu, sdv, mod_05, plot_on = "prob")
p4 <- make_partial_plot("termocline_intensity", mu, sdv, mod_05, plot_on = "prob")


mod_001 <- lista_models[["0.01"]]

beta0 <- mod_001$summary.fixed["intercept", "mean"]

make_partial_plot <- function(var_short, mu, sdv, mod_obj, n_groups = 20, plot_on = c("prob","logit")) {
  plot_on <- match.arg(plot_on)
  
  #TABLA RESUMEN DEL EFECTO ALEATORIO 
  nm_rnd <- paste0("inla.group(", var_short, ", n = ", n_groups, ")") #nombre con el que lo guarda inla
  re_tbl <- mod_obj$summary.random[[nm_rnd]]
  
  # los valores escalados originales 
  x_scaled <- df_monthly_1deg_scale[[var_short]]
  
  #intervalores en los valores escalados para simular lo que ha hecho inla (20 intervalos continuos)
  brks <- seq(min(x_scaled, na.rm = TRUE),
              max(x_scaled, na.rm = TRUE),
              length.out = n_groups + 1)
  centers_scaled <- (brks[-1] + brks[-length(brks)]) / 2 #calculo de la media 
  
  #Devolvemos a la escala original 
  mu_i <- mu[[ paste0(var_short, "_mean") ]]
  sd_i <- sdv[[ paste0(var_short, "_sd") ]]
  centers_orig <- centers_scaled * sd_i + mu_i
  
  #datafrmae en escala original
  dfp <- data.frame(
    x_orig    = centers_orig, #centra cada uno de los 20 puntos desescalados 
    effect_m  = re_tbl$mean, #media del efecto aleatorio en ese bin
    lower_m   = re_tbl$`0.025quant`, #cuantil inferior
    upper_m   = re_tbl$`0.975quant`#cuantil superior
  )
  # lo que buscamos = intercepto + efecto
  dfp$logit_par <- beta0 + dfp$effect_m #logit parcial sumando la pendiente
  dfp$logit_lo  <- beta0 + dfp$lower_m
  dfp$logit_hi  <- beta0 + dfp$upper_m
  
  #trasformacion si hemos seleccionado la option "prob" por que es la porbabilidad de 0 a 1
  #plogis aplica la funcion logistica inversa
  if (plot_on == "prob") {
    dfp$y_mean <- plogis(dfp$logit_par)
    dfp$y_lo   <- plogis(dfp$logit_lo)
    dfp$y_hi   <- plogis(dfp$logit_hi)
    y_label    <- "Probabilidad parcial"
  } else {
    dfp$y_mean <- dfp$logit_par
    dfp$y_lo   <- dfp$logit_lo
    dfp$y_hi   <- dfp$logit_hi
    y_label    <- "Logit parcial"
  }
  
  # graficamos 
  magma_col <- magma(1)
  p <- ggplot(dfp, aes(x = x_orig, y = y_mean)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = y_lo, ymax = y_hi),
                alpha = 0.3, fill = magma_col) +
    labs(
      x = var_short,
      y = "Posterior mean",
    ) +
    theme_classic(base_size = 11) +
    theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
           legend.position = "top",
           axis.line       = element_line(color = "black", linewidth = 0.4),
           axis.ticks      = element_line(color = "black", linewidth = 0.3),
           panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
           panel.grid      = element_blank(),
           text = element_text(family = "Times New Roman")
    )
  
  return(p)
}


p6 <- make_partial_plot("mnkc_epi",          mu, sdv, mod_001, plot_on = "prob")
p7 <- make_partial_plot("sst_front",         mu, sdv, mod_001, plot_on = "prob")
p8 <- make_partial_plot("chl_integrated",    mu, sdv, mod_001, plot_on = "prob")
p9 <- make_partial_plot("termocline_intensity", mu, sdv, mod_001, plot_on = "prob")





##Month no se ha ploteado 


mod_05 <- lista_models[["0.5"]]

beta0 <- mod_05$summary.fixed["intercept", "mean"]

#Para no crear cada partial plot a mano, hacemos una función

make_partial_plot <- function(var_short, mod_obj, n_groups = 20, plot_on = c("prob","logit")) {
  plot_on <- match.arg(plot_on)
  
  nm_rnd <- names(mod_obj$summary.random)[grep(var_short, names(mod_obj$summary.random))]
  re_tbl <- mod_obj$summary.random[[nm_rnd]]
  
  #datafrmae en escala original
  dfp <- data.frame(
    x_orig = re_tbl[,1],
    effect_m  = re_tbl$mean, #media del efecto aleatorio en ese bin
    lower_m   = re_tbl$`0.025quant`, #cuantil inferior
    upper_m   = re_tbl$`0.975quant`#cuantil superior
  )
  # lo que buscamos = intercepto + efecto
  dfp$logit_par <- beta0 + dfp$effect_m #logit parcial sumando la pendiente
  dfp$logit_lo  <- beta0 + dfp$lower_m
  dfp$logit_hi  <- beta0 + dfp$upper_m
  
  #trasformacion si hemos seleccionado la option "prob" por que es la porbabilidad de 0 a 1
  #plogis aplica la funcion logistica inversa
  if (plot_on == "prob") {
    dfp$y_mean <- plogis(dfp$logit_par)
    dfp$y_lo   <- plogis(dfp$logit_lo)
    dfp$y_hi   <- plogis(dfp$logit_hi)
    y_label    <- "Probabilidad parcial"
  } else {
    dfp$y_mean <- dfp$logit_par
    dfp$y_lo   <- dfp$logit_lo
    dfp$y_hi   <- dfp$logit_hi
    y_label    <- "Logit parcial"
  }
  
  # graficamos 
  magma_col <- magma(1)
  
  p <- ggplot(dfp, aes(x = x_orig, y = y_mean)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = y_lo, ymax = y_hi),
                alpha = 0.3, fill = magma_col) +
    labs(
      x = var_short,
      y = "Posterior mean",
    ) +
    theme_classic(base_size = 11) +
    theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
           legend.position = "top",
           axis.line       = element_line(color = "black", linewidth = 0.4),
           axis.ticks      = element_line(color = "black", linewidth = 0.3),
           panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
           panel.grid      = element_blank(),
           text = element_text(family = "Times New Roman")
    )
  
  return(p)
}


p_month05 <- make_partial_plot("month", mod_05, plot_on = "prob")


mod_001 <- lista_models[["0.1"]]

beta0 <- mod_001$summary.fixed["intercept", "mean"]

#Para no crear cada partial plot a mano, hacemos una función

make_partial_plot <- function(var_short, mod_obj, n_groups = 20, plot_on = c("prob","logit")) {
  plot_on <- match.arg(plot_on)
  
  nm_rnd <- names(mod_obj$summary.random)[grep(var_short, names(mod_obj$summary.random))]
  re_tbl <- mod_obj$summary.random[[nm_rnd]]
  
  #datafrmae en escala original
  dfp <- data.frame(
    x_orig = re_tbl[,1],
    effect_m  = re_tbl$mean, #media del efecto aleatorio en ese bin
    lower_m   = re_tbl$`0.025quant`, #cuantil inferior
    upper_m   = re_tbl$`0.975quant`#cuantil superior
  )
  # lo que buscamos = intercepto + efecto
  dfp$logit_par <- beta0 + dfp$effect_m #logit parcial sumando la pendiente
  dfp$logit_lo  <- beta0 + dfp$lower_m
  dfp$logit_hi  <- beta0 + dfp$upper_m
  
  #trasformacion si hemos seleccionado la option "prob" por que es la porbabilidad de 0 a 1
  #plogis aplica la funcion logistica inversa
  if (plot_on == "prob") {
    dfp$y_mean <- plogis(dfp$logit_par)
    dfp$y_lo   <- plogis(dfp$logit_lo)
    dfp$y_hi   <- plogis(dfp$logit_hi)
    y_label    <- "Probabilidad parcial"
  } else {
    dfp$y_mean <- dfp$logit_par
    dfp$y_lo   <- dfp$logit_lo
    dfp$y_hi   <- dfp$logit_hi
    y_label    <- "Logit parcial"
  }
  
  # graficamos 
  magma_col <- magma(1)
  
  p <- ggplot(dfp, aes(x = x_orig, y = y_mean)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = y_lo, ymax = y_hi),
                alpha = 0.3, fill = magma_col) +
    labs(
      x = var_short,
      y = "Posterior mean",
    ) +
    theme_classic(base_size = 11) +
    theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
           legend.position = "top",
           axis.line       = element_line(color = "black", linewidth = 0.4),
           axis.ticks      = element_line(color = "black", linewidth = 0.3),
           panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
           panel.grid      = element_blank(),
           text = element_text(family = "Times New Roman")
    )
  
  return(p)
}


p_month001 <- make_partial_plot("month", mod_001, plot_on = "prob")



titulo_columna <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "High", size = 6, fontface = "bold") +
  theme_void()

titulo_columna2 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Extremely low", size = 6, fontface = "bold") +
  theme_void()

p_final <-  (titulo_columna | titulo_columna2) / (p1 / p2 / p3 / p4 / p_month05 | p6 / p7 / p8 /p9/p_month001 ) +
  plot_layout(heights = c(0.05, 1))

p_final




