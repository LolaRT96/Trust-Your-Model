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


# Load necessary libraries
#Libraries
library(sp)
library(tidyr)
library(geoR)
library(ggridges)
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


#GAM-INLA ------------------------------------------------------

df_simulacion_fit <- df_simulacion_fit %>%
  dplyr::select(x, y, temp, time, bathy, pres)

glimpse(df_simulacion_fit)

# PERFORMANCE MODELOS 100 ITERACIONES  ------------------------------------

#Orden 

#Projector matrix
#SPDE + Spatial Field 
#Stack
#Formula 
#Model

#Set the mesh

#MESH

spp <- df_simulacion_fit[, c("x", "y", "pres")]
loc <- cbind(spp$x, spp$y)
loc <- as.data.frame(loc)
colnames(loc) <- c("x", "y")
coordinates(loc) <- ~x + y
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

convhull <- inla.nonconvex.hull(loc) 
mesh <- INLA::inla.mesh.2d(boundary = convhull,
                           max.edge = c(8,15),
                           cutoff = 0.4)

mesh$n 

plot(mesh)

#Definios los hiperpriors de los efectos aleatorios (random walk)
hyper_pc <- list(prec = list(prior = "pc.prec", param = c(3, 0.05)))

#cada iteraccion cambia el campo espacial por lo que debemos generar una formula para 
#que use cada spde en cada iteraccion 

formula_generator <- function(spde) {
  y ~ -1 + intercept +
    f(inla.group(temp, n = 20), model = "rw2", hyper = hyper_pc) +
    f(inla.group(bathy, n = 20), model = "rw2", hyper = hyper_pc) +
    f(time, model = "seasonal", season.length = 12) +
    f(spatial.field, model = spde)
}

#Bootstrapping de prevalencia 

sensitivity_prevalence_inla <- function(data, 
                                        y                 = "pres",
                                        coords            = c("x", "y"),
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
        temp, bathy, time #viene de lo que se defina en data original
      )
    
    # Pasamos a todos los argumentos que necesitamos para inla
    #MESH Y SPDE + INDICE 
    coords_mat <- as.matrix(dat_i[, coords]) #matriz para construir el mesh 
    hull       <- inla.nonconvex.hull(coords_mat) 
    mesh_i     <- inla.mesh.2d(boundary = hull, max.edge = c(8, 15), cutoff = 0.4) #conservamos los argumentos
    spde_i     <- inla.spde2.pcmatern(mesh_i,
                                      prior.range = c(2,0.1),
                                      prior.sigma = c(1,0.1))
    A_i        <- inla.spde.make.A(mesh_i, loc = coords_mat)
    s.index    <- inla.spde.make.index("spatial.field", spde_i$n.spde)
    
    stack_i <- inla.stack(
      data    = list(y = dat_i[[y]]), #y=variable respuesta 
      A       = list(A_i, 1),
      effects = list(
        s.index,
        data.frame(
          intercept = 1,
          temp      = dat_i$temp,
          bathy     = dat_i$bathy,
          time      = dat_i$time
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


prevalence <- c(0.50, 0.30, 0.25, 0.20, 0.15,
                0.10, 0.05, 0.025, 0.01)


all_res <- lapply(prevalence, function(p) {
  sensitivity_prevalence_inla(
    data               = df_simulacion_fit,
    y                  = "pres",
    coords             = c("x","y"),
    formula_generator  = formula_generator,
    hyper_pc           = hyper_pc,
    n_iter             = 100,
    prevalence         = p,
    verbose            = TRUE
  )
})

df_all <- bind_rows(all_res)

# Convertir prevalencia a factor ordenado para el eje

df_iteracciones_gam_v2 <- df_all

save(df_iteracciones_gam_v2, 
     file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/GAM_Simulated/df_iteracciones_gam_v2.RData")

summary(df_iteracciones_gam_v2)

df_iteracciones_gam_v2$Prevalence <- as.factor(df_iteracciones_gam_v2$Prevalence)

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
p_auc  <- plot_raincloud(df_iteracciones_gam , "AUC")
p_sens <- plot_raincloud(df_iteracciones_gam, "Sensitivity")
p_spec <- plot_raincloud(df_iteracciones_gam, "Specificity")
p_tss  <- plot_raincloud(df_iteracciones_gam, "TSS")
p_thr  <- plot_raincloud(df_iteracciones_gam, "Cutoff")
p_log <- plot_raincloud(df_iteracciones_gam, "LogLoss")

windows()
(p_auc | p_sens | p_thr) / (p_spec | p_tss | p_log) 


# Un modelo por cada ratio ------------------------------------------------

vars_keep <- c("pres","x","y",
               "time","bathy","temp")


loc_pred <- as.matrix(df_simulacion_predict[, c("x","y")])
n_pred <- nrow(loc_pred)
A.pred   <- inla.spde.make.A(mesh, loc = loc_pred)

#Definios los hiperpriors de los efectos aleatorios (random walk)
hyper_pc <- list(prec = list(prior = "pc.prec", param = c(3, 0.05)))

spde <- inla.spde2.pcmatern(mesh, prior.range=c(2, 0.1), prior.sigma=c(1, 0.1))
s.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

f <- y ~ -1 + intercept +
  f(inla.group(temp, n = 20), model = "rw2", hyper = hyper_pc) +
  f(inla.group(bathy, n = 20), model = "rw2", hyper = hyper_pc) +
  f(time, model = "seasonal", season.length = 12) +
  f(spatial.field, model = spde)


A <- inla.spde.make.A(mesh, loc = loc) #Para la inferencia

# Inicializar listas
lista_probs         <- list()
lista_fit           <- list()
lista_spatial_field <- list()
lista_models <- list()

# Ratios a probar
prevalences <- c(0.5, 0.2,0.1, 0.01)

for (p in prevalences) {
  message("Procesando prevalencia = ", p)
  
  # 1. Submuestra
  data_pres <- df_simulacion_fit %>% filter(pres == 1)
  data_abs  <- df_simulacion_fit %>% filter(pres == 0)
  n_pres    <- nrow(data_pres)
  n_abs     <- round(n_pres * (1 - p) / p)
  sampled_abs <- data_abs %>% slice_sample(n = n_abs)
  
  data_ratio <- bind_rows(data_pres, sampled_abs) %>%
    drop_na(all_of(vars_keep))
  
  # 2. Matriz A para inferencia
  # 2. Matriz A para inferencia (recalcular en cada ratio)
  loc_ratio <- as.matrix(data_ratio[, c("x","y")])
  A.inf     <- inla.spde.make.A(mesh, loc = loc_ratio)
  
  # 3. Stacks
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
      df_simulacion_predict  %>% transmute(
        intercept = 1,
        temp, time, bathy
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
  
  # 5. Predicciones sobre el grid (stack.pred)
  idx_pred <- inla.stack.index(stack.full, tag = "pred")$data
  res_pred <- data.frame(
    x = df_simulacion_predict$x,
    y  = df_simulacion_predict$y,
    link_mean = mod$summary.fitted.values[idx_pred, "mean"],
    link_sd   = mod$summary.fitted.values[idx_pred, "sd"],
    lower_q   = mod$summary.fitted.values[idx_pred, "0.025quant"],
    upper_q   = mod$summary.fitted.values[idx_pred, "0.975quant"]
  ) %>%
    mutate(
      prob_mean  = plogis(link_mean),
      prob_lower = plogis(lower_q),
      prob_upper = plogis(upper_q),
      prob_sd    = link_sd * prob_mean * (1 - prob_mean),
      # odds       = prob_mean / (1 - prob_mean),
      # ratio_glob = nrow(data_pres) / nrow(sampled_abs),
      # favorability = odds / (ratio_glob + odds),
      # ratio = as.character(p)
    )
  
  # 6. Predicciones sobre datos de entrenamiento (stack.fit)
  idx_fit <- inla.stack.index(stack.full, tag = "fit")$data
  res_fit <- data.frame(
    prob_mean = plogis(mod$summary.fitted.values[idx_fit, "mean"]),
    pres  = data_ratio$pres,
    ratio     = as.character(p)
  )
  
  proj <- inla.mesh.projector(
    mesh,
    xlim = range(mesh$loc[,1]),
    ylim = range(mesh$loc[,2]),
    dims = c(300, 300)
  )
  
  field_mean <- inla.mesh.project(proj,
                                  field = mod$summary.random$spatial.field$mean)
  
  df_field <- reshape2::melt(field_mean)
  names(df_field) <- c("x_id", "y_id", "mean")
  df_field$x <- proj$x[df_field$x_id]
  df_field$y <- proj$y[df_field$y_id]
  df_field$ratio <- as.character(p)
  
  # 8. Guardar resultados
  lista_probs[[as.character(p)]]         <- res_pred
  lista_fit[[as.character(p)]]           <- res_fit
  lista_spatial_field[[as.character(p)]] <- df_field
  lista_models[[as.character(p)]]        <- mod 
  
}

#Combinar todos los resultados

df_prob_ratios_sim <- bind_rows(lista_probs) %>%
  mutate(ratio = factor(ratio, levels = sort(unique(prevalences), decreasing = TRUE)))

head(df_prob_ratios_sim)

# save(df_prob_ratios_sim, 
#      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/GAM_Simulated/df_prob_ratios_sim.RData")


df_fit_all_sim <- bind_rows(lista_fit) %>%
  mutate(ratio = factor(ratio, levels = sort(unique(ratio), decreasing = TRUE)))
head(df_fit_all_sim)

# save(df_fit_all_sim, 
#      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/GAM_Simulated/df_fit_all_sim.RData")


df_spatial_field_sim <- bind_rows(lista_spatial_field) %>%
  mutate(ratio = factor(ratio, levels = sort(unique(ratio), decreasing = TRUE)))
head(df_spatial_field_sim)

# save(df_spatial_field_sim, 
#      file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/GAM_Simulated/df_spatial_field_sim.RData")


df_bins <- df_fit_all %>%
  group_by(ratio) %>%
  mutate(bin = ntile(prob_mean, 10)) %>%
  group_by(ratio, bin) %>%
  summarise(
    mean_pred = mean(prob_mean, na.rm = TRUE),
    obs_freq  = mean(pres, na.rm = TRUE),
    .groups   = "drop"
  )

df_bins <- df_bins %>%
  mutate(Prevalence = factor(ratio, levels = c("0.5","0.2","0.1","0.01"))) %>%
  filter(Prevalence %in% levels(Prevalence))

windows();ggplot(df_bins, aes(x = mean_pred, y = obs_freq)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_line(color = "black") +
  geom_point(size = 2, shape = 21, fill = "gray70", color = "black") +
  facet_wrap(~ratio, ncol = 5, scale = "free_x") +
  theme_classic() +
  labs(
    title = "Reliability diagram por ratio de prevalencia (GAM-INLA)",
    x = "Probabilidad media predicha (por bin)",
    y = "Frecuencia observada"
  )


#Maps prob mean

windows();ggplot(df_prob_ratios, aes(x = x, y = y, fill = prob_mean)) +
  geom_tile() +
  facet_wrap(~ratio, ncol = 5) +
  scale_fill_gradientn(
    name   = "Probabilidad",
    colors = c(
      "#08589e", "#2b8cbe", "#74a9cf", "#bfd3e6",
      "#fdbb84", "#e34a33", "#a50f15", "#5b0100"
    ),
    values = scales::rescale(c(0, 0.15, 0.30, 0.45, 0.55, 0.70, 0.85, 1)),
    guide  = "colourbar"
  ) +
  labs(
    title = "Mapas de probabilidad predicha por ratio de prevalencia (GAM-INLA)",
    x = "Longitud",
    y = "Latitud"
  ) +
  theme_minimal(base_size = 11)


#Spatial structure 

windows();ggplot(df_spatial_field, aes(x = x, y = y, fill = mean)) +
  geom_tile() +
  facet_wrap(~ratio, ncol = 5) +
  scale_fill_viridis_c(
    option = "C",
    name = "Campo espacial",
    guide = "colourbar"
  )  +
  labs(
    title = "Campo espacial estimado por ratio de prevalencia (GAM-INLA)",
    x = "Longitud",
    y = "Latitud"
  ) +
  theme_minimal(base_size = 11)

windows(); ggplot(df_spatial_field, aes(x = x, y = y, fill = mean)) +
  geom_tile() +
  facet_wrap(~ratio, ncol = 5) +
  scale_fill_viridis_c(option = "C", name = "Campo espacial", guide = "colourbar") +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    title = "Campo espacial estimado por ratio de prevalencia (GAM-INLA)",
    x = "Longitud",
    y = "Latitud"
  ) +
  theme_minimal(base_size = 11)

##densidad de probabilidades
df_prob_ratios %>%
  filter(!is.na(prob_mean)) %>%
  ggplot(aes(x = prob_mean, y = ratio, fill = ratio)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Probabilidad predicha", y = "Prevalencia",
       title = "Ridgelines de probabilidades por prevalencia") +
  theme_classic()


# PARTIAL PLOTS -------------------------------------------------

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

p1 <- make_partial_plot("temp", mod_05, plot_on = "prob")
p1
p2 <- make_partial_plot("time", mod_05, plot_on = "prob")
p2
p3 <- make_partial_plot("bathy", mod_05, plot_on = "prob")





mod_001 <- lista_models[["0.01"]]

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

p4 <- make_partial_plot("temp", mod_001, plot_on = "prob")

p5 <- make_partial_plot("time", mod_001, plot_on = "prob")

p6 <- make_partial_plot("bathy", mod_001, plot_on = "prob")


titulo_columna <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "High", size = 6, fontface = "bold") +
  theme_void()

titulo_columna2 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Extremely low", size = 6, fontface = "bold") +
  theme_void()

p_final <-  (titulo_columna | titulo_columna2) / (p1 / p2 / p3  | p4 / p6 / p5 ) +
  plot_layout(heights = c(0.05, 1))

p_final
