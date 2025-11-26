#Authors: Maria Dolores Riesgo (IEO-CSIC)
#Research paper: "Trust your model"
#General objective: fit the simulated data with GAM-INLA models 

#R version 4.4.2

#If necessary 

#R version 4.4.2

#If necessary 

rm(list=ls(all=TRUE)) 

# Load necessary libraries
#Libraries

library(INLA)
library(dismo)
library(hSDM)
library(spdep)
library(fields)
library(gridExtra)
library(ggplot2)
library(reshape)
library(patchwork)
library(viridis)
library(hrbrthemes)
library(INLA)
library(pROC)

inla.setOption(num.threads = 4)

set.seed(123456789)
print(.Random.seed[1:4])

# LOAD SIMULATION DATA -----------------------------------------

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

df_simulacion_fit <- df_simulacion_fit %>%
  dplyr::select(x, y, temp, time, bathy, pres)

glimpse(df_simulacion_fit)

# Assessing model performance with bootstrap resampling  ------------------------------------

#INLA coding framework 
#Build the mesh
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
convhull <- inla.nonconvex.hull(loc) 
mesh <- INLA::inla.mesh.2d(boundary = convhull,
                           max.edge = c(8,15),
                           cutoff = 0.4)

mesh$n 

plot(mesh)

#Random Walk hyperpriors 
hyper_pc <- list(prec = list(prior = "pc.prec", param = c(3, 0.05)))

#each iteration changes the spatial field, so we must generate a formula for 
#each spde to use in each iteration 

formula_generator <- function(spde) {
  y ~ -1 + intercept +
    f(inla.group(temp, n = 20), model = "rw2", hyper = hyper_pc) +
    f(inla.group(bathy, n = 20), model = "rw2", hyper = hyper_pc) +
    f(time, model = "seasonal", season.length = 12) +
    f(spatial.field, model = spde)
}


bootstrap_INLA <- function(data,
                           y  = "pres",
                           coords = c("x", "y"),
                           formula_generator,
                           hyper_pc,
                           n_iter = 100,
                           prevalence = 0.5,
                           verbose = TRUE) {
  
  data_pres <- filter(data, !!sym(y) == 1) #Subset of presences
  data_abs  <- filter(data, !!sym(y) == 0) #Subset of absence
  n_pres    <- nrow(data_pres) #total data
  # undersampling the zeros
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
  
  for (i in seq_len(n_iter)) { 
    t0 <- Sys.time() #computing itme of each iteration 
    if (verbose) cat("Prev =", prevalence, "| Iter", i, "of", n_iter, "...\n") #number of the iteration
    
    set.seed(i)
    
    sampled_abs <- slice_sample( 
      data_abs, 
      n       = n_abs_req,
      replace = (n_abs_req > nrow(data_abs)) 
    )
    
    dat_i <- bind_rows(data_pres, sampled_abs) %>%  
      drop_na(
        !!sym(y), all_of(coords), 
        temp, bathy, time 
      )
    
    #MESH + SPDE 
    coords_mat <- as.matrix(dat_i[, coords])
    hull       <- inla.nonconvex.hull(coords_mat) 
    mesh_i     <- inla.mesh.2d(boundary = hull, max.edge = c(8, 15), cutoff = 0.4) 
    spde_i     <- inla.spde2.pcmatern(mesh_i,
                                      prior.range = c(2,0.1),
                                      prior.sigma = c(1,0.1))
    A_i        <- inla.spde.make.A(mesh_i, loc = coords_mat)
    s.index    <- inla.spde.make.index("spatial.field", spde_i$n.spde)
    
    stack_i <- inla.stack(
      data    = list(y = dat_i[[y]]),
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
    
    idx_obs <- inla.stack.index(stack_i, tag = "model")$data 
    fitted  <- mod_i$summary.fitted.values$mean[idx_obs] 
    true    <- dat_i[[y]] 
    roc_obj <- roc(response = true, predictor = fitted, quiet = TRUE)
    best    <- coords(roc_obj, x = "best", best.method = "youden",
                      ret = c("threshold","sensitivity","specificity"))
    thr  <- best["threshold"]; sens <- best["sensitivity"]; spec <- best["specificity"]
    tss  <- sens + spec - 1
    aucv <- as.numeric(auc(roc_obj))
    brier <- mean((fitted - true)^2)
    
    results[i, c("AUC","Sensitivity","Specificity","TSS","Cutoff","BrierScore", "Prevalence")] 
    <- c(aucv, sens, spec, tss, thr, logloss, brier, prevalence)
    
    if (verbose) {
      dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
      cat("Duration iter", i, ":", dt, "mins\n\n")
    }
  }
  
  results$Prevalence <- prevalence 
  return(results) 
}


prevalence <- c(0.50, 0.20, 
                0.10, 0.01)

all_res <- lapply(prevalence, function(p) {
  bootstrap_INLA(
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

plot <- function(df, metric) {
  df <- df %>%
   
    mutate(Prevalence = fct_reorder(as.factor(Prevalence),
                                    as.numeric(as.character(Prevalence)),
                                    .desc = TRUE))
  
  ggplot(df, aes(x = Prevalence, y = .data[[metric]], fill = Prevalence)) +
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
p_auc  <- plot(df_all , "AUC")
p_sens <- plot(df_all, "Sensitivity")
p_spec <- plot(df_all, "Specificity")
p_tss  <- plot(df_all, "TSS")
p_thr  <- plot(df_all, "Cutoff")
p_log <- plot(df_all, "BrierScore")

# MODELS AND PREDICTION --------------------------------

vars_keep <- c("pres","x","y",
               "time","bathy","temp")


loc_pred <- as.matrix(df_simulacion_predict[, c("x","y")])
n_pred <- nrow(loc_pred)
A.pred   <- inla.spde.make.A(mesh, loc = loc_pred)


hyper_pc <- list(prec = list(prior = "pc.prec", param = c(3, 0.05)))

spde <- inla.spde2.pcmatern(mesh, prior.range=c(2, 0.1), prior.sigma=c(1, 0.1))
s.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

f <- y ~ -1 + intercept +
  f(inla.group(temp, n = 20), model = "rw2", hyper = hyper_pc) +
  f(inla.group(bathy, n = 20), model = "rw2", hyper = hyper_pc) +
  f(time, model = "seasonal", season.length = 12) +
  f(spatial.field, model = spde)


A <- inla.spde.make.A(mesh, loc = loc) #Para la inferencia

lista_probs         <- list()
lista_fit           <- list()
lista_spatial_field <- list()
lista_models <- list()

prevalences <- c(0.5, 0.2,0.1, 0.01)

for (p in prevalences) {
  message("Prevalence = ", p)
  
  # 1. Submuestra
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
      df_simulacion_predict  %>% transmute(
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
      ratio = as.character(p)
    )

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
  
  lista_probs[[as.character(p)]]         <- res_pred
  lista_fit[[as.character(p)]]           <- res_fit
  lista_spatial_field[[as.character(p)]] <- df_field
  lista_models[[as.character(p)]]        <- mod 
  
}


df_probabilities <- bind_rows(lista_probs) %>%
  mutate(ratio = factor(ratio, levels = sort(unique(prevalences), decreasing = TRUE)))

head(df_probabilities)

df_fit_probabilities <- bind_rows(lista_fit) %>%
  mutate(ratio = factor(ratio, levels = sort(unique(ratio), decreasing = TRUE)))

head(df_fit_probabilities)


df_spatial_field_sim <- bind_rows(lista_spatial_field) %>%
  mutate(ratio = factor(ratio, levels = sort(unique(ratio), decreasing = TRUE)))
head(df_spatial_field_sim)


#Reliability diagram 
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
  theme_classic() 


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
  theme_minimal(base_size = 11)


#Spatial structure 

windows();ggplot(df_spatial_field, aes(x = x, y = y, fill = mean)) +
  geom_tile() +
  facet_wrap(~ratio, ncol = 5) +
  scale_fill_viridis_c(
    option = "C",
    name = "Spatial field",
    guide = "colourbar"
  )  +
  theme_minimal(base_size = 11)


# PARTIAL PLOTS -------------------------------------------------

mod_05 <- lista_models[["0.5"]]

beta0 <- mod_05$summary.fixed["intercept", "mean"]


make_partial_plot <- function(var_short, mod_obj, n_groups = 20, plot_on = c("prob","logit")) {
  plot_on <- match.arg(plot_on)
  
  nm_rnd <- names(mod_obj$summary.random)[grep(var_short, names(mod_obj$summary.random))]
  re_tbl <- mod_obj$summary.random[[nm_rnd]]
  
  dfp <- data.frame(
    x_orig = re_tbl[,1],
    effect_m  = re_tbl$mean, 
    lower_m   = re_tbl$`0.025quant`, 
    upper_m   = re_tbl$`0.975quant`
  )
  
  dfp$logit_par <- beta0 + dfp$effect_m
  dfp$logit_lo  <- beta0 + dfp$lower_m
  dfp$logit_hi  <- beta0 + dfp$upper_m
  
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

make_partial_plot <- function(var_short, mod_obj, n_groups = 20, plot_on = c("prob","logit")) {
  plot_on <- match.arg(plot_on)
  
  nm_rnd <- names(mod_obj$summary.random)[grep(var_short, names(mod_obj$summary.random))]
  re_tbl <- mod_obj$summary.random[[nm_rnd]]
  
  dfp <- data.frame(
    x_orig = re_tbl[,1],
    effect_m  = re_tbl$mean, 
    lower_m   = re_tbl$`0.025quant`,
    upper_m   = re_tbl$`0.975quant`
  )
  
  dfp$logit_par <- beta0 + dfp$effect_m 
  dfp$logit_lo  <- beta0 + dfp$lower_m
  dfp$logit_hi  <- beta0 + dfp$upper_m
  
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
