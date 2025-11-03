#Marine BEACON Project (PhD Chapter One)
#Authors: Lola Riesgo 
#General objetive: evaluacion del poder predictivo de los modelos en el cambio de los ratios con datos 
#SIMULADOS 

#1. Generar los datos simulados (mismos que en GAM-INLA scripts)
#Datos prediccion 
#Datos fitted models 

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


#Set the seed

set.seed(123456789)
print(.Random.seed[1:4])

# AREA DE ESTUDIO (GRID) -------------------------------------------------------

campo <- function(c1, c2, c3, c4) {
  xy <- expand.grid(
    seq(c1, c2, length.out = 100),
    seq(c3, c4, length.out = 100)
  )
  cbind(xy[, 1], xy[, 2])
}

c1 <- 0; c2 <- 100; c3 <- 0; c4 <- 100 

loc_xy <- campo(c1, c2, c3, c4) 


# DEFINIMOS EL MESH PARA PROYECTAR LA ESTRUCTURA ESPACIAL ----------------------

with_seed(123456789, {mesh2d_sim <- fm_mesh_2d_inla(loc.domain = loc_xy, 
                                                    max.edge = c(8, 15))  # Creamos la malla 2D
})
# 
# mesh2d_sim <- fm_mesh_2d_inla(loc.domain = loc_xy, 
#                               max.edge = c(8, 15))  # Creamos la malla 2D

windows();plot(mesh2d_sim, main = "Malla", asp = 1, lwd = 0.5)
mesh2d_sim$n

# COVARIABLES 1 AND 2 -----------------------------------------------------

#fixed covariables (batimetría que cambia de Oeste a Este)

n_time <- 12
grid_list <- vector("list", n_time)

for (t in seq_len(n_time)) {
  
  loc_df <- data.frame(
    x <- loc_xy[,1],
    y <- loc_xy[,2]
  )
  
  depth_range <- c(0, 1000)
  
  loc_df$bathy <- depth_range[1] +
    ((loc_df$y - c3) / (c4 - c3)) * diff(depth_range)
  
  grid_list[[t]] <- data_frame (
    x =loc_df$x,
    y= loc_df$y,
    bathy = loc_df$bathy,
    time = t )
}

grid_bathy_df <- bind_rows(grid_list)
head(grid_bathy_df)

windows();ggplot(grid_bathy_df, aes(x =x, y = y, fill = bathy)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis_c(option = "G", direction =-1)+
  labs(
    x = "x",
    y = "y",
    fill = "Bathymetry (m)"  # ← aquí cambias el nombre de la leyenda
  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
         legend.position   = "right",
         axis.line         = element_line(color = "black", linewidth = 0.4),
         axis.ticks        = element_line(color = "black", linewidth = 0.3),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
         panel.grid        = element_blank(),
         strip.background  = element_blank(),                        # ← quita el fondo de la etiqueta
         strip.text        = element_text(face = "bold", size = 12)  # ← estilo de letra para A/B/C/D
  )

p_bathy <- ggplot(grid_bathy_df, aes(x =x, y = y, fill = bathy)) +
  geom_tile() +
  coord_equal() +
  scale_fill_viridis_c(option = "G", direction =-1)+
  labs(
    x = "x",
    y = "y",
    fill = "Bathymetry (m)"  # ← aquí cambias el nombre de la leyenda
  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
         legend.position   = "right",
         axis.line         = element_line(color = "black", linewidth = 0.4),
         axis.ticks        = element_line(color = "black", linewidth = 0.3),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
         panel.grid        = element_blank(),
         strip.background  = element_blank(),                        # ← quita el fondo de la etiqueta
         strip.text        = element_text(face = "bold", size = 12)  # ← estilo de letra para A/B/C/D
  )


#covariable that vary with time (simulation of temp)

n_time <- 12
temp_list <- vector("list", n_time)

#Función sinoidal con temperatura media de 20ºC
base_temp <- function(x, y) {
  20 + 3 * sin(pi * x / 100) + 2 * cos(pi * y / 100)
}

#Caculamos por fuera las anomalias temproales 

anomalies <- cumsum(rnorm(n_time, mean = 0, sd =0.3))

for (t in seq_len(n_time)) {
  
  loc_df <- data.frame(
    x <- loc_xy[,1],
    y <- loc_xy[,2]
  )
  
  loc_df$base_temp <- base_temp(loc_df$x, loc_df$y)
  
  anomaly_t <- anomalies[t]
  
  temp_list[[t]] <- data_frame (
    x = loc_df$x,
    y = loc_df$y,
    temp = loc_df$base_temp + anomaly_t + 
      rnorm(nrow(loc_df), mean = 0, sd = 0.5),
    time = t )
} 

grid_temp_df <- bind_rows(temp_list)
head(grid_temp_df)

windows();ggplot(grid_temp_df, aes(x =x, y = y, fill = temp)) +
  geom_tile() +
  facet_wrap(~time, ncol = 4) +
  coord_equal() +
  scale_fill_viridis_c(option = "H")+
  labs(
    x="x",
    y="y"
  )

#Evaluacion por unidad de tiempo


temp_time <- grid_temp_df %>%
  group_by(time) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE))

# Gráfico de líneas
p_tem <- ggplot(temp_time, aes(x = time, y = mean_temp)) +
  geom_line(color = "firebrick", size = 1.2) +
  geom_point(color = "firebrick", size = 2) +
  labs(
    x = "Time (unit)",
    y = "Temperature (°C)"
  ) +
  scale_x_continuous(breaks = seq(min(temp_time$time), max(temp_time$time), by = 1)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

library(patchwork)

(p_bathy | p_tem )


# ESTRUCTURA ESPACIAL  --------------------------------

##Estructura espacial aleatoria 
with_seed(123456789, {

#Precision
prec <- 1 / 4
#Rango + varianza 
rho <- 20
sigma <- sqrt(5)
phi <- 0.9 #parametro de control de la dependencia temporal entre dos unidades de tiempo
#cercano a 1: alta persistencia, ruido pequeño, menor variación local
#cercano a 0: patrones redefinidos a cada paso temporal 
#cercano a -1 o negativo: zonas altas en t1 son bajas en t2 y así sucesivamente. Efecto vaivén. 
#PROCESO ARIMA: ut = phi*ut-1 + errort, ∣ϕ∣<1
sigma_eps <- sigma * sqrt(1 - phi^2) #Desviacion del ruido, cada año añadimos un ruido de amplitud determinada 
#para que a pesar de la dependencia temporal el campo ni gane ni pierda con el paso del tiempo

u0_nodes <- fmesher::fm_matern_sample(mesh2d_sim, n = 1, rho = rho, sigma = sigma) #Campo aleatorio gaussiano 
#generamos una realización con función de covaraizna de Matern sobre los nodos de la malla generada 

u0_nodes <- u0_nodes - mean(u0_nodes) #centralizamos a media cero para eliminar desviaciones no deseadas del campo 
#es por teoría, los campos guassiano asumen media 0 

n_time <- 12 #unidades de tiempo 
A_grid <- fm_basis(mesh2d_sim, loc_xy) #calculamos la matriz que asocia a cada nodo de la malla cada punto del grid
#es union del mesh + grid (mesh tiene M nodos y grid tiene G celdas, dimensión de la matrix es MxG, y el conteido son los pesos 
#necesarios para interpolar el valor del campo en la celda i a partir de los valores en los nodos)

#el sigueinte paso es generar 12 dataframes guardando los locs y los valores del campo aleatorio espacial
#lo guardamos en una lista y ya lo unimos al final 
latent_list <-  vector("list", n_time) #una lista por cada unidad de tiempo si lo ponemos como vector facilitamos 
#que despues lo podamos unir en una sola lista
u_prev <- u0_nodes #punto de partida de nodos 

#Bucle 
for(t in seq(n_time)) {
  
  eps_t <-  fm_matern_sample(mesh2d_sim, n=1, rho = rho, sigma = sigma) #generamos el ruido
  eps_t <- eps_t - mean(eps_t) #media a 0
  
  u_t <- phi * u_prev + eps_t
  u_prev <- u_t #resultado en los nodos para tiempo t
  
  latent_t <- drop(A_grid %*% u_t) #interpolamos esos valores u_t a cada punto de grilla 
  
  latent_list[[t]] <- data_frame(
    x =loc_xy[,1],
    y=loc_xy[,2],
    latent = latent_t,
    time = t
  )
  
}

#unimos todo lo que se ha contenido en la lista dentro de un solo df

latent_time_df <- bind_rows(latent_list)

} )

windows();ggplot(latent_time_df, aes(x=x, y=y, fill=latent)) +
  geom_tile()+
  facet_wrap(~time, ncol = 4) +
  coord_equal()+
  scale_fill_viridis_c(option = "D") +
  labs(
    x = "x",
    y = "y",
   fill = "Spatial structure"  # ← aquí cambias el nombre de la leyenda
  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
         legend.position   = "right",
         axis.line         = element_line(color = "black", linewidth = 0.4),
         axis.ticks        = element_line(color = "black", linewidth = 0.3),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
         panel.grid        = element_blank(),
         strip.background  = element_blank(),                        # ← quita el fondo de la etiqueta
         strip.text        = element_text(face = "bold", size = 12)  # ← estilo de letra para A/B/C/D
  )


p_latent <- ggplot(latent_time_df, aes(x=x, y=y, fill=latent)) +
  geom_tile()+
  facet_wrap(~time,  nrow = 3) +
  coord_equal()+
  scale_fill_viridis_c() +
  labs(
    x = "x",
    y = "y",
    fill = "Spatial structure"  # ← aquí cambias el nombre de la leyenda
  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
         legend.position   = "bottom",
         axis.line         = element_line(color = "black", linewidth = 0.4),
         axis.ticks        = element_line(color = "black", linewidth = 0.3),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
         panel.grid        = element_blank(),
         strip.background  = element_blank(),                        # ← quita el fondo de la etiqueta
         strip.text        = element_text(face = "bold", size = 12)  # ← estilo de letra para A/B/C/D
  )



(p_bathy | p_tem ) / (p_latent)

# MUESTREO ESPECIE ALEATORIA  ---------------------------------------------

# pesos para las covariables y el latent
n_time        <- 12
n_pts         <- 3000     #3000 puntos de muestreo
beta_latent <- 4.5
beta_temp   <- 1.5
beta_bathy  <- 0.08
frac_ocup_media <- 0.20  
frac_ocup <- 0.20  
sd_frac_ocup    <- 0.01     # variabilidad (puedes ajustar el valor)

df_grid <- grid_temp_df %>%
  inner_join(latent_time_df, by = c("x","y","time")) %>%
  inner_join(grid_bathy_df,     by = c("x","y","time"))


df_list <- vector("list", n_time)

with_seed(123456789, {

for (t in seq_len(n_time)) {
  
  grid_t <- df_grid %>% filter(time == t)
  
  frac_ocup_t <- min(max(rnorm(1, mean = frac_ocup_media, sd = sd_frac_ocup), 0.10), 0.30)
  
  # Calcula el score de idoneidad
  grid_t <- grid_t %>%
    mutate(score = beta_latent * latent +
             beta_temp   * temp +
             beta_bathy  * bathy)
  
  # Determina el umbral para ocupar exactamente el % deseado
  umbral <- quantile(grid_t$score, probs = 1 - frac_ocup)
  
  grid_t <- grid_t %>%
    mutate(
      pres = ifelse(score > umbral, 1, 0)  # top 20% son presencias
    )
  
  df_list[[t]] <- grid_t
}

df_ocup <- bind_rows(df_list)
})

#Muestreo aleatorio estratificado por clases 

target_ratio_bounds <- c(0.01, 0.025) #que los ratios varien por tiempo entre 0.19 y 0.15
n_time              <- 12 #unidades de tiempo 
n_pts               <- 3000 #puntos a muestrear

# ratio aleatorio por cada t in n_time
ratio_vec <- runif(n_time, min = target_ratio_bounds[1], max = target_ratio_bounds[2])

df_muestreo_list <-  vector("list", n_time)

with_seed(123456789, {
for (t in seq_len(n_time)) {
  
  grid_t <- df_ocup %>% filter(time == t)
  pres_t <- df_ocup %>% filter(pres == 1)
  abs_t <- df_ocup %>% filter(pres == 0)
  
  ratio_t <- ratio_vec[t]
  
  #ausencias presencias de los 3000 puntos para cumplir el ratio fijado
  n_abs  <- floor(n_pts / (1 + ratio_t))
  n_pres <- n_pts - n_abs
  
  #SI no hay suficientes
  
  n_abs  <- min(n_abs,  nrow(abs_t))
  n_pres <- min(n_pres, nrow(pres_t))
  
  #Ahora sampleamos por separado para controlar esos ratios que buscamos 
  muestra_abs  <- abs_t  %>% sample_n(n_abs)
  muestra_pres <- pres_t %>% sample_n(n_pres)
  
  muestra_t <- bind_rows(muestra_abs, muestra_pres) #lo unimos en un solo df
  
  df_muestreo_list[[t]] <- muestra_t #lo guardamos de nuevo en la listae
  
}

df_muestra <- bind_rows(df_muestreo_list)
})

table(df_muestra$pres)

mean(df_muestra$pres)

df_muestra %>% 
  group_by(time) %>% 
  summarise(prev = mean(pres))

df_muestra %>%
  group_by(time) %>%
  summarise(
    n_pres = sum(pres == 1),
    n_abs  = sum(pres == 0),
    ratio  = n_pres / n_abs
  )


df_simulacion_predict <- df_ocup

df_simulacion_fit <- df_muestra

#Filtramos presencias
df_presencias <- df_simulacion_fit %>%
  filter(pres == 1)
df_ausencias <- df_simulacion_fit %>%
  filter(pres == 0)

ggplot(df_ocup, aes(x = x, y = y, fill = factor(pres))) +
  geom_tile(color = NA) +
  coord_equal() +
  scale_fill_manual(values = c("0" = "grey90", "1" = "red"),
                    labels = c("Absence", "Presence")) +
  facet_wrap(~ time, ncol = 4) +
  labs(
     x = "x", 
     y = "y") +
  # geom_tile(
  #   data = df_presencias,
  #  aes(x = x, y = y),
  # inherit.aes = FALSE,
  #  width  = 1,
  #    height = 1,
  #   fill   = NA,
  #    color  = "black",
  #   linewidth = 0.35
  #  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
         legend.position   = "bottom",
         axis.line         = element_line(color = "black", linewidth = 0.4),
         axis.ticks        = element_line(color = "black", linewidth = 0.3),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
         panel.grid        = element_blank(),
         strip.background  = element_blank(),                        # ← quita el fondo de la etiqueta
         strip.text        = element_text(face = "bold", size = 12)  # ← estilo de letra para A/B/C/D
  )


p_spp <- ggplot(df_ocup, aes(x = x, y = y, fill = factor(pres))) +
  geom_tile(color = NA) +
  coord_equal() +
  scale_fill_manual(
    values = c("0" = "grey90", "1" = "red"),
    labels = c("Absence", "Presence")
  ) +
  facet_wrap(~ time, nrow = 2) +
  labs(
    x = "x", 
    y = "y",
    fill = NULL,   # ← elimina título de leyenda de fill
    color = NULL   # ← elimina título de leyenda de color
  ) +
  geom_tile(
    data = df_presencias,
    aes(x = x, y = y, color = "Sampled presence"),  # ← mapea a categoría
    inherit.aes = FALSE,
    width  = 1,
    height = 1,
    fill   = NA,
    linewidth = 0.35
  ) +
  geom_tile(
    data = df_ausencias,
    aes(x = x, y = y, color = "Sampled absence"),   # ← mapea a categoría
    inherit.aes = FALSE,
    width  = 1,
    height = 1,
    fill   = NA,
    linewidth = 0.35
  ) +
  scale_color_manual(
    values = c("Sampled presence" = "purple",
               "Sampled absence"  = "black")
  ) +
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme(
    axis.text = element_text(size = 13), 
    axis.title = element_text(size = 14),
    legend.position   = "bottom",
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.grid        = element_blank(),
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 12)
  )


(p_bathy | p_tem ) / (p_latent | p_spp) +
  plot_annotation(tag_levels = "A") 

# save(df_simulacion_predict, file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/datos_simulacion_predict.RData")

# save(df_simulacion_fit, file = "D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/df_simulacion_fit.RData")
# 
# load("D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/df_simulacion_fit.RData")
# 
# load("D:/MarineBeacon/Rcodes/Exploration_analysis/GAM_INLA/datos_simulacion/datos_simulacion_predict.RData")

# MUESTREO ESPECIE ALEATORIA EN UNA ZONA DETERMINADA ----------------------
# 
# # pesos para las covariables y el latent
# n_time        <- 12
# n_pts         <- 3000     #3000 puntos de muestreo
# beta_latent <- 4.5
# beta_temp   <- 1.5
# beta_bathy  <- 0.08
# frac_ocup_media <- 0.20  
# frac_ocup <- 0.20  
# sd_frac_ocup    <- 0.01     # variabilidad (puedes ajustar el valor)
# 
# df_grid <- grid_temp_df %>%
#   inner_join(latent_time_df, by = c("x","y","time")) %>%
#   inner_join(grid_bathy_df,     by = c("x","y","time"))
# 
# 
# df_list <- vector("list", n_time)
# 
# with_seed(123456789, {
# 
# for (t in seq_len(n_time)) {
#   
#   grid_t <- df_grid %>% filter(time == t)
#   
#   frac_ocup_t <- min(max(rnorm(1, mean = frac_ocup_media, sd = sd_frac_ocup), 0.10), 0.30)
#   
#   # Calcula el score de idoneidad
#   grid_t <- grid_t %>%
#     mutate(score = beta_latent * latent +
#              beta_temp   * temp +
#              beta_bathy  * bathy)
#   
#   # Determina el umbral para ocupar exactamente el % deseado
#   umbral <- quantile(grid_t$score, probs = 1 - frac_ocup)
#   
#   grid_t <- grid_t %>%
#     mutate(
#       pres = ifelse(score > umbral, 1, 0)  # top 20% son presencias
#     )
#   
#   df_list[[t]] <- grid_t
# }
# 
# df_ocup <- bind_rows(df_list)
# 
# })
# 
# #Cuadrado superior izquierda
# side_len <- 60
# 
# df_ocup <- df_ocup %>% 
#   mutate(
#     accessible = (x >= 100 - side_len) &     # “derecha”:  x de 60 a 100
#       (y >= 100 - side_len)       # “arriba”:   y de 60 a 100
#   )
# 
# rect_df <- data.frame(
#   xmin = 100 - side_len,
#   xmax = 100,
#   ymin = 100 - side_len,
#   ymax = 100
# )
# 
# ggplot(df_ocup, aes(x = x, y = y, fill = factor(pres))) +
#   geom_tile(color = NA) +
#   geom_rect(data = rect_df,
#             aes(xmin = xmin, xmax = xmax,
#                 ymin = ymin, ymax = ymax),
#             inherit.aes = FALSE,
#             colour = "blue", fill = NA, linewidth = 0.6) +
#   coord_equal() +
#   scale_fill_manual(values = c("0" = "grey90", "1" = "red"),
#                     name   = "Presencia",
#                     labels = c("Ausencia", "Presencia")) +
#   facet_wrap(~ time, ncol = 4) +
#   labs(title = "Mapa de ocupación espacial por unidad de tiempo\n(con área de muestreo en azul)",
#        x = "X", y = "Y") +
#   theme_minimal()
# 
# target_ratio_bounds <- c(0.01, 0.025)   #El mínimo que vamos a simular 
# n_time  <- 12
# n_pts   <- 3000
# ratio_vec <- runif(n_time,
#                    min = target_ratio_bounds[1],
#                    max = target_ratio_bounds[2])
# 
# df_muestreo_list <- vector("list", n_time)
# 
# with_seed (123456789, {
#   
# for (t in seq_len(n_time)) {
#   
#   grid_t <- df_ocup %>% filter(time == t, accessible)
#   pres_t <- grid_t  %>% filter(pres == 1)
#   abs_t  <- grid_t  %>% filter(pres == 0)
#   
#   ratio_t <- ratio_vec[t]                
#   
#   # 1) ¿cuántas ausencias permite el ratio si usáramos TODAS las presencias?
#   max_abs_by_pres <- floor(nrow(pres_t) / ratio_t)
#   # 2) ¿cuántas presencias permite el ratio si usáramos TODAS las ausencias?
#   max_pres_by_abs <- floor(nrow(abs_t)  * ratio_t)
#   
#   #   (1+ratio_t) * k  no debe superar n_pts ni la disponibilidad
#   k <- min(
#     floor(n_pts / (1 + ratio_t)),
#     max_abs_by_pres,
#     floor(max_pres_by_abs / ratio_t)
#   )
#   
#   n_abs  <- k
#   n_pres <- floor(k * ratio_t)
#   
#   muestra_abs  <- abs_t  %>% sample_n(n_abs)
#   muestra_pres <- pres_t %>% sample_n(n_pres)
#   
#   df_muestreo_list[[t]] <- bind_rows(muestra_abs, muestra_pres)
# }
# 
# 
# df_muestra <- bind_rows(df_muestreo_list)
# })
# 
# table(df_muestra$pres) #ratio
# 
# mean(df_muestra$pres) #prevalencia
# 
# 
# df_muestra %>%
#   group_by(time) %>%
#   summarise(
#     n_pres = sum(pres == 1),
#     n_abs  = sum(pres == 0),
#     ratio  = n_pres / n_abs,
#     total_sample = n_pres + n_abs,
#     prevalence = n_pres/total_sample
#   )
# 
# df_simulacion_fit <- df_muestra
# 



# MUESTREO ALEATORIO DE UNA ESPECIE CON PROB DE DETECCION -----------------
# n_time        <- 12
# n_pts         <- 3000     #3000 puntos de muestreo
# target_prev   <- 0.20     # obligo a que tenga una prevalencia de menos de 0.2
# occ_frac_area <- 0.3     # área ocupada (especie rara)
# p_det         <- target_prev / occ_frac_area #probabilidad de detectar la especie en el area ocupada
# #es una probabilidad de deteccion condicional para que se obtenga el target de prevalencia que buscamos
# 
# # Efectos de cada componente en el logit
# beta_temp   <- 1.5        # coeficiente de la temperatura
# beta_latent <- 4.5        # coeficiente del campo latente (más fuerte)
# beta_bathy <- -0.08 #preferencia a fondos someros
# 
# #funcion logistica inversa para convertir escala logit a probabilidad entre 0 y 1
# #es la base de los modelos logisticos
# inv_logit <- function(x) 1/(1 + exp(-x))
# 
# #Unimos todo en un solo df
# 
# df_grid <- grid_temp_df %>%
#   inner_join(latent_time_df, by = c("x","y","time")) %>%
#   inner_join(grid_bathy_df,     by = c("x","y","time"))
# 
# # Ahora df_grid tiene x, y, temp, latent, bathy, time
# head(df_grid)
# summary(df_grid)
# 
# #lista 
# sim_list <- vector("list", n_time)
# 
# for (t in seq_len(n_time)) {
#   grid_t <- df_grid %>% filter(time == t) #Filtramos los datos de temp  y latent estructura de ese tiempo en el que el bucle esta operando
#   
#   # 40% del area ocupada SIN haber  muestreado, de lo que se parte
#   target_occ <- occ_frac_area
#   
#   #Medias de temperatura y de la estructura espacial aleatoria
#   #para que de media nos de la fraccion ocupada 
#   mean_T <- mean(grid_t$temp)
#   mean_L <- mean(grid_t$latent)
#   mean_B <- mean(grid_t$bathy)
#   
#   #Intercepto = logit(fraccion ocupada) - Beta tem * media temperatura - beta latente * media latente
#   #El intercepto compensa los valores para que la fracción ocupada sea 0.4
#   intercept_t <- qlogis(target_occ) -
#     beta_bathy * mean_B -
#     beta_temp   * mean_T -
#     beta_latent * mean_L 
#   
#   # Muestreo homogéneo
#   samp_idx  <- sample(nrow(grid_t), n_pts, replace = TRUE) #Simulamos el muestreo, n_pts aleatorios 
#   samp_data <- grid_t[samp_idx, ]
#   
#   #Prediccion lineal, core de cada punto 
#   
#   linpred <- intercept_t + 
#     beta_bathy   * samp_data$bathy +
#     beta_temp   * samp_data$temp +
#     beta_latent * samp_data$latent
#   
#   p_occ   <- inv_logit(linpred) #lo trasformamos en probabilidades
#   
#   #puntos binarios según la probabilidad dada = probabilidad de que la especie ocupe el sitio * probabilidad de que estando presente sea detectad
#   y_t <- rbinom(n_pts, size = 1, prob = p_occ * p_det)
#   #con la probabilidad de deteccion hacemo que algunos  sean ceros 
#   
#   sim_list[[t]] <- data.frame(
#     x      = samp_data$x,
#     y      = samp_data$y,
#     temp   = samp_data$temp,
#     bathy   = samp_data$bathy,
#     pres   = y_t,
#     time   = t
#     
#   )
# }
# 
# df_sampling <- bind_rows(sim_list)
# 
# df_sampling$pres_f <- factor(df_sampling$pres,
#                              levels = c(0, 1),
#                              labels = c("Ausencia", "Presencia"))
# 
# windows();ggplot(df_sampling, aes(x = x, y = y, color = pres_f)) +
#   geom_point(size = 1.2, alpha = 0.8) +
#   facet_wrap(~ time, ncol = 4) +
#   coord_equal() +
#   scale_color_manual(
#     values = c("Ausencia" = "lightgray", "Presencia" = "red")
#   ) +
#   labs(
#     title = "Presencia/Ausencia muestreada por unidad de tiempo",
#     x     = "Coordenada X",
#     y     = "Coordenada Y",
#     color = NULL
#   ) +
#   theme_minimal()
# 
# 
# pres_df <- df_sampling %>% filter(pres == 1)
# 
# # Plot combinado
# windows();ggplot(latent_time_df, aes(x = x, y = y, fill = latent)) +
#   geom_tile(color = NA) +
#   geom_point(
#     data        = pres_df,
#     aes(x = x, y = y),
#     inherit.aes = FALSE,
#     color       = "red",        # Color para los puntos de presencia
#     size        = 1.2,
#     alpha       = 0.7
#   ) +
#   facet_wrap(~ time, ncol = 4) +
#   coord_equal() +
#   scale_fill_viridis_c(option = "D") +   # <- Paleta continua para el campo latente
#   labs(
#     title = "Ocupación espacial con presencias muestreadas superpuestas",
#     x     = "Coordenada X",
#     y     = "Coordenada Y",
#     fill  = "Latente"
#   ) +
#   theme_minimal()
# 
# 
# windows();ggplot(grid_temp_df, aes(x = x, y = y, fill = temp)) +
#   geom_tile(color = NA) +
#   geom_point(
#     data        = pres_df,
#     aes(x = x, y = y),
#     inherit.aes = FALSE,
#     color       = "red",        # Color para los puntos de presencia
#     size        = 1.2,
#     alpha       = 0.7
#   ) +
#   facet_wrap(~ time, ncol = 4) +
#   coord_equal() +
#   scale_fill_viridis_c(option = "D") +   # <- Paleta continua para el campo latente
#   labs(
#     title = "Ocupación espacial con presencias muestreadas superpuestas",
#     x     = "Coordenada X",
#     y     = "Coordenada Y",
#     fill  = "Temperatura"
#   ) +
#   theme_minimal()
# 
# 
# windows();ggplot(grid_bathy_df, aes(x = x, y = y, fill = bathy)) +
#   geom_tile(color = NA) +
#   geom_point(
#     data        = pres_df,
#     aes(x = x, y = y),
#     inherit.aes = FALSE,
#     color       = "red",        # Color para los puntos de presencia
#     size        = 1.2,
#     alpha       = 0.7
#   ) +
#   facet_wrap(~ time, ncol = 4) +
#   coord_equal() +
#   scale_fill_viridis_c(option = "D") +   # <- Paleta continua para el campo latente
#   labs(
#     title = "Ocupación espacial con presencias muestreadas superpuestas",
#     x     = "Coordenada X",
#     y     = "Coordenada Y",
#     fill  = "Temperatura"
#   ) +
#   theme_minimal()
# 


# PLOTS  ------------------------------------------------------------------

df_simulacion_fit$pres_f <- factor(df_simulacion_fit$pres,
                            levels = c(0, 1),
                            labels = c("Ausencia", "Presencia"))

windows();ggplot(df_simulacion_fit, aes(x = x, y = y, color = pres_f)) +
  geom_point(size = 1.2, alpha = 0.8) +
  facet_wrap(~ time, ncol = 4) +
  coord_equal() +
  scale_color_manual(
    values = c("Ausencia" = "lightgray", "Presencia" = "red")
  ) +
  labs(
    title = "Presencia/Ausencia muestreada por unidad de tiempo",
    x     = "Coordenada X",
    y     = "Coordenada Y",
    color = NULL
  ) +
  theme_minimal()


pres_df <- df_simulacion_fit %>% filter(pres == 1)

# Plot combinado
windows();ggplot(df_simulacion_predict, aes(x = x, y = y, fill = latent)) +
  geom_tile(color = NA) +
  geom_point(
    data        = pres_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color       = "red",        # Color para los puntos de presencia
    size        = 1.2,
    alpha       = 0.7
  ) +
  facet_wrap(~ time, ncol = 4) +
  coord_equal() +
  scale_fill_viridis_c(option = "D") +   # <- Paleta continua para el campo latente
  labs(
    title = "Ocupación espacial con presencias muestreadas superpuestas",
    x     = "Coordenada X",
    y     = "Coordenada Y",
    fill  = "Latente"
  ) +
  theme_minimal()


windows();ggplot(df_simulacion_predict, aes(x = x, y = y, fill = temp)) +
  geom_tile(color = NA) +
  geom_point(
    data        = pres_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color       = "red",        # Color para los puntos de presencia
    size        = 1.2,
    alpha       = 0.7
  ) +
  facet_wrap(~ time, ncol = 4) +
  coord_equal() +
  scale_fill_viridis_c(option = "D") +   # <- Paleta continua para el campo latente
  labs(
    title = "Ocupación espacial con presencias muestreadas superpuestas",
    x     = "Coordenada X",
    y     = "Coordenada Y",
    fill  = "Temperatura"
  ) +
  theme_minimal()


windows();ggplot(df_simulacion_predict, aes(x = x, y = y, fill = bathy)) +
  geom_tile(color = NA) +
  geom_point(
    data        = pres_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color       = "red",        # Color para los puntos de presencia
    size        = 1.2,
    alpha       = 0.7
  ) +
  facet_wrap(~ time, ncol = 4) +
  coord_equal() +
  scale_fill_viridis_c(option = "G", direction = -1) +   # <- Paleta continua para el campo latente
  labs(
    title = "Ocupación espacial con presencias muestreadas superpuestas",
    x     = "Coordenada X",
    y     = "Coordenada Y",
    fill  = "Temperatura"
  ) +
  theme_minimal()



ggplot(df_simulacion_predict, aes(x = x, y = y, fill = factor(pres))) +
  geom_tile(color = NA) +
  # geom_point(
  #   data        = pres_df,
  #   aes(x = x, y = y),
  #   inherit.aes = FALSE,
  #   color       = "black",        # Color para los puntos de presencia
  #   size        = 1.2,
  #   alpha       = 0.7
  # ) +
  coord_equal() +
  scale_fill_manual(values = c("0" = "grey90", "1" = "red"),
                    name = "Presencia",
                    labels = c("Ausencia", "Presencia")) +
  facet_wrap(~ time, ncol = 4) +
  labs(title = "Mapa de ocupación espacial por unidad de tiempo",
       x = "X", y = "Y") +
  theme_minimal()


abs_df <- df_simulacion_fit %>% filter(pres == 0)

ggplot(df_simulacion_predict, aes(x = x, y = y, fill = factor(pres))) +
  geom_tile(color = NA) +
  geom_point(
    data        = abs_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color       = "black",        # Color para los puntos de presencia
    size        = 1.2,
    alpha       = 0.7
  ) +
  coord_equal() +
  scale_fill_manual(values = c("0" = "grey90", "1" = "red"),
                    name = "Presencia",
                    labels = c("Ausencia", "Presencia")) +
  facet_wrap(~ time, ncol = 4) +
  labs(title = "Mapa de ocupación espacial por unidad de tiempo",
       x = "X", y = "Y") +
  theme_minimal()