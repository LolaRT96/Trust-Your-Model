 
#Authors: Maria Dolores Riesgo (IEO-CSIC)
#Research paper: "Trust your model"
#General objective: generate simulated data for fit the models 

#R version 4.4.2

#If necessary 

rm(list=ls(all=TRUE)) 

#library

library(sp)  # spatial data classes
library(geoR)
library(spdep)  # spatial dependence / autocorrelation
library(fields) # spatial interpolation, kriging and plots
library(gridExtra) # arrange multiple ggplots
library(ggplot2) # data visualization
library(openxlsx)  # read/write Excel files (.xlsx)
library(INLA) # Bayesian spatial/temporal modeling 
library(GGally) # extensions for ggplot2 (e.g., ggpairs)
library(tidyverse) # data manipulation & visualization suite
library(gstat)  # geostatistical modeling & kriging
library(pROC) # ROC curves & AUC metrics
library(lobstr)  # memory usage inspection
library(fmesher) # mesh construction for INLA/SPDE
library(withr) # temporary environment changes
library(viridis) # colorblind-friendly color scales

#Set the seed

set.seed(123456789)
print(.Random.seed[1:4])

# Study area  -------------------------------------------------------

campo <- function(c1, c2, c3, c4) {
  xy <- expand.grid(
    seq(c1, c2, length.out = 100),
    seq(c3, c4, length.out = 100)
  )
  cbind(xy[, 1], xy[, 2])
}

c1 <- 0; c2 <- 100; c3 <- 0; c4 <- 100 

loc_xy <- campo(c1, c2, c3, c4) 


# Mesh  ----------------------

with_seed(123456789, {mesh2d_sim <- fm_mesh_2d_inla(loc.domain = loc_xy, 
                                                    max.edge = c(8, 15))  # Creamos la malla 2D
})

windows();plot(mesh2d_sim, main = "Malla", asp = 1, lwd = 0.5)

mesh2d_sim$n

# Predictors -----------------------------------------------------

#fixed variable (Bathymetry)

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
  theme( 
    axis.text = element_text(size = 13), axis.title = element_text(size = 14),
    legend.position   = "right",
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.grid        = element_blank(),
    strip.background  = element_blank(),                        
    strip.text        = element_text(face = "bold", size = 12)  
  )


#Temperature (vary with unit times)

n_time <- 12
temp_list <- vector("list", n_time)

#Sinusoidal function with an average temperature of 20 °C

base_temp <- function(x, y) {
  20 + 3 * sin(pi * x / 100) + 2 * cos(pi * y / 100)
}


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

#changes by time 

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


# Spatial structure  --------------------------------


with_seed(123456789, {

#Precision
prec <- 1 / 4
#Range and variance 
rho <- 20
sigma <- sqrt(5)
phi <- 0.9  #parameter controlling the temporal dependency between two time units
#close to 1: high persistence, low noise, less local variation
#close to 0: patterns redefined at each time step
#close to -1 or negative: high areas in t1 are low in t2 and so on. Pendulum effect. 
#ARIMA PROCESS: ut = phi*ut-1 + errort, ∣ϕ∣<1

sigma_eps <- sigma * sqrt(1 - phi^2) #Noise deviation: each time unit we add noise of a specific amplitude. 

u0_nodes <- fmesher::fm_matern_sample(mesh2d_sim, n = 1, rho = rho, sigma = sigma) #Guassian Random Field

u0_nodes <- u0_nodes - mean(u0_nodes) #center at zero to eliminate unwanted deviations from the field. 

n_time <- 12 
A_grid <- fm_basis(mesh2d_sim, loc_xy) #matrix that associates each node of the mesh with each point of the grid

latent_list <-  vector("list", n_time)

u_prev <- u0_nodes #node starting point

for(t in seq(n_time)) {
  
  eps_t <-  fm_matern_sample(mesh2d_sim, n=1, rho = rho, sigma = sigma) #noise
  eps_t <- eps_t - mean(eps_t) #mean to zero
  
  u_t <- phi * u_prev + eps_t
  u_prev <- u_t #result at the nodes for time t
  
  latent_t <- drop(A_grid %*% u_t) #we interpolate those values u_t at each grid point 
  
  latent_list[[t]] <- data_frame(
    x =loc_xy[,1],
    y=loc_xy[,2],
    latent = latent_t,
    time = t
  )
  
}

latent_time_df <- bind_rows(latent_list)

} )

p_latent <- ggplot(latent_time_df, aes(x=x, y=y, fill=latent)) +
  geom_tile()+
  facet_wrap(~time,  nrow = 3) +
  coord_equal()+
  scale_fill_viridis_c() +
  labs(
    x = "x",
    y = "y",
    fill = "Spatial structure"  
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
         strip.background  = element_blank(),                        
         strip.text        = element_text(face = "bold", size = 12)  
  )

(p_bathy | p_tem ) / (p_latent)

#AREA OCCUPIED BY THE SIMULATED SPECIES  ---------------------------------------------

# weights of the predictors
n_time        <- 12
beta_latent <- 4.5
beta_temp   <- 1.5
beta_bathy  <- 0.08
frac_ocup_media <- 0.20  #proportion of the average occupied area (20%)
frac_ocup <- 0.20  #proportion of the occupied area (20%)
sd_frac_ocup    <- 0.01 #variability of the occupied area

df_grid <- grid_temp_df %>%
  inner_join(latent_time_df, by = c("x","y","time")) %>%
  inner_join(grid_bathy_df,     by = c("x","y","time"))

df_list <- vector("list", n_time)

with_seed(123456789, {

for (t in seq_len(n_time)) {
  
  grid_t <- df_grid %>% filter(time == t)
  
  frac_ocup_t <- min(max(rnorm(1, mean = frac_ocup_media, sd = sd_frac_ocup), 0.10), 0.30)
  
  #Suitability score
  grid_t <- grid_t %>%
    mutate(score = beta_latent * latent +
             beta_temp   * temp +
             beta_bathy  * bathy)
  
  #  Determine the threshold to occupy exactly the desired %
  umbral <- quantile(grid_t$score, probs = 1 - frac_ocup)
  
  grid_t <- grid_t %>%
    mutate(
      pres = ifelse(score > umbral, 1, 0)  # top 20% son presencias
    )
  
  df_list[[t]] <- grid_t
}

df_ocup <- bind_rows(df_list)
})


# Random sapling of the simulated species ---------------------------------

target_ratio_bounds <- c(0.01, 0.025) #desired prevalences
n_time              <- 12 
n_pts               <- 3000 #points of sampling

# random prevalence for each t in n_time
ratio_vec <- runif(n_time, min = target_ratio_bounds[1], max = target_ratio_bounds[2])

df_muestreo_list <-  vector("list", n_time)

with_seed(123456789, {
for (t in seq_len(n_time)) {
  
  grid_t <- df_ocup %>% filter(time == t)
  pres_t <- df_ocup %>% filter(pres == 1)
  abs_t <- df_ocup %>% filter(pres == 0)
  
  ratio_t <- ratio_vec[t]
  
  #absences attendance of 3,000 points to meet the set ratio
  n_abs  <- floor(n_pts / (1 + ratio_t))
  n_pres <- n_pts - n_abs
  
  n_abs  <- min(n_abs,  nrow(abs_t))
  n_pres <- min(n_pres, nrow(pres_t))
  
  muestra_abs  <- abs_t  %>% sample_n(n_abs)
  muestra_pres <- pres_t %>% sample_n(n_pres)
  
  muestra_t <- bind_rows(muestra_abs, muestra_pres) 
  
  df_muestreo_list[[t]] <- muestra_t 
  
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



# Plotting results --------------------------------------------------------

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
  scale_x_continuous(expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme( axis.text = element_text(size = 13), axis.title = element_text(size = 14),
         legend.position   = "bottom",
         axis.line         = element_line(color = "black", linewidth = 0.4),
         axis.ticks        = element_line(color = "black", linewidth = 0.3),
         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
         panel.grid        = element_blank(),
         strip.background  = element_blank(),                       
         strip.text        = element_text(face = "bold", size = 12)  
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
    fill = NULL,   
    color = NULL   
  ) +
  geom_tile(
    data = df_presencias,
    aes(x = x, y = y, color = "Sampled presence"),  
    inherit.aes = FALSE,
    width  = 1,
    height = 1,
    fill   = NA,
    linewidth = 0.35
  ) +
  geom_tile(
    data = df_ausencias,
    aes(x = x, y = y, color = "Sampled absence"),   
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

