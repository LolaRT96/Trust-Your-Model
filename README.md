# "TRUST YOUR MODEL: UNBALANCED DATA FROM RARE SPECIES DO NOT IMPLY UNRELIABLE PREDICTIONS"

This repository contains the R code used in the simulation case study for our manuscript submitted to [insert journal].
The objective of the study is to evaluate concerns regarding resampling procedures and the potential loss of ecological realism when modelling rare species. To do so, we generate different prevalence scenarios through an undersampling strategy and fit them into two modelling approaches (Random Forest and GAM–INLA) using 100 bootstrap iterations.

Please note that exploratory analyses and the Python scripts used for data visualisation (figures and maps) are not included in this repository.

#Pipeline 
The workflow to reproduce the simulated data and fit the models is organised as follows:

1. Generate the simulated data (01_SimulationData.R) 
2. Fit the Random Forest models with the generated data (02_SimulationRF.R)
3. Fit the GAM-INLA models with the generated data (02_SimulationGAM.R)

#Data Availability
The real case study data and associated scripts cannot be shared due to data availability restrictions and usage policies.
However, the scripts used for the empirical case study can be provided upon request.

#References:
a) Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA. Elias T. Krainski, Virgilio Gómez-Rubio, Haakon Bakka, Amanda Lenzi, Daniela Castro-Camilo, Daniel Simpson, Finn Lindgren and Håvard Rue. CRC Press/Taylor and Francis Group, 2019. (https://becarioprecario.bitbucket.io/spde-gitbook/ch-stapp.html#sec:hgst)

b) Breiman, L., 2001. Random Forests.

