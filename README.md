# "TRUST YOUR MODEL: UNBALANCED DATA FROM RARE SPECIES DO NOT IMPLY UNRELIABLE PREDICTIONS"

In this repository we present the R code from the simulation case study from our manuscript at [insert journal]. 
The focus is to answer the hypothesis about the concernings of calibration techniques under rare species contexts by using a undersampling approach to create different prevalence scenarios, fitting them into two different models (Random Forest and GAM-INLA) by using 100 bootstrap iterations. 
No exploratory analysis will be provided as well as the Phyton code for the data visualization of the graphs and maps. 
The pipeline to reproduce the data and fit into the models is:
1. Generate the simulated data (01_SimulationData.R) 
2. Fit the Random Forest models with the generated data (02_SimulationRF.R)
3. Fit the GAM-INLA models with the generated data (02_SimulationGAM.R)

The real case study data plus the scripts will be not provide due to data availability and usage policies. However, scripts for the empirical case study will be provided under request.  

References:

a) Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA. Elias T. Krainski, Virgilio Gómez-Rubio, Haakon Bakka, Amanda Lenzi, Daniela Castro-Camilo, Daniel Simpson, Finn Lindgren and Håvard Rue. CRC Press/Taylor and Francis Group, 2019. (https://becarioprecario.bitbucket.io/spde-gitbook/ch-stapp.html#sec:hgst)

b) Breiman, L., 2001. Random Forests.

