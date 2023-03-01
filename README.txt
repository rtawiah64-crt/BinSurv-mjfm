######## Readme file ##############
#Title: Multilevel joint frailty model for hierarchically clustered survival and binary response data

The compressed file contains a readme file "README.txt", data  layout file "data_format"  and an R-code file "MultBinSurv.R".

1. data_format    # gives a tabular overview of how to prepare the data for analysis. 

3. MultBinSurv.R  # R code for the proposed multilevel joint frailty model
   
   # Notes: 
   # Set work directory: setwd("C:/Users/...")
   # Load the data into R.
   # thetau0, thetaU0, rho0, sigv0 & sigV0 are initial values for thetau^2, thetaU^2, rho, sigmav^2 & sigmaV^2, respectively. 
   # The main function is "jointBinSurv"; it calls the output of the model. 
   # Maximum iteration = 300
   # eps.reg and eps.var corresponds to epsilon1 and epsilon2, repectively: tolerance for the optimizers of regression coefficients and variance components.    
   

