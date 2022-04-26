# This code generates "Figure 5 (a,b,c,d) calibrating disease model using the clinical placebo survival curve in mild and sever cases"

## Running the code---------------------------------------------------------
All simulations are run in R-4.0.2   
This code uses the following R packages:  
library(ggplot2   [version 3.3.3])  
library(gridExtra [version 2.3]  )  
library(grid      [version 4.0.2])  
library(optparse  [version 1.6.6])  
library(deSolve   [version 1.28] )  

## Folders and Files-------------------------------------------------------
"Data_Covid"            includes the CTT clinical data.  
"Fig5"                  includes the output figures.  
"Inputs_parameters"     includes the input parameters and sub-population IDs.  
"Main.R" 			    is the main code for generating the figures 5 a,b,c,d.  
"Calibration.R" 		is a subroutine that calculate the calibration (placebo)    data.  
"Validation.R" 			is a subroutine that calculate the validation  (remdesivir) data.  
"eventRMD.R"            is a subroutine that used in "Validation.R" for remdesivir injection scheme.  
## Running-----------------------------------------------------------------
*To generate Figure 5,  
   In severe cases (figures b,d) by using an *R version of the model:  
		Rscript Main.R -d severe -m R  

   In mild cases   (figures a,c) by using an *R version of the model:  
		Rscript Main.R -d mild -m R  

   In severe cases (figures b,d) by using an *so version of the model:  
		Rscript Main.R -d severe -m C  

   In mild cases   (figures a,c) by using an *so version of the model:  
		Rscript Main.R -d mild -m C  

####-----------------------------------------------------
The "readme" file and codes are written by:  
Mohammadreza Samieegohar
