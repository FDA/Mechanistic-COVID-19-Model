# This code generates "Figure 2 Disease model calibration." 

## Running the code---------------------------------------------------------
All simulations are run in R-4.0.2 
This code uses the following R packages:
library(ggplot2)
library(ggplot2)
library(gridExtra)
library(grid)
library(optparse)
library(deSolve)

## Folders and Files-------------------------------------------------------
"Data_Covid"        	includes the experimental data
"Data"       	 		includes the experimental data
"Figs"    	 			includes the output figures
"Inputs_parameters" 	includes the input data
"funs" 					includes the helper functions
"Main.R" 				is the main code for generating the figure 4A

## Running-----------------------------------------------------------------
Two different formats are provided for running the code to generate Figure 2:
II) There is an *R   version of the model that is Windows and Linux compatible, but slower to run.
PAHT: Data/delaymodelfun.R
I)  There is a  *so  version of the model that is 		      Linux compatible, but Faster to run.
PAHT: Data/delaymymod.so

Note: the *so format of model can be compiled/generated by 
cd Data
R CMD SHLIB delaymymod.c
###---------------------------------------------------------------------------
The default is to generate Figure 2, Optimum curve in mild cases by using an R version of the model:
Rscript Main.R

To generate Figure 2, Optimum curve in severe cases by using an *R version of the model:
Rscript Main.R -d severe

To generate Figure 2, Optimum curve in mild cases by using an *so version of the model:
Rscript Main.R -m C

To generate Figure 2, Optimum curve in severe cases by using an *so version of the model:
Rscript Main.R -d severe -m C

To generate Figure 2, Optimum curve + uncertainty band (2000 virtual patients) in severe cases by using an *R version of the model:
Rscript Main.R -d severe -c Population

To generate Figure 2, Optimum curve + uncertainty band (2000 virtual patients) in mild cases by using an *R version of the model:
Rscript Main.R -c Population

To generate Figure 2, Optimum curve + uncertainty band (2000 virtual patients) in severe cases by using an *so version of the model:
Rscript Main.R -d severe -c Population -m C


To generate Figure 2, Optimum curve + uncertainty band (2000 virtual patients) in mild cases by using an *so version of the model:
Rscript Main.R -c Population -m C


To generate Figure 2 with a optional population (for example 113) in severe cases and *so version of model:
Rscript Main.R -d severe -c Population -m C -p 113

####-----------------------------------------------------
The "readme" file and codes are written by:
Mohammadreza Samieegohar