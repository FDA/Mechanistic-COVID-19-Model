# This code generates "Figure 4b Intracellular nucleoside triphosphate (TP) concentration following in vitro incubation with the parent drug remdesivir" 

## Running the code---------------------------------------------------------
All simulations are run in R-4.0.2   
This code uses the following R packages:  
library(ggplot2   [version 3.3.3])  
library(gridExtra [version 2.3]  )  
library(optparse  [version 1.6.6])  
library(deSolve   [version 1.28] )  

## Folders and Files-------------------------------------------------------
"Fig4b"     			includes the output figures  
"Inputs_parameters" 	includes the input data  
"Main.R" 				is the main code for generating the figure 4b  

## Running-----------------------------------------------------------------
The default is to generate  Figures 4b Optimum curve.  
Rscript Main.R  

For generating the figure 4b, Optimal and uncertainty band (2000 virtual population), the following option needs to be used:  
 
Rscript Main.R -c Population  

For performing a faster run, the virtual population (for example, 113) can be adjusted by:   
Rscript Main.R -c Population -p 113  

Or  

The options can be adjusted in Main.R  
parser<-add_option(parser, c("-c", "--case"), default="Optimum",type="character", help="Population or Optimum")  
parser<-add_option(parser, c("-p", "--popNum"), default="2000", help="population number")  
####-----------------------------------------------------
The "readme" file and codes are written by:  
Mohammadreza Samieegohar
