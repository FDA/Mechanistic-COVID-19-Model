# This code generates "Figure 4 (a) Calibrating of the PK model." 

## Running the code---------------------------------------------------------
All simulations are run in R-4.0.2 
This code uses the following R packages:
library(ggplot2)
library(gridExtra)
library(optparse)
library(deSolve)

## Folders and Files-------------------------------------------------------
"data"        			includes the experimental data
"results"     			includes the output figures
"input_parameters" 		includes the input data
"funs" 					includes the helper functions
"Main.R" 				is the main code for generating the figure 4a

## Running-----------------------------------------------------------------
The default is to generate  Figures 4a by using a 2000 virtual population 
Rscript Main.R
for generating figure by different virtual population (for example virtual population=113), following option is used:

Rscript Main.R -p 113

Or
change the default value from "2000" to "113" in Main.R
parser<-add_option(parser, c("-p", "--popNum"), default="113", help="population number")
####-----------------------------------------------------
The "readme" file and codes are written by:
Mohammadreza Samieegohar