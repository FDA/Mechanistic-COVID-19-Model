# This code generates "Figure 5 (a,b,c,d) calibrating disease model using the clinical placebo survival curve in mild and sever cases"

## Running the code---------------------------------------------------------
All simulations are run in R-4.0.2 
This code uses the following R packages:
library(ggplot2)
library(grid)
library(optparse)

## Folders and Files-------------------------------------------------------
"Data_Covid"        includes the experimental data
"Fig5 includes"     includes the output figures
"Inputs_parameters" includes the input data
"Main.R" 			is the main code for generating the figures 5 a,b,c,d

## Running-----------------------------------------------------------------
The default is to generate "mild" cases which are Figures 5a-b
for generating the Severe, figures(5c-d), following option is used:
Rscript Main.R -c severe
Or
change the default value from "mild" to "severe" in Main.R
parser<-add_option(parser, c("-c", "--Runcase"), default="severe",type="character", help="severe or mild")
####-----------------------------------------------------
The "readme" file and codes are written by:
Mohammadreza Samieegohar