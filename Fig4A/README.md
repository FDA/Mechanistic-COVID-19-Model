# This code generates "Figure 4 (a) Calibration of the PK model." 
**Author:** Mohammadreza Samieegohar

## Running
Executing **Main.R** generates Figure 4a by simulating a virtual population of 2000 subjects. 

<br/>
For generating the figure for a virtual population of a different size (for example, 113 subjects), following option can be used: _-p 113_
 
or, by changing the default virtual population size from _2000_ to _113_ in **Main.R** in the following way:

_parser<-add_option(parser, c("-p", "--populationSize"), default="113", help="population size")_


## Running the code
All simulations are run in R-4.0.2.

This code uses the following R packages: ggplot2, gridExtra, optparse, deSolve

## Folders and Files
"Main.R" 				is the main code for generating the figure

"results"     			includes the output figures

"input_parameters" 		includes the input data

"data"        			includes the experimental data

"funs" 					includes the helper functions


