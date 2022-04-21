## Figure 4 (b) Intracellular nucleoside triphosphate (TP) concentration following in vitro incubation with the parent drug remdesivir
**Author:** Mohammadreza Samieegohar

## Execution
The default is to generate  Figures 4b Optimum curve.
Rscript Main.R

For generating the figure 4b , Optimal and uncertainty band (2000 virtual population) , the following option needs to use: 
Rscript Main.R -c Population

For performing a faster run, the virtual population (for example, 113) can be adjusted by: 
Rscript Main.R -c Population -p 113

or the options can be adjusted in **Main.R** as well in the following way:

_parser<-add_option(parser, c("-c", "--case"), default="Optimum",type="character", help="Population or Optimum")_

_parser<-add_option(parser, c("-p", "--popSize"), default="2000", help="population size")_

## Folders and Files
**Main.R** is the main code for generating the figure

**results** includes the output figures

**inputs_parameters** includes the input data

**funs** includes the helper functions

**models** includes the model parameters and states

## Software
All simulations are run in R version 4.0.2.

This code uses the following R packages: ggplot2, gridExtra, optparse, deSolve.
