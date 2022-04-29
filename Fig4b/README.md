## This code generates "Figure 4b Intracellular nucleoside triphosphate (TP) concentration following in vitro incubation with the parent drug remdesivir" 

## Running the code---------------------------------------------------------
All simulations are run in R-4.0.2   
This code uses the following R packages:
<pre>
ggplot2   [version 3.3.3]  
gridExtra [version 2.3] 
optparse  [version 1.6.6] 
deSolve   [version 1.28]  
</pre>
## Folders and Files-------------------------------------------------------
<pre>
"Fig4b"                 includes the output figures  
"Inputs_parameters"     includes the input data  
"Main.R"                is the main code for generating the figure 4b  
</pre>
## Running simulations-----------------------------------------------------
To generate Figure 2, Just Optimum case (the single best-fit patient):    
  The default is to generate  Figures 4b Optimum curve. 
  <pre>
    Under Linux, one simply run:    
    Rscript Main.R
     
   Under other operation systems (windows, Mac, etc.) after launching a R console:    
    source("Main.R")
 </pre>  
To generating the figure 4b, Optimal and uncertainty band (2000 virtual population):
   <pre>
    Under Linux: 
    Rscript Main.R -c Population  

    Under other operation systems (windows, Mac, etc.), one needs to change the default value from "Optimum" to "Population" in line 12 of Main.R below, then     execute source("Main.R") in a console:    
    parser<-add_option(parser, c("-c", "--case"), default="Optimum",type="character", help="Population or Optimum")
  </pre>
 
####-----------------------------------------------------  
The "readme" file and codes are written by:   
Mohammadreza Samieegohar
