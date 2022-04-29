## This code generates "Figure 4 (a) Calibrating of the PK model." 

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
"data"        			includes the experimental data  
"Fig4a"     			includes the output figures  
"input_parameters" 		includes the input data  
"Main.R" 			is the main code for generating the figure 4a  
</pre>
## Running-----------------------------------------------------------------
The default is to generate  Figures 4a by using a 2000 virtual population
<pre>
Under Linux, one simply run:
	Rscript Main.R  
	
Under other operation systems (windows, Mac, etc.) after launching a R console:    
        source("Main.R")
</pre>
####-----------------------------------------------------  
The "readme" file and codes are written by:   
Mohammadreza Samieegohar
