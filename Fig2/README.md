## This set of code generates Figure 2 (Disease model calibration) in the manuscript-------
It runs simulations for a population of 2000 virtual patients to produce time courses of various COVID-19 disease biomarkers (viral titer, cytokines, antibodies, etc.) and overlay them with clinical data. Two populations (mild and severe) of virtual patients are provided. In addition to simulating the whole population, an option is also given in the code to simulate a single best-fit (optimum) patient, which is much faster.

## Requirement for running the code---------------------------------------------------------
All simulations are run in R-4.0.2   
This code uses the following R packages:  
ggplot2   [version 3.3.3]  
gridExtra [version 2.3]  
grid      [version 4.0.2]  
optparse  [version 1.6.6]   
deSolve   [version 1.28] 

## Description of Folders and Files-------------------------------------------------------
"Data_Covid"        	includes the experimental data  
"Fig2"    	 			includes the output figures  
"Inputs_parameters" 	includes the input data  
"Model" 				includes the models in different formats  
"Main.R" 				is the main code for generating the figure 2  
"Plot_Optimum R" 		is a subroutine that is used in Main.R for plotting the Optimum case.  
"Plot_Population.R"     is a subroutine that is used in Main.R for plotting the Population case.  
## Description of the model -----------------------------------------------------------------
Two different formats of the COVID-19 disease model are provided for running the code to generate Figure 2:  
I)  There is a  *so  version of the model that is 		      compiled for linux, but Faster to run.   
PAHT: Model/delaymodelfun.so  

II) There is an *R   version that is a pure R implementation of the model. It can be run across platforms (operatign systems), but slower to run.  
PAHT: Model/delaymymod.R  

Note: If you use non-linux operating system (e.g., windows, Mac, or a version of linux that cannot use the provided *so file directly), but you want to take advantage of the faster running for option I, you can compile the provided c source code for your own system. For details please see https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/SHLIB  
Some modifications to Main.R will then be needed to tell the script to use the compiled model file you just created.

## Running simulations-----------------------------------------------------------------
*To generate Figure 2, Just Optimum case (the single best-fit patient):    
 
   1)The default is to generate Optimum case for mild patient by using an R version of the model:  
   
   Under Linux, one simply run:    
		Rscript Main.R  
		
   Under other operation systems (windows, Mac, etc.) after launching a R console:    
                source("Main.R")

   2)To generate Optimum case for severe patient by using an *R version of the model:  
   
   Under Linux:    
		Rscript Main.R -d severe  
		
   Under other operation systems (windows, Mac, etc.), one needs to change the default value from "mild" to "severe" in line 15 of Main.R below, then execute source("Main.R") in a console:    
   parser<-add_option(parser, c("-d", "--severity"), default="severe",type="character", help="severity group,options:'mild' or 'severe'")  

    
*To generate Figure 2, Optimum curve + uncertainty band (2000 virtual patients), which is the format of Figure 2 of the manuscript:  

   1)To generate population simulaitons for mild cases by using an *R version of the model:  
   
   Under Linux:    
		Rscript Main.R -d severe -c Population  
		
   Under other operation systems (windows, Mac, etc.), one needs to change the default value from "mild" to "severe" in line 15, and from "Optimum" to "Population" in line 16, of Main.R below, then execute source("Main.R") in a console:    
   
   parser<-add_option(parser, c("-d", "--severity"), default="severe",type="character", help="severity group,options:'mild' or 'severe'")   
   parser<-add_option(parser, c("-c", "--case"), default="Population", help="options: 'Population' or 'Optimum'")  
   
   2)To generate population simulations for mild cases by using an *R version of the model:  
   
   Under Linux:    
		Rscript Main.R -c Population  
		
   Under other operation systems (windows, Mac, etc.), one needs to change the default value from "severe" to "mild" in line 15, and from "Optimum" to "Population" in line 16, of Main.R below, then execute source("Main.R") in a console:    
   parser<-add_option(parser, c("-d", "--severity"), default="mild",type="character", help="severity group,options:'mild' or 'severe'")   
   parser<-add_option(parser, c("-c", "--case"), default="Population", help="options: 'Population' or 'Optimum'")

   
   Of note, for any of the runs above, to use the compiled model (*.so version of the model) for faster running under linux, one should use the command line option -m. For example:    
		Rscript Main.R -m C  



####-----------------------------------------------------  
The "readme" file and codes are written by:    
Mohammadreza Samieegohar
