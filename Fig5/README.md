## This set of code generates Figure 5 (calibrating and validating the full model based on the clinical endpoint of time to recovery) in the manuscript.
It simulates the clinical trial of ACTT. 100 subpopulations (each subpopulation 100 virtual patients) were sampled from the 2000 virtual patients (used in Fig 2) to match the recovered% on each day in the placebo arm as reported by ACTT (model calibraiton). Then the same subpopulations went through remdesivir therapy simulation to predict the recovered% on each day, which is subsequently compared to the ACTT-reported values in the remdesivir arm (model validation). Similar to ACTT, patients were stratified into mild and severe strata in Figure 5.

## Running the code---------------------------------------------------------
All simulations are run in R-4.0.2   
This code uses the following R packages:
<pre>
ggplot2   [version 3.3.3]  
gridExtra [version 2.3]    
grid      [version 4.0.2]  
optparse  [version 1.6.6]  
deSolve   [version 1.28]  
</pre>
## Folders and Files-------------------------------------------------------
<pre>
"Data_Covid"            includes the CTT clinical data.  
"Fig5"                  includes the output figures.  
"Inputs_parameters"     includes the input parameters and sub-population IDs.  
"Main.R" 		is the main code for generating the figures 5 a,b,c,d.  
"Calibration.R" 	is a subroutine that calculate the calibration (placebo)    data.  
"Validation.R" 		is a subroutine that calculate the validation  (remdesivir) data.  
"eventRMD.R"            is a subroutine that used in "Validation.R" for remdesivir injection scheme.  
</pre>

## Description of the model -----------------------------------------------------------------

Two different formats of the COVID-19 full model are provided for running the code to generate Figure 5:  

I)  There is a  *so  version of the model that is 		      compiled for linux, but Faster to run.   
PAHT: Model/delaymodelfun.so  

II) There is an *R   version that is a pure R implementation of the model. It can be run across platforms (operatign systems), but slower to run.  
PAHT: Model/delaymymod.R  

Note: If you use non-linux operating system (e.g., windows, Mac, or a version of linux that cannot use the provided *so file directly), but you want to take advantage of the faster running for option I, you can compile the provided c source code for your own system. For details please see https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/SHLIB  
Some modifications to Main.R will then be needed to tell the script to use the compiled model file you just created.  

Also note that the figures in our manuscript were created by running the *so version of the model compiled on CentOS Linux Version 7. If you run the *so version of the model on other versions of linux, or run the *R version of the model on any other operation systems, the results may be slightly different due to machine-to-machine and platform-to-platform differences.

## Running simulations -----------------------------------------------------------------
*To generate Figure 5, 
The default is to generate plots of mild patient by using an R version of the model:
<pre>
Under Linux, one simply run: 

   In mild cases   (figures a,c) by using an *R version of the model:  
		Rscript Main.R
		
   In severe cases (figures b,d) by using an *R version of the model:  
		Rscript Main.R -d severe 
</pre>	
Under other operation systems (windows, Mac, etc.), one needs to change the default value from "mild" to "severe" in line 15 of Main.R below, then execute source("Main.R") in a console: 
<pre>	
   parser<-add_option(parser, c("-d", "--severity"), default="severe",type="character", help="severity group,options:'mild' or 'severe'")  
</pre>		
Of note, for any of the runs above, to use the compiled model (*.so version of the model) for faster running under linux, one should use the command line option -m. For example:
<pre>
   In severe cases (figures b,d) by using an *so version of the model:  
		Rscript Main.R -d severe -m C  

   In mild cases   (figures a,c) by using an *so version of the model:  
		Rscript Main.R -m C  
</pre>
####-----------------------------------------------------   
The "readme" file and codes are written by:  
Mohammadreza Samieegohar

