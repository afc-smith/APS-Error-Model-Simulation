#=========================================================================================================================#
#============================================ ERROR MODEL SIMULATION =====================================================#
#======== EXAMPLE CODE FOR CASE STUDY: FAECAL CALPROTECTIN FOR DISTINGUISHING BETWEEN IBS AND IBD IN PRIMARY CARE ========#
#=========================================================================================================================#

#----------------------------- 0.0 Notes on code ---------------------------------

# This code was written by Ms. Alison F Smith to accompany the publication:
#'Smith AF, Messenger MP, Hulme CT, Turvil J, Hall PS and Shinkins B. Setting outcome-based analytical performance specifications using simulation: a case study of faecal calprotectin. Clin Chem 2019; in press'

# The code provides an illustrative example of how to run an error model simulation, using parametric or bootstrap sampling
# To avoid data sharing issues, the bootstrap method is based on a mock dataset for the YFCCP data (see notes under this method section)
# The results are therefore not intended to inform clinical practice, but rather provide a methodological framework

# DISCLAIMERS:
# Alison Smith is funded by a National Institute for Health Research (NIHR) Doctoral Research Fellowship for this research project.
# Bethany Shinkins and Mike Messenger are also supported by the NIHR Leeds In Vitro Diagnostics Co-operative. 
# This publication presents independent research. 
# The views expressed are those of the author(s) and not necessarily those of the NHS, the NIHR or the Department of Health and Social Care. 

# CODE SECTIONS:
# The code is split into code sections which can be expanded/collapsed by clicking on the arrow icon to the left of the headings
# Code sections are created by inserting "====" or "----" at the end of header lines, and can be removed by deleting this bit of code

# COMMENTS:
# All lines beginning with a "#" are comments. 
# Sections of code can be commented/uncommented by selecting the required lines of code and pressing Ctrl+Shift+C

# ABBREVIATIONS: 
# FC = faecal calprotectin; YFCCP = York FC Care Pathway; IBD = Inflammatory Bowel Disease; 
# noIBD = Irritable Bowel Syndrome (IBS) population (used instead of IBS to avoid coding errors between IBD and IBS)


#----------------------------- 0.1 Initial set up: global parameters ---------------------------------

#---- Set working directory
#setwd("Input\\Directory\\Pathway\\name")     #Set working directory as required, to where the accompanying data files have been stored (data_YFCCP_mock.csv, data_params.csv and data_outer_loop) 

#---- Install and load packages 
#This code only installs packages if not already installed, then uploads packages
if (!require("fitdistrplus"))   install.packages("fitdistrplus");    library("fitdistrplus")
if (!require("zoo"))            install.packages("zoo");             library("zoo")

#---- Remove all objects currently stored in R 
rm(list=ls())  

#---- Set test cut-off threshold (100 ug/g)
FC_cutoff     <- 100        

#---- Set number of simulations 
Nsim          <- 10000        #Base case analyses = 10000 (recommend using at least 5,000 to minimise noise)
#Nsim         <- 100000       #Sensitivity analyses [1.7 and 2.8]: Nsim = 100000

#---- Set random number seed (ensures future runs of the code produce the same output)
seed <<- 10
set.seed(seed)

#---- Specify sampling method (change as required i.e. put hash in front of analysis NOT being run)
method         <- "Bootstrap" 
#method        <- "Parametric"

#---- Specify range of imprecision and bias values to explore
Imp_values   <- seq(0, 1, by=0.005)
Bias_values  <- seq(-100, 100, length.out=length(Imp_values))   

#---- Smoothing algorithm
analysis_smooth <- "Yes"   #Base case: smoothing algorithm applied
#analysis_smooth <- "No"   #sensitivity analyses [1.6] and [2.7]: no smoothing algorithm applied

#----------------------------- 0.2 Initial set up: bootstrap method parameters ---------------------------------

if (method == "Bootstrap"){

#---- Specify method for dealing with censored data (unhash required row)
analysis_censored_values <- "replace with limits"        #Base case analysis   [1.0]: replace left-censored data ("<10") with 10, right-censored data (">600") with 600
# analysis_censored_values <- "replace with 750"         #Sensitivity analysis [1.1]: replace left-censored data with 5, right-censored data with 750
# analysis_censored_values <- "replace with 900"         #Sensitivity analysis [1.2]: replace left-censored data with 5, right-censored data with 900
# analysis_censored_values <- "replace with 1200"        #Sensitivity analysis [1.3]: replace left-censored data with 5, right-censored data with 1200
# analysis_censored_values <- "replace with 1800"        #Sensitivity analysis [1.4]: replace left-censored data with 5, right-censored data with 1800

#--- Specify whether or not to run bootstrap sampling
run_bootstrap <-  "Yes"     #Base case (runs Nsim bootstrap samples)
#run_bootstrap <-  "No"     #Sensitivity analysis [1.5]: no bootstrapping conducted (uses YFCCP n=951 data alone)
}

#----------------------------- 0.3 Initial set up: parametric method parameters ---------------------------------

if (method == "Parametric"){
  
#---- Specify upper bound for right-censored data region in fitdistrplus code
FC_maximum      <- 1000          #Base case 
#FC_maximum     <- 2000          #Sensivity analysis  [2.1]: Upper bound set to 2000 (keep base case parameterisation)
#FC_maximum     <- 3000          #Sensitvity analysis [2.2]: Upper bound set to 3000 (keep base case parameterisation)
  
#---- Specify Paramaterisation
param          <- "base case"    #Base case [2.0]: lognormal distributions used for noIBD subgroups (FC1 and FC2), Weibull used for IBD subgroups
#param         <- "lnorm"        #Sensitivity analysis [2.3]: lognormal distributions used for all subgroups
#param         <- "weib"         #Sensitivity analysis [2.4]: Weibull distributions used for all subgroups
#param         <- "gamma"        #Sensitivity analysis [2.5]: gamma distributions used for all subgroups
#param         <- "norm"         #Sensitivity analysis [2.6]: normal distributions used for all subgroups
  
#---- Specify whether to run an outer loop to account for uncertainty in the parametric parameter estimates
outer_loop <- "NO"               #Base case: no outer loop run
#outer_loop <- "YES"             #Sensitivity analysis [2.9]: Outer loop n=1000 
  
if (outer_loop=="YES"){
    Nsim_outer_loop <- 1000       #For [2.9] specify number of outer loop simulations to run (warning - this increases the time taken to run the simulation significantly - may take several days to run outer loop of 1,000 depending on computer specification)
}
}


#----------------------------- 0.4 Load YFCCP dataset -----------------------------------------------

#---- Note on Bootstrap method
#Bootstrap sampling can be conducted when using an empirical dataset, such as the YFCCP dataset in our case study
#However, to avoid issues around the sharing of individual-level data, we do not provide the YFCCP dataset here
#Instead, a mock YFCCP dataset has been generated, by drawing from the parametric distributions used in the parametric method
#Values below 10 have been replaced with "<10" and values above 600 have been replaced with ">600" to mimick the censored data seen in the original dataset
#All results of the Bootstrap example in this code will therefore be different than those reported in the manuscript due to using the mock dataset, and should be used for illustrative purposes only

#---- Note on Parametric method
#For the parametric sampling method, we provide the parameterisations derived from the YFCCP dataset using the fitdiscens code, directly, and run the error model using these specifications
#i.e. since we are not providing the YFCCP dataset, we do not provide the derivation of the parametric distributions for this data, as this requires inputting the data
#However, for illustrative purposes we provide an example of how the fitdistrplus code is run

#---- Load the mock YFCCP dataset
data <- read.csv("data_YFCCP_mock.csv")


#----------------------------- 0.5 Data set up: censored values ------------------------------------

#Temporarily set left-censored data to 5 and right-censored data to 650 to keep these values distinguishable 
#These values will be overwritten in subsequent code depending on analysis being run

temp  <- data
temp$FC1 <- as.character(temp$FC1); temp$FC2 <- as.character(temp$FC2)
temp$FC1[temp$FC1=="<10"]    <- 5   
temp$FC2[temp$FC2=="<10"]    <- 5   
temp$FC1[temp$FC1==">600"]   <- 650   
temp$FC2[temp$FC2==">600"]   <- 650   

#Set FC1 and FC2 as numeric variables
temp$FC1 <- as.numeric(temp$FC1)  
temp$FC2 <- as.numeric(temp$FC2)  

#Overwright data with new version
data <- temp      
rm(temp)  


#----------------------------- 1.0 ANALYSIS: determine the clinical accuracy of the testing pathway -------

#YFCCP pathway:
#If first test (FC1) < 100 -> treat as IBS (remain in primary care); If FC1 >= 100 -> retest.
#If FC2 <100 -> treat as IBS; If FC2 >250 -> treat as IBD (urgent); If 100<= FC2 <= 250 -> treat as IBD (routine referral)

#---- Determine "diagnosis" according to FC test result(s) 

data$FC_diagnosis <-  ifelse((data$FC1>=FC_cutoff & data$FC2>=FC_cutoff), "IBD",                
                      ifelse((data$FC1<FC_cutoff | (data$FC1>=FC_cutoff & data$FC2<FC_cutoff)),"noIBD", NA))  
#table(data$FC_diagnosis) 

#---- Determine accuracy of YFCCP vs. final clinical diagosis
data$FC_accuracy  <-  ifelse(data$FC_diagnosis=="noIBD" & data$Diagnosis=="noIBD", "TN", 
                      ifelse(data$FC_diagnosis=="noIBD" & data$Diagnosis=="IBD",   "FN",
                      ifelse(data$FC_diagnosis=="IBD"   & data$Diagnosis=="IBD",   "TP",
                      ifelse(data$FC_diagnosis=="IBD"   & data$Diagnosis=="noIBD", "FP", NA))))
#table(data$FC_accuracy)   
sens  <- sum(data$FC_accuracy == "TP", na.rm=TRUE)/(sum(data$FC_accuracy == "TP", na.rm=TRUE)+ sum(data$FC_accuracy == "FN", na.rm=TRUE))
spec  <- sum(data$FC_accuracy == "TN", na.rm=TRUE)/(sum(data$FC_accuracy == "TN", na.rm=TRUE)+ sum(data$FC_accuracy == "FP", na.rm=TRUE))
print(sens); print(spec)   
#The YFCCP mock data produces a baseline sensitivity of 0.9487179 and specificity of 0.9106529 (compared to 0.94 & 0.92 in the original dataset)


#----------------------------- 2.0 ANALYSIS: error model simulation using bootstrap sampling method ----------------    
if(method=="Bootstrap"){
  
#---- Time the simulation
time_start <- proc.time()
  
#---- create arrays to store simulation loop results
Loop_sens  <- array(NA, dim=c(length(Imp_values),length(Bias_values))) 
Loop_spec  <- array(NA, dim=c(length(Imp_values),length(Bias_values))) 

#---- Copy data for simulation
data_sim <- data

#---- Replace censored values with selected option
  
if(analysis_censored_values=="replace with limits"){
    data_sim$FC1[data_sim$FC1==650]   <- 600 
    data_sim$FC1[data_sim$FC1==5]     <- 10
    data_sim$FC2[!(is.na(data_sim$FC2)) & data_sim$FC2==650]   <- 600 
    data_sim$FC2[!(is.na(data_sim$FC2)) & data_sim$FC2==5]     <- 10
}

#Nb: left-censored data already set to 5 for all thee sensitivity analyses:
if(analysis_censored_values=="replace with 750"){  
    data_sim$FC1[data_sim$FC1==650]   <- 750 
    data_sim$FC2[!(is.na(data_sim$FC2)) & data_sim$FC2==650]   <- 750 
}
  
if(analysis_censored_values=="replace with 900"){
    data_sim$FC1[data_sim$FC1==650]   <- 900 
    data_sim$FC2[!(is.na(data_sim$FC2)) & data_sim$FC2==650]   <- 900 
}
  
if(analysis_censored_values=="replace with 1200"){
    data_sim$FC1[data_sim$FC1==650]   <- 1200 
    data_sim$FC2[!(is.na(data_sim$FC2)) & data_sim$FC2==650]   <- 1200 
}
  
if(analysis_censored_values=="replace with 1800"){
    data_sim$FC1[data_sim$FC1==650]   <- 1800 
    data_sim$FC2[!(is.na(data_sim$FC2)) & data_sim$FC2==650]   <- 1800 
}
  
#----- Bootstrap sampling
#Use 'sample' function to sample dataset rows with replacement
if(run_bootstrap=="Yes"){
data_sim <- data_sim[sample(nrow(data_sim),size=Nsim, replace=TRUE),] 
}
  
#------ Error model simulation 
for (i in 1:length(Imp_values)){
for (b in 1:length(Bias_values)) {
      
#-- First error model application: FC1 values                
#Generate FC1 values including additional imprecision and bias (g_FC1) from the baseline FC1 values
data_sim$g_FC1 <- data_sim$FC1 + data_sim$FC1*rnorm(length(data_sim$FC1),0,1)*Imp_values[i] + Bias_values[b]
      
#-- FC2 data: get complete FC2 data by supplementing with random sampling from YFCCP FC2 data where necessary  
data_sim$FC2_full <- NA
#Populate with available data:
data_sim$FC2_full[data_sim$g_FC1>=FC_cutoff & !is.na(data_sim$FC2)] <- data_sim$FC2[data_sim$g_FC1>=FC_cutoff & !is.na(data_sim$FC2)]
#Resample to fill in data gaps:
data_sim$FC2_full[data_sim$g_FC1>=FC_cutoff & is.na(data_sim$FC2) & data_sim$Diagnosis=="noIBD"] <-  sample(data_sim$FC2[complete.cases(data_sim$FC2) & data_sim$Diagnosis=="noIBD"], 
                                                                                                            size= length(data_sim$FC2_full[data_sim$g_FC1>=FC_cutoff & is.na(data_sim$FC2) & data_sim$Diagnosis=="noIBD"]),
                                                                                                            replace=TRUE)
data_sim$FC2_full[data_sim$g_FC1>=FC_cutoff & is.na(data_sim$FC2) & data_sim$Diagnosis=="IBD"]   <-  sample(data_sim$FC2[complete.cases(data_sim$FC2) & data_sim$Diagnosis=="IBD"], 
                                                                                                            size= length(data_sim$FC2_full[data_sim$g_FC1>=FC_cutoff & is.na(data_sim$FC2) & data_sim$Diagnosis=="IBD"]),
                                                                                                            replace=TRUE)
#--Second error model application: FC2 values
data_sim$g_FC2 <- ifelse(data_sim$g_FC1<FC_cutoff, NA, (data_sim$FC2_full + data_sim$FC2_full*rnorm(length(data_sim$FC2_full),0,1)*Imp_values[i] + Bias_values[b]))
        
#--Calculate FC diagnosis based on simulated results, as per YFCCP protocol
data_sim$g_FC_diagnosis <-  ifelse(data_sim$g_FC1>=FC_cutoff & data_sim$g_FC2>=FC_cutoff, "IBD",           
                            ifelse(data_sim$g_FC1<FC_cutoff | (data_sim$g_FC1>=FC_cutoff & data_sim$g_FC2<FC_cutoff),"noIBD", NA))    

#--Calculate sensitivity and specificity of simulated FC values     
data_sim$g_FC_accuracy   <- ifelse(data_sim$g_FC_diagnosis=="noIBD" & data_sim$Diagnosis=="noIBD", "TN", 
                            ifelse(data_sim$g_FC_diagnosis=="noIBD" & data_sim$Diagnosis=="IBD",   "FN",
                            ifelse(data_sim$g_FC_diagnosis=="IBD"   & data_sim$Diagnosis=="IBD",   "TP",
                            ifelse(data_sim$g_FC_diagnosis=="IBD"   & data_sim$Diagnosis=="noIBD", "FP", NA))))
#table(data_sim$g_FC_accuracy)
sens  <- sum(data_sim$g_FC_accuracy == "TP", na.rm=TRUE)/(sum(data_sim$g_FC_accuracy == "TP", na.rm=TRUE)+ sum(data_sim$g_FC_accuracy == "FN", na.rm=TRUE))
spec  <- sum(data_sim$g_FC_accuracy == "TN", na.rm=TRUE)/(sum(data_sim$g_FC_accuracy == "TN", na.rm=TRUE)+ sum(data_sim$g_FC_accuracy == "FP", na.rm=TRUE))

#--Populate results array
Loop_sens[i,b]  <- sens     
Loop_spec[i,b]  <- spec

}}

#--Assign row and column names
rownames(Loop_sens) <- Imp_values; colnames(Loop_sens) <- Bias_values
rownames(Loop_spec) <- Imp_values; colnames(Loop_spec) <- Bias_values
  
#--Transpose matrices so bias is on the x axis for plots
Loop_sens <- t(Loop_sens)  
Loop_spec <- t(Loop_spec)
  
#--Print baseline results (at bias and CV% = 0) [YFCCP baseline sensitivity =94%, specificity = 92%]
print(Loop_sens[match(0,Bias_values),match(0,Imp_values)])    
print(Loop_spec[match(0,Bias_values),match(0,Imp_values)])      
  
time_elapsed <- proc.time() - time_start
} 
#end of code for Bootstrap/raw data model

#print(time_elapsed)  #Run to see how long the simulation took



#----------------------------- 3.0 ANALYSIS: error model simulation using parametric sampling method ----------------    

#For guidance on fitdistrplus package used in this code see file:///C:/Users/hssasm/Downloads/v64i04%20(2).pdf

##---------- Example of distribution derivation (for reference only, not needed to run the simulation) ----
# #Example provided for FC1 noIBD population
# 
# #-- Define "left" and "right" vectors specifying the region for left- and right-censored data
# #For complete data, data$left and data$right are equal to the data point:
# data$left   <- data$FC1
# data$right  <- data$FC1
# 
# #For left-censored data, data$left=0 and data$right=10 (i.e. the region spans from 0 to 10)
# data$left[data$left==5]   <- 0               
# data$right[data$left==0]  <- 10
# 
# #For right-censored data, data$left=600 and data$right=FC_maximum (1000 in the base case)
# data$right[data$left==650] <- FC_maximum   
# data$left[data$left==650]  <- 600
# 
# #-- Subset data to no_IBD population
# data_noIBD       <- data[data$Diagnosis=="noIBD",]     #overwright previous dataframe to include accuracy and left/right data
# temp             <- data.frame(data_noIBD$left, data_noIBD$right)
# names (temp)[c(1,2)] <- c("left", "right")
# 
# #Plot CDF of censored data (turnbull plot)
# #plotdistcens(temp)
# 
# #-- Explore different parameterizations (use summary to obtain the mean distribution parameters and AIC and BIC criteria)
# fit_noIBD_FC1_gamma   <- fitdistcens(temp, "gamma")    #;summary(fit_noIBD_FC1_gamma)
# fit_noIBD_FC1_weib    <- fitdistcens(temp, "weibull")  #;summary(fit_noIBD_FC1_weib)
# fit_noIBD_FC1_ln      <- fitdistcens(temp, "lnorm")    #;summary(fit_noIBD_FC1_ln)
# fit_noIBD_FC1_norm    <- fitdistcens(temp, "norm")     #;summary(fit_noIBD_FC1_norm)
# #Lognormal performs best (lowest AIC & BIC)
# #Note these values will not exactly match those in the data_params file due to using the mock database
# 
# #PLot parametric fits onto censored data CDF plot:
# #cdfcompcens(list(fit_noIBD_FC1_gamma, fit_noIBD_FC1_weib, fit_noIBD_FC1_ln, fit_noIBD_FC1_norm), legendtext = c("gamma", "weibull", "lognormal", "normal"), xlab="FC1 result")
# 
# #-- Simulate frequency distributions 
# sim_noIBD_FC1_gamma   <- rgamma  (Nsim,  shape=   fit_noIBD_FC1_gamma$estimate[1], rate=  fit_noIBD_FC1_gamma$estimate[2])
# sim_noIBD_FC1_weib    <- rweibull(Nsim,  shape=   fit_noIBD_FC1_weib$estimate[1],  scale= fit_noIBD_FC1_weib$estimate[2])
# sim_noIBD_FC1_lnorm   <- rlnorm  (Nsim,  meanlog= fit_noIBD_FC1_ln$estimate[1],    sdlog= fit_noIBD_FC1_ln$estimate[2])
# sim_noIBD_FC1_norm    <- rnorm   (Nsim,  mean=    fit_noIBD_FC1_norm$estimate[1],  sd=    fit_noIBD_FC1_norm$estimate[2])
# 
# #Plot frequency density distributions
# # plot (density(data_noIBD$FC1), lwd=2, xlim=c(-200,1000), main="", xlab="Initial FC Result (FC1)", ylab="Density")
# # lines(density(sim_noIBD_FC1_gamma), col="seagreen2",   lty="longdash", lwd=3)
# # lines(density(sim_noIBD_FC1_weib),  col="purple",      lty="twodash",  lwd=3)
# # lines(density(sim_noIBD_FC1_lnorm), col="deepskyblue", lty="dotdash",  lwd=3)
# # lines(density(sim_noIBD_FC1_norm),  col="gold2"      , lty="dotted",   lwd=3)
# # legend(600,0.012, legend=c("Raw data", "Gamma fit", "Weibull fit", "Lognormal fit", "Normal fit"), col=c("black", "seagreen2", "purple", "deepskyblue", "gold2"), lty=c("solid", "longdash", "twodash","dotdash","dotted"), lwd=c(2,2,2,2,2), cex=1.2)
# 
# #Remove simulation data
# rm(temp, sim_noIBD_FC1_gamma, sim_noIBD_FC1_weib, sim_noIBD_FC1_lnorm, sim_noIBD_FC1_norm)


#----------- Simulation code ----

if(method=="Parametric"){

#---- Subset data into IBD and noIBD groups
data_IBD   <- data[data$Diagnosis=="IBD",]
data_noIBD <- data[data$Diagnosis=="noIBD",]

#---- Time the simulation
time_start <- proc.time()

#----- Specify number of samples for each population
prev        <- nrow(data_IBD)/nrow(data)   #prevalence of IBD in dataset 
Nsim_IBD    <- round(Nsim*prev)
Nsim_noIBD  <- round(Nsim*(1-prev))
#Nsim_IBD ;Nsim_noIBD

#----- create arrays to store simulation loop results
Loop_sens  <- array(NA, dim=c(length(Imp_values),length(Bias_values)))
Loop_spec  <- array(NA, dim=c(length(Imp_values),length(Bias_values)))

#----- Error model simulation (base case - no outer loop)
if(outer_loop=="NO"){
  
#---- Read in parametric specifications data
data_params <- read.csv("data_params.csv")

#------ Define FC1 parametric distribution (FC2 specified within simulation loop)
#Using the parameter specifications provided in the data_params file, simulate distributions for IBD and noIBD subgroups

if (param=="base case") {
#Base case = lnorm parameterization for noIBD populations, and Weibull for IBD
temp_noIBD  <- rlnorm(Nsim_noIBD, data_params$noIBD_FC1_param1[data_params$Parameterization=="Lnorm" & data_params$FC_upper==FC_maximum], 
                                  data_params$noIBD_FC1_param2[data_params$Parameterization=="Lnorm" & data_params$FC_upper==FC_maximum])
temp_IBD    <- rweibull(Nsim_IBD, data_params$IBD_FC1_param1[data_params$Parameterization=="Weib"    & data_params$FC_upper==FC_maximum],   
                                  data_params$IBD_FC1_param2[data_params$Parameterization=="Weib"    & data_params$FC_upper==FC_maximum])
}
if (param=="lnorm") {
temp_noIBD  <- rlnorm(Nsim_noIBD, data_params$noIBD_FC1_param1[data_params$Parameterization=="Lnorm" & data_params$FC_upper==FC_maximum], 
                                  data_params$noIBD_FC1_param2[data_params$Parameterization=="Lnorm" & data_params$FC_upper==FC_maximum])
temp_IBD    <- rlnorm(Nsim_IBD,   data_params$IBD_FC1_param1[data_params$Parameterization=="Lnorm"   & data_params$FC_upper==FC_maximum],   
                                  data_params$IBD_FC1_param2[data_params$Parameterization=="Lnorm"   & data_params$FC_upper==FC_maximum])
}
if (param=="weib") {
temp_noIBD  <- rlnorm(Nsim_noIBD, data_params$noIBD_FC1_param1[data_params$Parameterization=="Weib" & data_params$FC_upper==FC_maximum], 
                                  data_params$noIBD_FC1_param2[data_params$Parameterization=="Weib" & data_params$FC_upper==FC_maximum])
temp_IBD    <- rweibull(Nsim_IBD, data_params$IBD_FC1_param1[data_params$Parameterization=="Weib"   & data_params$FC_upper==FC_maximum],   
                                  data_params$IBD_FC1_param2[data_params$Parameterization=="Weib"   & data_params$FC_upper==FC_maximum])
}
if (param=="gamma") {
temp_noIBD  <- rgamma(Nsim_noIBD, data_params$noIBD_FC1_param1[data_params$Parameterization=="Gamma" & data_params$FC_upper==FC_maximum], 
                                  data_params$noIBD_FC1_param2[data_params$Parameterization=="Gamma" & data_params$FC_upper==FC_maximum])
temp_IBD    <- rgamma(Nsim_IBD,   data_params$IBD_FC1_param1[data_params$Parameterization=="Gamma"   & data_params$FC_upper==FC_maximum],   
                                  data_params$IBD_FC1_param2[data_params$Parameterization=="Gamma"   & data_params$FC_upper==FC_maximum])
}
if (param=="norm") {
temp_noIBD  <- rnorm(Nsim_noIBD,  data_params$noIBD_FC1_param1[data_params$Parameterization=="Norm" & data_params$FC_upper==FC_maximum], 
                                  data_params$noIBD_FC1_param2[data_params$Parameterization=="Norm" & data_params$FC_upper==FC_maximum])
temp_IBD    <- rnorm(Nsim_IBD,    data_params$IBD_FC1_param1[data_params$Parameterization=="Norm"   & data_params$FC_upper==FC_maximum],   
                                  data_params$IBD_FC1_param2[data_params$Parameterization=="Norm"   & data_params$FC_upper==FC_maximum])
}
    
#---- Create dataframe of FC1 and Diagnosis (combine IBD and noIBD data)
temp1  <- data.frame(temp_noIBD, rep("noIBD", times=Nsim_noIBD))   
temp2  <- data.frame(temp_IBD,   rep("IBD",   times=Nsim_IBD))
names(temp1)[c(1,2)]  <- c("FC1","Diagnosis")
names(temp2)[c(1,2)]  <- c("FC1","Diagnosis")
data_sim <- rbind(temp1, temp2)

rm(temp_noIBD, temp_IBD, temp1, temp2)

#------ Error model simulation 
for (i in 1:length(Imp_values)){
for (b in 1:length(Bias_values)) {
        
#-- First error model application: FC1
data_sim$g_FC1 <- data_sim$FC1 + data_sim$FC1*rnorm(length(data_sim$FC1),0,1)*Imp_values[i] + Bias_values[b]
        
#-- Generate FC2 values for all g_FC1 values>100 by population
if(param=="base case"){
temp_noIBD    <- rlnorm  (length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="noIBD"))), 
                          data_params$noIBD_FC2_param1[data_params$Parameterization=="Lnorm" & data_params$FC_upper==FC_maximum], 
                          data_params$noIBD_FC2_param2[data_params$Parameterization=="Lnorm" & data_params$FC_upper==FC_maximum])
temp_IBD      <- rweibull(length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="IBD"))),   
                          data_params$IBD_FC2_param1[data_params$Parameterization=="Weib" & data_params$FC_upper==FC_maximum], 
                          data_params$IBD_FC2_param2[data_params$Parameterization=="Weib" & data_params$FC_upper==FC_maximum])
}

if(param=="lnorm"){
temp_noIBD    <- rlnorm   (length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="noIBD"))), 
                          data_params$noIBD_FC2_param1[data_params$Parameterization=="Lnorm" & data_params$FC_upper==FC_maximum], 
                          data_params$noIBD_FC2_param2[data_params$Parameterization=="Lnorm" & data_params$FC_upper==FC_maximum])
temp_IBD      <- rlnorm   (length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="IBD"))),   
                          data_params$IBD_FC2_param1[data_params$Parameterization=="Lnorm" & data_params$FC_upper==FC_maximum], 
                          data_params$IBD_FC2_param2[data_params$Parameterization=="Lnorm" & data_params$FC_upper==FC_maximum])
}           

if(param=="weib"){
temp_noIBD    <- rweibull(length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="noIBD"))), 
                          data_params$noIBD_FC2_param1[data_params$Parameterization=="Weib" & data_params$FC_upper==FC_maximum], 
                          data_params$noIBD_FC2_param2[data_params$Parameterization=="Weib" & data_params$FC_upper==FC_maximum])
temp_IBD      <- rweibull(length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="IBD"))),   
                          data_params$IBD_FC2_param1[data_params$Parameterization=="Weib" & data_params$FC_upper==FC_maximum], 
                          data_params$IBD_FC2_param2[data_params$Parameterization=="Weib" & data_params$FC_upper==FC_maximum])
}

if(param=="gamma"){
temp_noIBD    <- rgamma  (length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="noIBD"))), 
                          data_params$noIBD_FC2_param1[data_params$Parameterization=="Gamma" & data_params$FC_upper==FC_maximum], 
                          data_params$noIBD_FC2_param2[data_params$Parameterization=="Gamma" & data_params$FC_upper==FC_maximum])
temp_IBD      <- rgamma  (length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="IBD"))),   
                          data_params$IBD_FC2_param1[data_params$Parameterization=="Gamma" & data_params$FC_upper==FC_maximum], 
                          data_params$IBD_FC2_param2[data_params$Parameterization=="Gamma" & data_params$FC_upper==FC_maximum])
}

if(param=="norm"){
temp_noIBD    <- rnorm   (length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="noIBD"))), 
                          data_params$noIBD_FC2_param1[data_params$Parameterization=="Norm" & data_params$FC_upper==FC_maximum], 
                          data_params$noIBD_FC2_param2[data_params$Parameterization=="Norm" & data_params$FC_upper==FC_maximum])
temp_IBD      <- rnorm   (length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="IBD"))),   
                          data_params$IBD_FC2_param1[data_params$Parameterization=="Norm" & data_params$FC_upper==FC_maximum], 
                          data_params$IBD_FC2_param2[data_params$Parameterization=="Norm" & data_params$FC_upper==FC_maximum])
}
            
#-- Fill data_sim FC2 values with fitted data
data_sim$FC2  <- ifelse(data_sim$g_FC1<FC_cutoff, NA, ifelse(data_sim$Diagnosis=="noIBD", temp_noIBD, temp_IBD))
            
#-- Second error model application: FC2 
data_sim$g_FC2 <- NA
data_sim$g_FC2 <- as.numeric(data_sim$g_FC2)
#We want to ignore the empty FC2 cells, so conduct error model only on the complete cases for FC2 using the 'complete.cases()' function:
data_sim$g_FC2[complete.cases(data_sim$FC2)] <- data_sim$FC2[complete.cases(data_sim$FC2)] + data_sim$FC2[complete.cases(data_sim$FC2)]*rnorm(length(data_sim$FC2[complete.cases(data_sim$FC2)]),0,1)*Imp_values[i] + Bias_values[b]
            
#-- Calculate FC diagnoses as per YFCCP protocol
data_sim$g_FC_diagnosis <-  ifelse(data_sim$g_FC1>=FC_cutoff & data_sim$g_FC2>=FC_cutoff, "IBD",           
                            ifelse((data_sim$g_FC1<FC_cutoff |(data_sim$g_FC1>=FC_cutoff & data_sim$g_FC2<FC_cutoff)),"noIBD", NA))     
        
#-- Calculate clinical accuracy based on FC simulated values
data_sim$g_FC_accuracy   <- ifelse(data_sim$g_FC_diagnosis=="noIBD" & data_sim$Diagnosis=="noIBD", "TN",
                            ifelse(data_sim$g_FC_diagnosis=="noIBD" & data_sim$Diagnosis=="IBD",   "FN",
                            ifelse(data_sim$g_FC_diagnosis=="IBD"   & data_sim$Diagnosis=="IBD",   "TP",
                            ifelse(data_sim$g_FC_diagnosis=="IBD"   & data_sim$Diagnosis=="noIBD", "FP", NA))))
#table(data_sim$g_FC_accuracy)
sens  <- sum(data_sim$g_FC_accuracy == "TP", na.rm=TRUE)/(sum(data_sim$g_FC_accuracy == "TP", na.rm=TRUE)+ sum(data_sim$g_FC_accuracy == "FN", na.rm=TRUE))
spec  <- sum(data_sim$g_FC_accuracy == "TN", na.rm=TRUE)/(sum(data_sim$g_FC_accuracy == "TN", na.rm=TRUE)+ sum(data_sim$g_FC_accuracy == "FP", na.rm=TRUE))
Loop_sens[i,b]  <- sens     
Loop_spec[i,b]  <- spec
        
}}
    
}
#end of error model loop (no outer loop)

#----- Error model simulation (sensitivity analysis [2.9]: with outer loop simulation to capture parametric uncertainty)
  
if(outer_loop=="YES"){

#----- create arrays to store loop results with outer loop
Loop_sens_outer  <- array(NA, dim=c(length(Imp_values),length(Bias_values),Nsim_outer_loop))
Loop_spec_outer  <- array(NA, dim=c(length(Imp_values),length(Bias_values),Nsim_outer_loop))

#---- Read in outer loop data
data_outer_loop <- read.csv("data_outer_loop.csv")
    
#---- Define FC1 parametric distribution 
for(z in 1:Nsim_outer_loop){
temp_noIBD  <- rlnorm( Nsim_noIBD, data_outer_loop$noIBD_FC1_est1[z], data_outer_loop$noIBD_FC1_est2[z])
temp_IBD    <- rweibull(Nsim_IBD,  data_outer_loop$IBD_FC1_est1[z],   data_outer_loop$IBD_FC1_est2[z])
      
#Create dataframe of FC1 and diagnosis
temp1  <- data.frame(temp_noIBD, rep("noIBD", times=Nsim_noIBD))   
temp2  <- data.frame(temp_IBD,   rep("IBD",   times=Nsim_IBD))
names(temp1)[c(1,2)]  <- c("FC1","Diagnosis")
names(temp2)[c(1,2)]  <- c("FC1","Diagnosis")
data_sim <- rbind(temp1, temp2)
      
rm(temp_noIBD, temp_IBD, temp1, temp2)
      
#---- Inner loop to apply error model 
for (i in 1:length(Imp_values)){
for (b in 1:length(Bias_values)) {
          
#---- First error model application: FC1
data_sim$g_FC1 <- data_sim$FC1 + data_sim$FC1*rnorm(length(data_sim$FC1),0,1)*Imp_values[i] + Bias_values[b]
          
#---- Generate FC2 values for all g_FC1 values>100 for noIBD and IBD populations:
temp_noIBD    <- rlnorm(length(which(  (data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="noIBD"))), data_outer_loop$noIBD_FC2_est1[z], data_outer_loop$noIBD_FC2_est2[z]) 
temp_IBD      <- rweibull(length(which((data_sim$g_FC1>=FC_cutoff) & (data_sim$Diagnosis=="IBD"))),   data_outer_loop$IBD_FC2_est1[z],   data_outer_loop$IBD_FC2_est2[z])
            
#---- Fill data_sim FC2 values with fitted data
data_sim$FC2  <- ifelse(data_sim$g_FC1<FC_cutoff, NA, ifelse(data_sim$Diagnosis=="noIBD", temp_noIBD, temp_IBD))
            
#---- Second error model appliction: FC2
data_sim$g_FC2 <- NA
data_sim$g_FC2 <- as.numeric(data_sim$g_FC2)
data_sim$g_FC2[complete.cases(data_sim$FC2)] <- data_sim$FC2[complete.cases(data_sim$FC2)] + data_sim$FC2[complete.cases(data_sim$FC2)]*rnorm(length(data_sim$FC2[complete.cases(data_sim$FC2)]),0,1)*Imp_values[i] + Bias_values[b]
            
#---- Calculate FC diagnoses as per YFCCP protocol
data_sim$g_FC_diagnosis <-  ifelse(data_sim$g_FC1>=FC_cutoff & data_sim$g_FC2>=FC_cutoff, "IBD",           
                            ifelse((data_sim$g_FC1<FC_cutoff |(data_sim$g_FC1>=FC_cutoff & data_sim$g_FC2<FC_cutoff)),"noIBD", NA))     

#---- Calculate clinical accuracy based on FC simulated values
data_sim$g_FC_accuracy   <- ifelse(data_sim$g_FC_diagnosis=="noIBD" & data_sim$Diagnosis=="noIBD", "TN",
                            ifelse(data_sim$g_FC_diagnosis=="noIBD" & data_sim$Diagnosis=="IBD",   "FN",
                            ifelse(data_sim$g_FC_diagnosis=="IBD"   & data_sim$Diagnosis=="IBD",   "TP",
                            ifelse(data_sim$g_FC_diagnosis=="IBD"   & data_sim$Diagnosis=="noIBD", "FP", NA))))
#table(data_sim$g_FC_accuracy)
sens  <- sum(data_sim$g_FC_accuracy == "TP", na.rm=TRUE)/(sum(data_sim$g_FC_accuracy == "TP", na.rm=TRUE)+ sum(data_sim$g_FC_accuracy == "FN", na.rm=TRUE))
spec  <- sum(data_sim$g_FC_accuracy == "TN", na.rm=TRUE)/(sum(data_sim$g_FC_accuracy == "TN", na.rm=TRUE)+ sum(data_sim$g_FC_accuracy == "FP", na.rm=TRUE))
Loop_sens_outer[i,b,z]  <- sens     
Loop_spec_outer[i,b,z]  <- spec   
}}}

#Take mean of outer loops   
Loop_sens <- apply(Loop_sens_outer, c(1,2), mean)
Loop_spec <- apply(Loop_spec_outer, c(1,2), mean)

}#end of simulation loop (with outer loop)
   
#---- Add row and column names
rownames(Loop_sens) <- Imp_values; colnames(Loop_sens) <- Bias_values
rownames(Loop_spec) <- Imp_values; colnames(Loop_spec) <- Bias_values

#---- Transpose matrices so bias is on the x axis for plots
Loop_sens <- t(Loop_sens)  
Loop_spec <- t(Loop_spec)

print(Loop_sens[match(0,Bias_values),match(0,Imp_values)])    
print(Loop_spec[match(0,Bias_values),match(0,Imp_values)])    

time_elapsed <- proc.time() - time_start    

}#end of code for parametric method

#print(time_elapsed)



#----------------------------- 4.0 RESULTS:  smoothing algorithm --------

if(analysis_smooth == "Yes"){
  
#----- Define windows over which to calculate moving average values
window_main  <- 10
window_edges <- 5

#----- Clinical sensitivity
#Set up several temp arrays over which different moving average windows will be applied:
temp_sens   <- Loop_sens; temp_sens_1 <- Loop_sens; temp_sens_2 <- Loop_sens; temp_sens_3 <- Loop_sens; temp_sens_4 <- Loop_sens

#-- For main data apply moving average over window of 10:
for(i in 1:ncol(temp_sens_1)){ temp_sens_1[,i] <- rollmean(temp_sens_1[,i],window_main , fill=NA)}

#-- For values unaddressed by the above code, calculate moving average over of window of 5, starting with central alignment, then left, then right:
for(i in 1:ncol(temp_sens_2)){ temp_sens_2[,i] <- rollmean(temp_sens_2[,i],window_edges, fill=NA)}
for(i in 1:ncol(temp_sens_3)){ temp_sens_3[,i] <- rollmean(temp_sens_3[,i],window_edges, align="left",fill=NA)}
for(i in 1:ncol(temp_sens_4)){ temp_sens_4[,i] <- rollmean(temp_sens_4[,i],window_edges, align="right",fill=NA)}

#-- Fill in values using temp_sens_1 first, then followed by the other moving average calculations where necessary:
temp      <- ifelse(is.na(temp_sens_1),temp_sens_2,temp_sens_1)
temp_sens <- ifelse(is.na(temp),ifelse(is.na(temp_sens_3),temp_sens_4, temp_sens_3), temp)

#-- Transpose the data and re-do the moving average calculations:
temp_sens <- t(temp_sens)
temp_sens_1 <- temp_sens; temp_sens_2 <- temp_sens; temp_sens_3 <- temp_sens; temp_sens_4 <- temp_sens
for(i in 1:ncol(temp_sens_1)){ temp_sens_1[,i] <- rollmean(temp_sens_1[,i],window_main , fill=NA)}
for(i in 1:ncol(temp_sens_2)){ temp_sens_2[,i] <- rollmean(temp_sens_2[,i],window_edges, fill=NA)}
for(i in 1:ncol(temp_sens_3)){ temp_sens_3[,i] <- rollmean(temp_sens_3[,i],window_edges, align="left",fill=NA)}
for(i in 1:ncol(temp_sens_4)){ temp_sens_4[,i] <- rollmean(temp_sens_4[,i],window_edges, align="right",fill=NA)}
temp <- ifelse(is.na(temp_sens_1),temp_sens_2,temp_sens_1)
temp_sens <- ifelse(is.na(temp),ifelse(is.na(temp_sens_3),temp_sens_4, temp_sens_3), temp)
temp_sens <- t(temp_sens)

#----- Clinical specificity (same process as for sensitivity)
temp_spec   <- Loop_spec; temp_spec_1 <- Loop_spec; temp_spec_2 <- Loop_spec; temp_spec_3 <- Loop_spec; temp_spec_4 <- Loop_spec
for(i in 1:ncol(temp_spec_1)){ temp_spec_1[,i] <- rollmean(temp_spec_1[,i],window_main , fill=NA)}
for(i in 1:ncol(temp_spec_2)){ temp_spec_2[,i] <- rollmean(temp_spec_2[,i],window_edges, fill=NA)}
for(i in 1:ncol(temp_spec_3)){ temp_spec_3[,i] <- rollmean(temp_spec_3[,i],window_edges, align="left",fill=NA)}
for(i in 1:ncol(temp_spec_4)){ temp_spec_4[,i] <- rollmean(temp_spec_4[,i],window_edges, align="right",fill=NA)}
temp <- ifelse(is.na(temp_spec_1),temp_spec_2,temp_spec_1)
temp_spec <- ifelse(is.na(temp),ifelse(is.na(temp_spec_3),temp_spec_4, temp_spec_3), temp)

temp_spec <- t(temp_spec)
temp_spec_1 <- temp_spec; temp_spec_2 <- temp_spec; temp_spec_3 <- temp_spec; temp_spec_4 <- temp_spec; 
for(i in 1:ncol(temp_spec_1)){ temp_spec_1[,i] <- rollmean(temp_spec_1[,i],window_main , fill=NA)}
for(i in 1:ncol(temp_spec_2)){ temp_spec_2[,i] <- rollmean(temp_spec_2[,i],window_edges, fill=NA)}
for(i in 1:ncol(temp_spec_3)){ temp_spec_3[,i] <- rollmean(temp_spec_4[,i],window_edges, align="left",fill=NA)}
for(i in 1:ncol(temp_spec_4)){ temp_spec_4[,i] <- rollmean(temp_spec_4[,i],window_edges, align="right",fill=NA)}
temp <- ifelse(is.na(temp_spec_1),temp_spec_2,temp_spec_1)
temp_spec <- ifelse(is.na(temp),ifelse(is.na(temp_spec_3),temp_spec_4, temp_spec_3), temp)
temp_spec <- t(temp_spec)

#----- Define smooth arrays for loop results
Loop_sens_smooth <- temp_sens
Loop_spec_smooth <- temp_spec

rm(temp, temp_sens, temp_spec, temp_sens_1, temp_sens_2, temp_sens_3, temp_sens_4, temp_spec_1, temp_spec_2, temp_spec_3, temp_spec_4)
#End of smoothing code
}

#----------------------------- 5.0 RESULTS:  save files ---------------------------------------

#Set working directory (if want to save files to a different location as current directory)
#setwd("Select\\Directory\\Pathway")

#----- Save results 
#write.csv(Loop_sens, "Loop_sens.csv")                #Nb: these are the raw "noisy" versions of the results [Sensitivity analyses 1.6 and 2.7]
#write.csv(Loop_spec, "Loop_spec.csv")
#write.csv(Loop_sens_smooth, "Loop_sens_smooth.csv")  #Nb: these are the "smooth" versions of the results (base case)
#write.csv(Loop_spec_smooth, "Loop_spec_smooth.csv")
##Note: remember to rename/move these files as necessary when running sensitivity analyses so existing files do not get overwritten


#----------------------------- 6.0 RESULTS:  contour plots --------------------------------

#For notes on contour function used for these plots see: 
#https://www.math.ucla.edu/~anderson/rw1001/library/base/html/contour.html

#---- Import files (if composing plots at a later date to analyses) 
#Loop_sens_smooth <- read.csv("File\\Pathway", header=TRUE, row.names=1)  #Change file pathway(s) as required
#Loop_spec_smooth <- read.csv("File\\Pathway", header=TRUE, row.names=1)
#Loop_sens        <- read.csv("File\\Pathway", header=TRUE, row.names=1)
#Loop_spec        <- read.csv("File\\Pathway", header=TRUE, row.names=1)


#----- Get data into appropriate format for plotting  

#-- Convert required results to matrix format
if(analysis_smooth == "Yes"){
Loop_sens <- as.matrix(Loop_sens_smooth)
Loop_spec <- as.matrix(Loop_spec_smooth) }

if(analysis_smooth =="No") {
Loop_sens <- as.matrix(Loop_sens)
Loop_spec <- as.matrix(Loop_spec) }

#-- Add column names (required for contour plots) & 
#-- Convert imprecision to % (0-100) scale
Imp_values  <- seq(0, 100, by=0.5)
colnames(Loop_sens) <- Imp_values; colnames(Loop_spec) <- Imp_values


#----- Contour figure ----

#-- Run 'tiff(...)' and 'dev.off()' lines of code to save plot as tiff file in current working directory
tiff(filename="Contour_plot_temp.tiff", units="in", width=10, height=7, res=300) #width and height defined in inches, resolution in dots per inch

par(mfrow = c(1, 1))  #Specify one plot in the figure (can change to have mutliple plots in one figure)
par(mar=c(4.5,5,2,2)) #Specify space around the plot (bottom, left, top, right) (extra space may be needed for axis labels)

#-- Plot initial contour of sensitivity results
contour(x=as.numeric(rownames(Loop_sens)), y=as.numeric(colnames(Loop_sens)), Loop_sens, levels=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95) , xlab=expression(paste(Delta," Bias")), ylab=expression(paste(Delta," Imprecision (CV%)")), col="deepskyblue", labcex=1.2, cex.lab=1.2, cex.axis=1.2, lwd=2 ,drawlabels=FALSE)  

#-- Add grid lines to plot
abline(h=seq(from=0,to=100,by=20), v=seq(from=-100,to=100,by=25),col="grey91",lty="solid")  #Use abline to specify exact placement of grid lines

#-- Reinstate contour lines over grid lines
contour(x=as.numeric(rownames(Loop_sens)), y=as.numeric(colnames(Loop_sens)), Loop_sens, levels=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95) , add=TRUE, col="deepskyblue",  labcex=1.2, cex.lab=1.2, cex.axis=1.2, lwd=2, drawlabels=TRUE, method="flattest")               # change to drawlabels=FALSE if inserting labels via manual text additions
contour(x=as.numeric(rownames(Loop_spec)), y=as.numeric(colnames(Loop_spec)), Loop_spec, levels=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),  add=TRUE, col="midnightblue", labcex=1.2, cex.lab=1.2, cex.axis=1.2, lwd=2, lty="dashed", drawlabels=TRUE, method="flattest") # change to drawlabels=FALSE if inserting labels via manual text additions   

#-- Legend 
legend("topleft", inset=c(0.04,0.1), horiz=FALSE, c("Sensitivity", "Specificity"), lty=c("solid","dashed"), col=c("deepskyblue","midnightblue"), lwd=c(3,3), cex=1.2)

#-- Add manual contour line labels as required (the user will need to update these values depending on the specific plot as the line positions change):
#text(-42,99,  "0.95",cex=1, pos=3, col="midnightblue") #Example
#text(12,99,   "0.9",cex=1, pos=3, col="midnightblue")  #Etc. 
#Alternatively within the contour function, add "drawlabels = TRUE" (but adding labels as text allows greater control over their placement & improves readability)

#-- Save file
dev.off()


#----- Contour figure + acceptable region + TE bands ----  

#-- Specify minimum outcome requirements for acceptable region (AR)
#If these values are changed, the plot legend will also require updating 
AR_min_sens <- 0.85    #BC = 0.85 (lower CI from data). 
AR_min_spec <- 0.90    #BC = 0.90 (lower CI from data). 

#-- Highlight sens and spec regions meeting requirements
AR_cont_spec <- Loop_spec; AR_cont_spec[which(AR_cont_spec<AR_min_spec)] <- NA
AR_cont_sens <- Loop_sens; AR_cont_sens[which(AR_cont_sens<AR_min_sens)] <- NA  

cont_AR_sens <- ifelse(is.na(AR_cont_sens) | is.na(AR_cont_spec), NA, AR_cont_sens)  #Leaves sensitivity values
cont_AR_spec <- ifelse(is.na(AR_cont_sens) | is.na(AR_cont_spec), NA, AR_cont_spec)  #Leaves specificity values

#-- TEa bands calculation (if required)
data_TE <- Loop_sens
for (i in 1:length(Imp_values)){
  for (b in 1:length(Bias_values)) {
    data_TE[b,i] <- abs(Bias_values[b])+1.96*Imp_values[i]
  }}

#-- Add stop function to provide a warning if the acceptable region is empty
if (all(is.na(cont_AR_sens))){
  stop('Contour plot acceptable region cannot be printed: the acceptable region for the specified minimum sensitivity and specificity is empty')
} else {

#-- Tiff code to save plot as a tiff file:  
tiff(filename="Cont_AR_temp.tiff", units="in", width=10, height=7, res=300) 

par(mfrow = c(1, 1))  #Specify one plot in the figure
par(mar=c(4.5,5,2,2)) #Specify space around the plot (bottom, left, top, right) (extra space needed for axis labels)

#-- Contour:  
contour(x=as.numeric(rownames(Loop_spec)), y=as.numeric(colnames(Loop_spec)), Loop_spec, levels=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95), col="midnightblue", xlab=expression(paste(Delta," Bias")), ylab=expression(paste(Delta," Imprecision (CV%)")), cex.axis=1.2, labcex=1.2, cex.lab=1.2, lwd=2, cex.axis=1.2, lty="dashed",drawlabels=FALSE)
contour(x=as.numeric(rownames(cont_AR_sens)), y=as.numeric(colnames(cont_AR_sens)), cont_AR_sens, add=TRUE, nlevels=10000, col="grey",drawlabels=FALSE, lwd=4) 
contour(x=as.numeric(rownames(cont_AR_spec)), y=as.numeric(colnames(cont_AR_spec)), cont_AR_spec, add=TRUE, nlevels=10000, col="grey",drawlabels=FALSE, lwd=4)

#-- Add grid lines to plot:
abline(h=seq(from=0,to=100,by=20), v=seq(from=-100,to=100,by=25),col="grey91",lty="solid")  #Use abline to specify exact placement of grid lines

#-- Re-add the contour lines so they appear over the AR:
contour(x=as.numeric(rownames(Loop_sens)), y=as.numeric(colnames(Loop_sens)), Loop_sens, levels=c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95), add=TRUE, col="deepskyblue",  labcex=1.2, cex.lab=1.2, lwd=2, drawlabels=TRUE, method="flattest")
contour(x=as.numeric(rownames(Loop_spec)), y=as.numeric(colnames(Loop_spec)), Loop_spec, levels=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),       add=TRUE, col="midnightblue", labcex=1.2, cex.lab=1.2, lwd=2, drawlabels=TRUE, method="flattest", cex.axis=1.2, lty="dashed")

#-- Add TE bands using data_TE (if required)
contour(x=as.numeric(rownames(data_TE)), y=as.numeric(colnames(data_TE)), data_TE, levels=c(10, 20, 30, 40, 50, 60), add=TRUE, col="magenta4", labcex=1, lwd=3, cex.axis=1, lty="dotted", drawlabels=TRUE, method="flattest")   #drawlabels=FALSE

#-- Legend (with TE band)
legend("topleft", inset=c(0.03,0.08), horiz=FALSE, c("Sensitivity", "Specificity", expression(paste("Acceptable region: sensitivity" >=0.85,"; specificity" >=0.90)), "TE% bands"), lty=c("solid","dashed","solid", "dotted"), col=c("deepskyblue","midnightblue","grey", "magenta4"), lwd=c(2,2,6,2))

#-- Legend (no TE bands)
#legend("topleft", inset=c(0.03,0.08), horiz=FALSE, c("Sensitivity", "Specificity", expression(paste("Acceptable region: sensitivity" >=0.85,"; specificity" >=0.90))), lty=c("solid","dashed","solid"), col=c("deepskyblue","midnightblue","grey"), lwd=c(2,2,6))

#-- Add manual contour line labels as desired (single example provided - user will need to add and update additional text labels as required dependant on the plot)
#text(-46,99,   "0.95",cex=1, pos=3, col="midnightblue")

#-- Save tiff file:
  dev.off()
}

#----- Acceptability Region information (Table data) -----

#-- Find the range of Bias allowed in the acceptatibility region when Imprecision ==0
temp <- cont_AR_sens[,"0"]              #Limit region to bias=0
temp <- temp[complete.cases(temp)]      #Remove NAs
AR_bias_atImp0 <- as.numeric(names(temp))
range(AR_bias_atImp0)   #Note this will return Inf to -Inf if the array is empty

#-- Find the range of Imprecision allowed when Bias ==0
temp <- cont_AR_sens["0",]
temp <- temp[complete.cases(temp)]
AR_Imp_atBias0 <- as.numeric(names(temp))
range(AR_Imp_atBias0)   






