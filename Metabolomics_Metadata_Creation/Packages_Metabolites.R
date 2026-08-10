#renv::activate()

library(here) #to be able to use this script given the project file
library(gtools) #to create combination from a vector
#library(psych) 
library(tidyverse) #for data wrangling
library(rstatix) # for different additional functions for correlation matrix handling and to extract eta-squared as effect size estimate for aov (anova) objects
library(future) #to enable smooth metadeconfoundR function
#library(lmerTest) # to perform the step function on a mixed model (lmer) object
library(lme4) #to model metabolomics data
#library(performance) #to check distribution of meta variables and (not implemented) metabolite abundance - for check_distribution() the random.forest package is needed.
#library(bestNormalize) #to normalize metabolite abundance data
library(mice) # to impute missing values for PCA analysis
#library(ggfortify) # to plot a nice score plot
library(factoextra) # for a nice loading plot
#library(formula.tools) # to restructure formulas in a loop
library(reshape2) # to create a plot for metadeconfoundR
library(rcompanion) # to calculate Man-Whitney-U test effect sizes
library(wesanderson) # for pretty colors in plots
#library(merTools) # to be able to scale columns within dplyr (strips attributes and thus retains structure after using base::scale())
library(metadeconfoundR) #see respective section
library(caret) #for machine learning apporaches such as random forest
library(pls) # function needed by caret to perform plsda
library(corpcor) # for outlier calculation

## Quality Control
library(factoextra) #for PCAs
library(metaproc) #for calibration
#library(metacalibr) # for metaproc
library(statTarget) # for QCRLSC batch correction
library(CovTools) # For Outlier filtering
library(emmeans) # To analyse the models using posthoc tests
library(ggbreak) # To enable breaks in a boxplot