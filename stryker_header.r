#script to load necessary libraries and scripts to run the stryker analysis

#data manipulation
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(reshape2)))

#fitting distributions
suppressMessages(suppressWarnings(library(fitdistrplus)))

#bayesian
#suppressMessages(suppressWarnings(library(runjags)))
suppressMessages(suppressWarnings(library(rstan)))

#plotting
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(gridExtra))) #plot arrangement
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(colorspace)))

#load custom functions
source("lib/general_utility_functions.R")
source("lib/stryker_functions.R")          #Stryker analysis from R script file
source("lib/custom_plotting_functions.R")  #custom plotting specs

#set the default ggplot2 theme to custom specifications
theme_set(custom_theme) 
