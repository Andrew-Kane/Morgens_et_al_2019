## Morgens_fig_S4-2E_analysis.R ##

#Red:Green ratios were calculated in FlowJo prior to data export

#load dependencies
library(readr)
library(plyr)
library(dplyr)
library(drc)

#Load functions used in this processing script. A brief description of each function is provided below:

#The import function takes the name of a csv file in the active directory and imports it, assigning them variables based on the original filename from their export from FlowJo. This returns a dataframe containing all flow cytometry data with appropriate variable names.
import <- function(f){
  data <- read_csv(f)
  split_fnames <- strsplit(data$file, '_')
  split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
  data$replicate <- as.numeric(substr(strsplit(f,'_')[[1]][3],4,4))
  data$strain <- split_fnames$X4
  data$reporter <- split_fnames$X5
  data$drug <- split_fnames$X6
  data$concentration <- as.character(split_fnames$X7)
  data$concentration <- as.numeric(data$concentration)
  data
}

#Normalize takes a dataframe containing data to be normalized (d), numeric replicate number (rep), and a dataframe containing averaged control data to be normalized to (ctrl). This subsets the d dataframe into a dataframe containing only one replicate and normalizes it to the control of the same replicate. It returns a subset of the dataframe d containing only data points of replicate rep, normalized to the control of their respective replicate.
normalize <- function(d,rep,ctrl){
  reps <- subset(data,replicate == rep)
  reps$normalized <- reps$`mCherry.A.FITC.A`/ctrl$mean_ratio[which(ctrl$replicate == rep)]
  reps
}

#Process takes a dataframe consisting of flow cytometry data (data) and the numeric number of replicates (r) and summarizes the data. A control dataframe (ctrl) is created from the data dataframe, which is then used to normalize all observations to their respective controls. The function returns the returndata dataframe. This is the summarized mean of each replicate for each strain, reporter, drug treatment, and concentration. 
process <- function(data,r){
  strains <- unique(data$strains)
  ctrl <- ddply(.data = subset(data,concentration==0 & strain == 'K562WT'),.variables = 'replicate', summarize,
                mean_ratio = mean(`mCherry.A.FITC.A`),
                ratio_sd = sd(`mCherry.A.FITC.A`))
  normdata <- lapply(1:r,normalize,d=data,ctrl=ctrl) %>% bind_rows()
  returndata <- ddply(.data=normdata, .variables = c('drug','strain','reporter','concentration','replicate'), summarize,
                      mean_normalized = mean(normalized))
  returndata$logconc <- log10(returndata$concentration)
  returndata
}

#Fitter takes a dataframe (vals) containing summarized and normalized flow cytometry data and a list of strains (s) to be individually fitted. It uses the dose response curve package to fit a curve to a subset of the data corresponding to the strain indicated by s.
fitter <- function(vals,s){
  fitted <- drm(formula = mean_normalized ~ concentration, data = subset(vals,strain == s), fct = LL.4(names = c("Slope", "Lower", "Upper", "EC50")))
  fitted
}

#Change the following path to point to your clone of the rep1 data
setwd("/Users/akane/R_testing/Morgens_et_al_2019/figure_s4-2/")

#Creates a list of csv files to be processed in this script.
files = list.files(pattern="*s4-2e_rep*")

#Imports the files and gives appropriate variable names, as described in the import function. Binds all fo these files into one dataframe.
data <- lapply(files,import) %>% bind_rows()

#Processes the data as described in the process function. Creates a summary of the data broken down by replicate or as mean and standard deviation (elements 1 and 2 of processed object)
processed <- process(data,length(unique(data$replicate)))

#Creates a list of all strains in the processed dataset.
strains <- unique(processed$strain)

#Creates a list of drc fits for each strain in the processed dataset.
fits <- lapply(strains,fitter,vals=processed)

#Plots the WT fit.
index <- which(strains=='K562WT')
plot(fits[[index]], type = c("bar"),ylim=c(0.35,1.05))
plot(fits[[index]], add=TRUE, type = c("average"),col="black",pch=21)

#Plots the A149V fit.
index <- which(strains=='K562A149V')
plot(fits[[index]], add=TRUE, type = c("bar"), col="red")
plot(fits[[index]], add=TRUE, type = c("average"),col="red",pch=22)

#Plots the A149A fit.
index <- which(strains=='K562A149A')
plot(fits[[index]], add=TRUE, type = c("bar"), col="blue")
plot(fits[[index]], add=TRUE, type = c("average"),col="blue",pch=22)

#ggsave('~/path/to/Fig/S4-2E/plot.pdf', device = cairo_pdf(width = 10, height = 5))
# Further aesthetic changes were made in Adobe Illustrator.