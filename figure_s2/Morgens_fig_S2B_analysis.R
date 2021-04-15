## Morgens_fig_S2B_analysis.R ##

#Red:Green ratios were calculated in FlowJo prior to data export

#load dependencies
library(readr)
library(ggplot2)
library(plyr)
library(dplyr)

#Load functions used in this processing script. A brief description of each function is provided below:

#The import function takes the name of a csv file in the active directory and imports it, assigning them variables based on the original filename from their export from FlowJo. This returns a dataframe containing all flow cytometry data with appropriate variable names.
import <- function(f){
  data <- read_csv(f)
  split_fnames <- strsplit(data$file, '_')
  split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
  data$replicate <- as.numeric(substr(strsplit(f,'_')[[1]][3],4,4))
  data$strain <- split_fnames$X4
  data$drug <- split_fnames$X6
  data
}

#Normalize takes a dataframe containing data to be normalized (d), numeric replicate number (rep), and a dataframe containing averaged control data to be normalized to (ctrl). This subsets the d dataframe into a dataframe containing only one replicate and normalizes it to the control of the same replicate. It returns a subset of the dataframe d containing only data points of replicate rep, normalized to the control of their respective replicate.
normalize <- function(d,rep,ctrl){
  reps <- subset(data,replicate == rep)
  reps$normalized <- reps$`mCherry.A.FITC.A`/ctrl$mean_ratio[which(ctrl$replicate == rep)]
  reps
}

#Process takes a dataframe consisting of flow cytometry data (data) and the numeric number of replicates (r) and summarizes the data. A control dataframe (ctrl) is created from the data dataframe, which is then used to normalize all observations to their respective controls. The function returns a list of two dataframes. The first element of this list, statdata,  is the summarized mean of each replicate for each strain and drug treatment. The second element of this list, plotdata, takes the mean of these replicate means and their standard deviations. 
process <- function(data,r){
  ctrl <- ddply(.data = subset(data,drug=='mock' & strain == 'WT'),.variables = 'replicate', summarize,
                mean_ratio = mean(`mCherry.A.FITC.A`),
                ratio_sd = sd(`mCherry.A.FITC.A`))
  normdata <- lapply(1:r,normalize,d=data,ctrl=ctrl) %>% bind_rows()
  statdata <- ddply(.data=normdata, .variables = c('drug','strain','replicate'), summarize,
                    mean_normalized = mean(normalized))
  plotdata <- ddply(.data=statdata, .variables = c('drug','strain'),summarize,
                    mean = mean(mean_normalized),
                    sd = sd(mean_normalized))
  list(statdata,plotdata)
}

#Change the following path to point to your clone of the rep1 data
setwd("/Users/akane/R_testing/Morgens_et_al_2019/figure_s2/")

#Creates a list of csv files to be processed in this script.
files = list.files(pattern="*s2b_rep*")

#Imports the files and gives appropriate variable names, as described in the import function. Binds all fo these files into one dataframe.
data <- lapply(files,import) %>% bind_rows()

#Processes the data as described in the process function. Creates a summary of the data broken down by replicate or as mean and standard deviation (elements 1 and 2 of processed object)
processed <- process(data,length(unique(data$replicate)))

#Takes elements 1 and 2 of processed object and makes them their own dataframe objects for ease of plotting and statistical analysis.
statdata <- as.data.frame(processed[1])
plotdata <- as.data.frame(processed[2])

#Creates a ggplot as a bar chart with standard deviations as error bars. Proportions can be changed using coor_fixed
ggplot(plotdata, aes(x = factor(drug,level=c('mock','Retro2','DHQZ36.1','DHQZ5','Retro1','hydroR1')), y = mean, fill = factor(strain,level=c('WT','ASNA1KO')))) +
  coord_fixed(ratio = 3.5)+
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                width=.2,                    
                position=position_dodge(.9))+
  scale_y_continuous(expand=c(0,0), name ='Normalized Mean RFP-Stx5 Transmembrane Domain:GFP') +
  scale_x_discrete(name = '')+
  theme(panel.background = element_rect(fill=NA, color = NA),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.line = element_line(color='black'))

#ggsave('~/path/to/Fig/S2B/plot.pdf', device = cairo_pdf(width = 10, height = 5))
# Further aesthetic changes were made in Adobe Illustrator.