## Morgens_fig_2D_analysis.R ##

#Red:Green ratios were calculated in FlowJo prior to data export

#load dependencies
library(ggplot2)
library(plyr)
library(dplyr)

#Load functions used in this processing script. A brief description of each function is provided below:

#The import function takes the name of a csv file in the active directory and imports it, assigning them variables based on the original filename from their export from FlowJo. This returns a dataframe containing all flow cytometry data with appropriate variable names.
import <- function(f){
  data <- read.csv(f)
  split_fnames <- strsplit(data$file, '_')
  split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
  data$replicate <- as.numeric(substr(strsplit(f,'_')[[1]][3],4,4))
  data$strain <- split_fnames$X4
  data$reporter <- split_fnames$X5
  data$drug <- split_fnames$X6
  data
}

#Normalize takes a dataframe containing data to be normalized (d), numeric replicate number (rep), dataframe containing averaged control data to be normalized to (ctrl), and a list consisting of the unique reporters used in this experiment (reporters). This subsets the d dataframe into a dataframe containing only one type of reporter and replicate, normalizes it to the control of the same reporter in the ctrl dataframe and binds all reporter normalized dataframes into a single dataframe. It returns a subset of the dataframe d containing only data points of replicate rep, normalized to the control of their respective reporter.
normalize <- function(d,rep,ctrl,reporters){
  reps <- subset(data,replicate == rep)
  temprep <- lapply(reporters,function(x) {
    report <- subset(reps, reporter==x)
    report$normalized <- report$`mCherry.A.FITC.A`/ctrl$mean_ratio[which(ctrl$reporter == x & ctrl$replicate==rep)]
    report}) %>% bind_rows()
  temprep
}

#Process takes a dataframe consisting of flow cytometry data (data) and the numeric number of replicates (r) and summarizes the data. A control dataframe (ctrl) is created from the data dataframe, which is then used to normalize all observations to their respective controls. The function returns a list of two dataframes. The first element of this list, statdata,  is the summarized mean of each replicate for each strain, reporter and drug treatment. The second element of this list, plotdata, takes the mean of these replicate means and their standard deviations. 
process <- function(data,r){
  reporters <- unique(data$reporter)
  ctrl <- ddply(.data = subset(data,drug=='mock' & strain == 'WT'),.variables = c('reporter','replicate'), summarize,
                mean_ratio = mean(`mCherry.A.FITC.A`),
                ratio_sd = sd(`mCherry.A.FITC.A`))
  normdata <- lapply(1:r,normalize,d=data,ctrl=ctrl,reporters=reporters) %>% bind_rows()
  statdata <- ddply(.data=normdata, .variables = c('drug','strain','reporter','replicate'), summarize,
                    mean_normalized = mean(normalized))
  plotdata <- ddply(.data=statdata, .variables = c('drug','strain','reporter'),summarize,
                    mean = mean(mean_normalized),
                    sd = sd(mean_normalized))
  list(statdata,plotdata)
}

#Change the following path to point to your clone of the rep1 data
setwd("/Users/akane/R_testing/Morgens_et_al_2019/figure_2/")

#Creates a list of csv files to be processed in this script.
files = list.files(pattern="*2d_rep*")

#Imports the files and gives appropriate variable names, as described in the import function. Binds all fo these files into one dataframe.
data <- lapply(files,import) %>% bind_rows()

#Processes the data as described in the process function. Creates a summary of the data broken down by replicate or as mean and standard deviation (elements 1 and 2 of processed object)
processed <- process(data,length(unique(data$replicate)))

#Takes elements 1 and 2 of processed object and makes them their own dataframe objects for ease of plotting and statistical analysis.
statdata <- as.data.frame(processed[1])
plotdata <- as.data.frame(processed[2])

#Sets levels to order the reporters on the plot.
plotdata$reporter <- factor(plotdata$reporter, levels = c('GFP2ARFPSec61','GFP2ARFPSec61BFP'))

#Creates a ggplot as a bar chart with standard deviations as error bars. Proportions can be changed using coor_fixed
ggplot(plotdata, aes(x = factor(drug,level=c('mock','DHQZ36.1')), y = mean, fill = factor(strain,level=c('WT','ASNA1KO')))) +
  facet_wrap(~reporter, ncol = 2)+
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                width=.2,                    
                position=position_dodge(.9))+
  scale_y_continuous(expand=c(0,0), name ='Normalized Mean RFP-Sec61 Transmembrane Domain:GFP') +
  scale_x_discrete(name = '')+
  theme(panel.background = element_rect(fill=NA, color = NA),
        panel.grid =  element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.line = element_line(color='black'))

#ggsave('~/path/to/Fig/2B/plot.pdf', device = cairo_pdf(width = 10, height = 5))
# Further aesthetic changes were made in Adobe Illustrator.