# Filippo Abbondanza, 25-02-2019
# Script adapted from Else Eising
# Script to get descriptions of a phenotypes
# This script is part of the Genome Wide Association Meta-Analysis of Quantitative Math traits - Analysis plan part 1 - version version 25 February 2019 

# --------------------------------------
# Manually adapt this section to match your phenotype data file
# Then run the remainder of the script
# --------------------------------------

phenotypeDirectory <- "~/Documents/PhD_data/GWAS_Math/NeuroDys/" # Directory with the data
outputDirectory <- "~/Documents/PhD_data/GWAS_Math/NeuroDys/" # Directory for the output
phenotypeData <- "~/Documents/PhD_data/GWAS_Math/NeuroDys/mzs_filippo.csv" # Name of the file
separator <- "," #separator in phenotype file. Put "," if working with .csv file
na_string <- "NA" #NA string in phenotype file
cohort = "MZS" # paste cohort name
measure_type = "scaled-age-normalized-Filippo" # describe the type of measure; e.g. "raw", "scaled-age-normalized"
 # Note: please run one type of measure at the time; run the script multiple times if different types of measures are available.

# fill in column name that contains age information
age_column <- "age"
# fill in column names of the phenotypes of interest as tehy are stated in the phenotype file
phenotype_columns <- c("readd_m_z_age_adj","remul_m_z_age_adj","ph_rechn1_Filippo")
# fill in type of phenotype data, in same order as columns
phenotype_descriptions <-  c("readd_m_z_age_adj","remul_m_z_age_adj","ph_rechn1_Filippo")
# fill in column name that contains gender information
gender <- "SEX"
# specify how males and females are coded
males <- "1"
females <- "2"


# --------------------------------------
# Open libraries and the multiplot function 
# --------------------------------------

if (!require("psych",character.only = TRUE))
    {install.packages("psych",dep=TRUE)}
if (!require("ggplot2",character.only = TRUE))
{install.packages("ggplot2",dep=TRUE)}
if (!require("dplyr",character.only = TRUE))
  {install.packages("dplyr",dep=TRUE)}
if (!require("corrplot",character.only = TRUE))
    {install.packages("corrplot",dep=TRUE)}
	
library(psych)
library(dplyr)
library(tidyverse)
library(corrplot)
library(grid)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
        ncol = cols, nrow = ceiling(numPlots/cols))
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# --------------------------------------
# Make and store overview files 
# --------------------------------------

setwd(phenotypeDirectory)
# pheno <- readxl::read_excel(phenotypeData)
pheno <- read.table(phenotypeData, sep=separator, header=TRUE, fill=T)
gender2 <- pheno[,which(colnames(pheno) == gender)]
columns <- c(age_column, phenotype_columns)
phenotypes <- c("Age", phenotype_descriptions)

# set working directory to location where phenotype info will be stored
setwd(outputDirectory)

# get descriptives for total cohort and for males and females separately
pheno_subset <- subset(pheno, select=columns)
colnames(pheno_subset) <- phenotypes
overview <- as.data.frame(describe(pheno_subset, quant=c(.25,.75), check=TRUE))
overview_males <- as.data.frame(describe(pheno_subset[which(gender2 == males),], quant=c(.25,.75), check=TRUE))
overview_females <- as.data.frame(describe(pheno_subset[which(gender2 == females),], quant=c(.25,.75), check=TRUE))

# combine descriptives for total cohort and for males and females separately
overview$subset <- "all"
overview_males$subset <- "males"
overview_females$subset <- "females"
overview_total <- rbind(overview, overview_males, overview_females)

# store the overview file
filename = paste(cohort, measure_type, "description_phenotypes.txt", sep="_")	
write.table(overview_total, filename, sep="\t", col.names=TRUE, row.names=TRUE)	

# --------------------------------------
# Make and store histograms of the phenotypes 
# --------------------------------------

# first of total cohort
plots <- list()
for(i in 1:length(columns)) {
	local({
	i <- i
	# get phenotype measures
	measures <- pheno_subset[,i]
	measures2 <- measures[!is.na(measures)]
	if((max(measures2)-min(measures2))<2) {bins <- 0.01} else if((max(measures2)-min(measures2))<10) {bins <- 0.1} else {bins <- 1}
	histo <- qplot(measures2, geom="histogram", binwidth=bins, xlab=phenotypes[i])
	# now store result in a list for the overview plot
	plots[[i]] <<- histo 
	})
}

# now with males and females having different colors
gender2 <- gsub(males, "males", gender2)
gender2 <- gsub(females, "females", gender2)

plots_males_and_females <- list()
for(i in 1:length(columns)) {
	local({
	i <- i
	# get phenotype measures
	measures <- pheno_subset[,i]
	measures2 <- measures[!is.na(measures)]
	if((max(measures2)-min(measures2))<2) {bins <- 0.01} else if((max(measures2)-min(measures2))<10) {bins <- 0.1} else {bins <- 1}
	histo <- ggplot(pheno, aes(x=measures, fill = gender2)) + geom_histogram(binwidth = bins, alpha=0.5,position="identity") + theme(legend.title=element_blank(), axis.text=element_text(size=5), axis.title=element_text(size=5)) +  labs(title=phenotypes[i], x=columns[i])
	# now store result in a list for the overview plot
	plots_males_and_females[[i]] <<- histo 
	})
}

# Now save all histograms
# Including an overview for all phenotypes together made using the multiplot function 
plotname = paste(cohort, measure_type, "histogram_phenotypes_all.pdf", sep="_")
pdf(plotname)
multiplot(plotlist=plots, cols=2)
plots
dev.off()

plotname = paste(cohort, measure_type, "histogram_phenotypes_males_and_females.pdf", sep="_")
pdf(plotname)
multiplot(plotlist=plots_males_and_females, cols=2)
plots_males_and_females
dev.off()

# --------------------------------------
# Make and store correlations between phenotypes and with age
# --------------------------------------

# calculate and store correlations between phenotypes
correlations <- cor(pheno_subset, use = "complete.obs", method="pearson")

plotname = paste(cohort, measure_type, "correlation_phenotypes.pdf", sep="_")
pdf(plotname)
corrplot(correlations,order="AOE",type="lower",tl.pos="tp", tl.col="black", tl.cex=0.7)
corrplot(correlations,add=TRUE, type="upper", method="number",order="AOE", col="black", diag=FALSE,tl.pos="n", cl.pos="n", number.cex = 0.7)
dev.off()

# plot and store relation between age and phenotypes
plots_age <- list()
age_info <- pheno[,which(colnames(pheno) == age_column)]

for(i in 1:length(columns)) {
	local({
	i <- i
	# get phenotype measures
	measures <- pheno_subset[,i]
	histo <- qplot(measures, age_info, geom="point", alpha=I(0.2), xlab=phenotypes[i], ylab="Age")
	plots_age[[i]] <<- histo 
	})
}

plotname = paste(cohort, measure_type, "correlation_phenotypes_with_age.pdf", sep="_")
pdf(plotname)
multiplot(plotlist=plots_age, cols=2)
plots_age
dev.off()


# QQ plot
plotname = paste(cohort, measure_type, "QQ_plot.pdf", sep="_")
plots_QQ <- list()
qq_info <- pheno_subset[,which(colnames(pheno_subset) != 'Age')]

for(i in 1:length(phenotype_columns)) {
  local({
    i <- i
    # get phenotype measures
    # measures <- qq_info[i]
    # print(colnames(measures))
    plot <- ggplot(qq_info, aes_string(sample = colnames(qq_info[i]) )) +
                             stat_qq() +
                             stat_qq_line() +
                             ggtitle(paste0("QQ plot for: ", colnames(qq_info[i])))
    plots_QQ[[i]] <<- plot
  })
}

pdf(plotname)
multiplot(plotlist=plots_QQ, cols=2)
dev.off()
message('Script finished')

