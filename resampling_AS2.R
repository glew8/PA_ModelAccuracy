
#Script calculates AS2 across human samples to set the AS2 baseline of a microbe in situ
#Approach repeatedly calculates AS2, using the majority of human transcriptomes as the gold-standard human datasets, and using a few random transcriptomes as the comparison model to see how similar the datasets are to each other

#set working directory
wd = "C:/Users/Gina/Documents/PA/"
setwd(wd)

library(tidyverse)
library(cowplot)
library(zeallot)
library(ggsunburst)

#input data
inputversion <- "input1"

#I input data that I have already normalized using vst or rlog. Otherwise, you can just load the raw data and there is code to normalize below. You need to make sure that you always normalize all of your data together (even if you are only analyzing a small subset) because the normalization varies depending on the entire dataset
#load normalized counts file here (e.g., VST, rlog). Want file where columns are samples and rows are genes.
counts_normalized <- read.table(paste(inputversion, "vst.dds2.csv", sep = "."), sep = ",", header = TRUE)
names(counts_normalized)[1] <- "locus_tag"
rownames(counts_normalized) <- NULL

#load raw counts file here. Want file where columns are samples and rows are genes.
#final_raw <- read.table("counts.txt", sep = "\t", header = TRUE)
#names(final_raw)[1] <- "locus_tag"
#will then need to normalize using rlog or VST in DESeq2
# Load raw counts file here. Want file where columns are samples and rows are genes, then normalize with DESeq2 with VST
# final_raw %>% select(locus_tag, all_tested_samples_list) %>% column_to_rownames("locus_tag")
# dds = DESeqDataSetFromMatrix(countData = final_raw %>% select(locus_tag, all_tested_samples_list) %>% column_to_rownames("locus_tag"), colData =    data.frame(condition=conditions_list), design = ~condition)
# dds2 <- varianceStabilizingTransformation(dds, blind=TRUE)

#set metadata: change list names and filter condition to match datasets you want to compare. Here you want the human samples to be clearly defined.
metadata_file <- read_csv(file = "SunburstMetadata.csv")

#define sample lists
human_list <- metadata_file %>% filter(type1 == "human") %>% .$filename %>% str_replace_all("-", "_")

conditions_list <- names(counts_normalized)[-1]

#setsession
iterations <- 1200 #set number of times to resample
leaveout <- 3 #set number of samples to use as test model
output_version <- "run1"

#################################################################################################################

#run multiple iterations of AS2 using X samples as the gold-standard and Y samples as the "model." For example, for my 12 human samples, I run 10 as the human "gold-standard" and 2 as the test across 100 iterations.
#define functions below
counts_normalized_double <- c(1:iterations) %>% map(function(x) {dummy1 <- model_selfvalidation_leaveout(self_validation, leaveout)
dummy1$penalty <- abs(dummy1$penalty)
names(dummy1)[2] <- paste0("penalty_", x)
return(dummy1)})  %>% purrr::reduce(left_join)

write.table(counts_normalized_double, file = paste(output_version, "_resamplingAS2_iterations", iterations, "_leaveout", leaveout, ".txt", sep = ""), sep ="\t")
#output is penalties (z-scores). In R or excel, calculate AS2 for each iteration then average iterations

#################################################################################################################
#functions to calculate AS2
score_target_vs_model <- function(Target_samplenames, Model_samplenames)
{
  DF_target <- counts_normalized  %>% select(locus_tag, Target_samplenames) %>% gather(key=sample_name, value=expression, -locus_tag) %>% group_by(locus_tag) %>% dplyr::summarize(target_mean = mean(expression), target_SD = sd(expression))
  
  df_allzscores <- DF_target %>% inner_join(counts_normalized %>% select(locus_tag, Model_samplenames)) %>% mutate_at(.vars=vars(-locus_tag, -target_mean, -target_SD), .funs=funs((.-target_mean)/target_SD))
  
  mean_modelZscoreDF <- df_allzscores %>% select(locus_tag, Model_samplenames)  %>% transmute(locus_tag, penalty_temp = pmap_dbl(.[c(-1)], function(...)  mean(c(...)))) %>% mutate(penalty = round(penalty_temp, digits = 4)) %>% select(locus_tag, penalty) 
  
  # below should get rid of NAs
  
  if(length(which(is.na(mean_modelZscoreDF$penalty)))>0)
  {
    mean_modelZscoreDF[which(is.na(mean_modelZscoreDF$penalty)),]$penalty <- 0
  }  
  final_DF <- mean_modelZscoreDF  
  
  return(final_DF)
}

model_selfvalidation_leaveout <- function(sample_DF1, left_out_number)
{
  test_sputum <- sample(1:length(human_list), left_out_number, replace = FALSE)
  train_sputum <- setdiff(1:length(human_list), test_sputum)
  
  dummy1 <- score_target_vs_model(human_list[train_sputum], human_list[test_sputum])
  return(dummy1)
}