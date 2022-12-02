#script calculates AS2 to compare a model to human transcriptome and then graphs output using sunburst plots
#   1. Calculates AS2 or zscores for all samples in a model and plots sunburst graph
#   2. Calculates AS2 and zscores after subsampling to keep input numbers consistent across different conditions
#   3. Calculates zscores (penalties) for each individual replicate of a model

##########################################################################################
#set parameters for all analyses
#also need to run functions at the end of the script

#set working directory
wd = "C:/Users/Gina/Documents/PA/"
setwd(wd)

library(tidyverse)
library(cowplot)
library(zeallot)
library(ggsunburst)
library(scico)

#input data. 
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

#load TIGR annotations
TIGR_ann <- read.csv("AnnotationTIGRFAM.txt", sep = "\t")
SunburstInput <- sunburst_data("SunburstCategories.txt", sep = "\t", type = "node_parent", node_attributes = "color")

#set metadata: change list names and filter condition to match datasets you want to compare
metadata_file <- read_csv(file = "SunburstMetadata.csv")

#setsession
invitro_dataset <- "CP" #set based on model set you want to analyze
outputversion <- "run1"

human_list <- metadata_file %>% filter(type1 == "human") %>% .$filename %>% str_replace_all("-", "_")
invitro_list <- metadata_file %>% filter(type1 == invitro_dataset) %>% .$filename %>% str_replace_all("-", "_") #adjust filter (currently "type1") to match correct column of metadata_file

#########################################################################################
#1. calculate AS2s and plot sunbursts for groups of samples, based on parameters set above (could also work for individual samples depending on parameters)
#run functions below first

AS_2 <- draw_sunburst_STANDARD_PA(counts_normalized, TIGR_ann, SunburstInput, human_list, invitro_list, 2, TRUE)

#output AS2
penalty_output <- AS_2[[2]]
AS_2_score <- penalty_output %>% mutate(score_2 = ifelse((penalty > 2 | penalty < -2), 0, 1)) %>% mutate(AS2 = mean(score_2))
write.table(AS_2_score, paste(outputversion, invitro_dataset, "AS2penalty.txt", sep = "."), sep = "\t", row.names = FALSE) 
head(AS_2_score)

#plot sunburst
AS_2[[1]] 
pdf(file = paste(outputversion, invitro_dataset, "AS2sunburst.pdf", sep = "."), width = 8.5, height = 11)
AS_2[[1]]
dev.off()
#we finish editing this sunburst image manually in Adobe Illustrator



#######################################################################################
#2. to calculate AS2 from subsampling X replicates at once (set below), doing all possible combinations.
#run functions below first

#set session
reps <- 2 #number of replicates to run at once if subsampling
outputversion <- paste("run1.sub", reps, sep = "")

#JUST OUTPUT AS2: outputs average AS2 value across iterations
c(1:choose(length(invitro_list), reps)) %>% map(function(.x) {dummy1 <- score_target_vs_model_PA(counts_normalized, human_list, combn(invitro_list, reps)[,.x]) %>% .$penalty
dummy2 <- abs(dummy1) < 2
dummy3 <- mean(dummy2)
return(dummy3)}) %>% purrr::reduce(c) %>% mean

#outputs table of each iteration with avg zscores across iterations. Can edit then use as input into sunburst to directly graph this version of the analysis 
zscores <- c(1:choose(length(invitro_list), reps)) %>% map(function(.x) {dummy1 <- score_target_vs_model_PA(counts_normalized, human_list, combn(invitro_list, reps)[,.x]) %>% .$penalty
return(dummy1)})
zscores2 <- as.data.frame(do.call(cbind, zscores))
zscores3 <- cbind(counts_normalized$locus_tag, zscores2)
colnames(zscores3)[1] <- "locus_tag"
zscores4 <- zscores3 %>% mutate(xmean = rowMeans(across(starts_with("V"))))
write.table(zscores4, paste(outputversion, invitro_dataset, "zscores.txt", sep = "."), sep = "\t", row.names = FALSE) 


#########################################################################################
#3. to calculate zscores for individual replicates
#run functions below first

datasets <- as.character(read.csv("PA_individual_input_all.txt", header = FALSE)$V1) #text file with individual sample names listed in single column
#also must have a column called "ind_condition" in metadata file that lists the name of each replicate

AS2_bygene <- data.frame(locus_tag=c(0))

for(i in 1:length(datasets)) {
  invitro_list <- metadata_file %>% filter(ind_condition == datasets[i]) %>% .$filename %>% str_replace_all("-", "_")
  dummy1 <- score_target_vs_model_PA(counts_normalized, human_list, invitro_list) #%>% .$penalty
  colnames(dummy1) <- c("locus_tag", invitro_list)
  AS2_bygene <- merge(AS2_bygene, dummy1, all=T)
}


AS2_bygene <- tail(AS2_bygene, n=-1)
write.table(AS2_bygene, paste(version, "AS2_bygene_individual.txt", sep = "."), sep = "\t", row.names = FALSE) #input into excel and calculate AS2 from penalties


#########################################################################################
#functions used above, run each
score_target_vs_model_PA <- function(counts_DF, Target_samplenames, Model_samplenames) # given a counts DF, lists of target samples, and lists of model samples, will return a score for each gene
{
  
  all_tested_samples_list <- names(counts_DF)[-1]
  names(counts_DF)
  
  conditions_list <- all_tested_samples_list
  dim(counts_DF %>% select(locus_tag, all_tested_samples_list) %>% column_to_rownames("locus_tag"))
  
  df1 <- counts_DF
  setdiff(Target_samplenames, names(df1))
  DF_target <- df1  %>% select(locus_tag, Target_samplenames) %>% gather(key=sample_name, value=expression, -locus_tag) %>% group_by(locus_tag) %>% dplyr::summarize(target_mean = mean(expression), target_SD = sd(expression))
  
  df_allzscores <- DF_target %>% inner_join(df1 %>% select(locus_tag, Model_samplenames)) %>% mutate_at(.vars=vars(-locus_tag, -target_mean, -target_SD), .funs=funs((.-target_mean)/target_SD))
  mean_modelZscoreDF <- df_allzscores %>% select(locus_tag, Model_samplenames)  %>% transmute(locus_tag, penalty_temp = pmap_dbl(.[c(-1)], function(...)  mean(c(...)))) %>% mutate(penalty = round(penalty_temp, digits = 4)) %>% select(locus_tag, penalty) 
  
  # below should get rid of NAs
  
  if(length(which(is.na(mean_modelZscoreDF$penalty)))>0)
  {
    mean_modelZscoreDF[which(is.na(mean_modelZscoreDF$penalty)),]$penalty <- 0
  }  
  final_DF <- mean_modelZscoreDF  
  
  return(final_DF)
}

draw_sunburst_STANDARD_PA <- function(PA_df_counts, TIGR_PA, sb, target_list, model_list, maxlim, node_labs)
{
  
  TIGR_PA <- TIGR_PA  %>% mutate(sub_role = if_else(sub_role == "Other", paste0(sub_role, "@", main_role), sub_role))
  PA_DE <- score_target_vs_model_PA(PA_df_counts, target_list, model_list)
  
  subrole_penalty <- TIGR_PA %>% left_join(PA_DE)  %>% group_by(sub_role) %>%  select(locus_tag, sub_role, penalty) %>% distinct %>% dplyr::summarize(penalty = mean(abs(penalty)<maxlim, na.rm=TRUE))
  main_role_penalty <- TIGR_PA %>% left_join(PA_DE) %>% group_by(main_role) %>% select(locus_tag, main_role, penalty) %>% distinct %>% dplyr::summarize(penalty = mean(abs(penalty)<maxlim, na.rm=TRUE))
  Meta_penalty <- TIGR_PA %>% left_join(PA_DE) %>% group_by(Meta) %>% select(locus_tag, Meta, penalty) %>% distinct %>% dplyr::summarize(penalty = mean(abs(penalty)<maxlim, na.rm=TRUE))
  subrole_penalty %>% arrange(penalty)
  sum(abs(PA_DE$penalty) < 2)
  TIGR_PA %>% left_join(PA_DE)  %>% arrange(abs(penalty)) #%>% filter(sub_role == "Sulfur metabolism")
  
  subrole_penalty %>% filter()
  
  names(Meta_penalty)[1] <- "role"
  names(subrole_penalty)[1] <- "role"
  names(main_role_penalty)[1] <- "role"
  
  all_penalty_levels <- rbind(Meta_penalty, main_role_penalty, subrole_penalty)
  all_penalty_levels$role <- all_penalty_levels$role %>% as.character
  
  sb$rects$name <- sb$rects$name %>% as.character
  temp1 <- sb$rects %>% left_join(all_penalty_levels, by = c("name"="role"))
  
  sb$rects$color <- temp1$penalty
  sb$rects$leaf <- TRUE
  
  
  tempinnerplot <- sunburst(sb, node_labels = node_labs, node_labels.min = .2, node_labels.size = .8, leaf_labels = FALSE, rects.fill.aes = "color", rects.size  = 0.1) + scico::scale_fill_scico(palette = "romaO", limits=c(-.1,1))
  d=data.frame(x1=.5, x2=1, y1=.5, y2=1)
  backplot_1 <- ggplot() + geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill= sum(PA_DE$penalty %>% abs < maxlim, na.rm = TRUE)/length(PA_DE$penalty))) + scico::scale_fill_scico(palette = "romaO", limits=c(-.1,1)) + theme(legend.position="none")
  prow <- plot_grid(tempinnerplot,
                    backplot_1,
                    align = 'v',
                    labels = c("inner", "middle"),
                    ncol = 1
  )
  
  plot_grid(prow)
  p <- plot_grid(prow)
  
  return(list(p, PA_DE))
}
