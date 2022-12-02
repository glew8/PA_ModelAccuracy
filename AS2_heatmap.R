#script plots AS2 scores in heat map (ex. to graph color for unknown genes in sunburst plots)

##########################################################################################
#set parameters for analyses
wd = "C:/Users/Gina/Documents/PA/"
setwd(wd)

heat <- read.table("AS2_heatmap.csv", sep = ",", header = TRUE)

library(ggplot2)
library(reshape)
library(scico)

heat2 <- melt(heat)

ggplot(heat2, aes(variable, Condition, fill= value)) + 
  geom_tile() +
  scico::scale_fill_scico(palette = "romaO", limits=c(-.1,1))

pdf(file = "AS2Heatmap.pdf", width = 8.5, height = 11)
ggplot(heat2, aes(variable, Condition, fill= value)) + 
  geom_tile() +
  scico::scale_fill_scico(palette = "romaO", limits=c(-.1,1))
dev.off()