# PA_ModelAccuracy
Scripts to assess model accuracy using metatranscriptomic and transcriptomic data for Pseudomoans aeruginosa
Example input files are also included

AS2.R:
   1. Calculates AS2 or zscores for all samples in a model and plots sunburst graph
   2. Calculates AS2 and zscores after subsampling to keep input numbers consistent across different conditions
   3. Calculates zscores (penalties) for each individual replicate of a model

SunburstMetadata.csv: example metadata file. Can add additional columns to analyze data in different combinations.

SunburstCategories.txt: example input to make sunburst graphs using TIGRFAM categories

AnnotationTIGRFAM.txt: example input to make sunburst graphs using TIGRFAM categories

PA_individual_input_all.txt: example text file to calculate zscores for a set of individual replicates

################################

INSTALL the following packages:

use install.packages command in R unless otherwise noted

1. tidyverse (confirm installation of ggplot2, readr, tibble, purrr, dplyr, tidyr)
2. cowplot
3. zeallot
4. scico
5. ggsunburst: see below

ggsunburst install (This is what worked for me):

>install.packages(c("devtools", "reticulate", "reshape2", "rappdirs", "backports"))
>
>library(devtools)
>
>install_github("didacs/ggsunburst")
>
>library(reticulate)
>
>install.miniconda()
>
>py_install("six")
