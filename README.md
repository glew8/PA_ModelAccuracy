# PA_ModelAccuracy
Scripts to assess model accuracy using metatranscriptomic and transcriptomic data for Pseudomoans aeruginosa
Example input files are also included

AS2.R:
   1. Calculates AS2 or zscores for all samples in a model and plots sunburst graph
   2. Calculates AS2 and zscores after subsampling to keep input numbers consistent across different conditions
   3. Calculates zscores (penalties) for each individual replicate of a model
   
AS2_heatmap.R: plots AS2 scores in heat map (ex. to graph color for unknown genes in sunburst plots)

Input files
   1. SunburstMetadata.csv: example metadata file. Can add additional columns to analyze data in different combinations.
   2. SunburstCategories.txt: example input to make sunburst graphs using TIGRFAM categories
   3. SunburstCategories_withUnknowns.txt: example input to make sunburst graphs using TIGRFAM categories, includes genes of unknown function
   4. AnnotationTIGRFAM.txt: example input to make sunburst graphs using TIGRFAM categories
   5. AnnotationTIGRFAM_withUnknowns.txt: example input to make sunburst graphs using TIGRFAM categories, includes genes of unknown function
   6. PA_individual_input_all.txt: example text file to calculate zscores for a set of individual replicates
   7. AS2_heatmap.csv: example input to make heatmap

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


################################
scripts modified from
   1. D. M. Cornforth, F. L. Diggle, J. A. Melvin, J. M. Bomberger, M. Whiteley, Quantitative framework for model evaluation in microbiology research using Pseudomonas aeruginosa and cystic fibrosis infection as a test case. mBio 11, e03042-03019 (2020). https://journals.asm.org/doi/10.1128/mBio.03042-19
   2. G. R. Lewin, K. S. Stocke, R. J. Lamont, M. Whiteley, A quantitative framework reveals traditional laboratory growth is a highly accurate model of human oral infection. Proc Natl Acad Sci U S A 119, e2116637119 (2022). https://www.pnas.org/doi/10.1073/pnas.2116637119; https://github.com/glew8/Pgingivalis_Metatranscriptome_Analyses
