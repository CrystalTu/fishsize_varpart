# Run the analysis
# Crystal Tu 20170929
# 1st revise 20180213

setwd("~/Desktop/fishsize_varpart/")
dir.create("~/Desktop/fishsize_varpart/output")

# Prepare temperature and fishing mortality data 
#source("~/Desktop/fishsize_varpart/script/exploitation2F.r") #fishing mortality conversion
#source("~/Desktop/fishsize_varpart/data/Alaska/envi/presetAlaska.r") #EBS, GOA temperature
#source("~/Desktop/fishsize_varpart/data/NorthSea/Oceanography/temperature.r") #IBTS Q1 temperature

# Load dataset for analysis
setwd("~/Desktop/fishsize_varpart/")
source("~/Desktop/fishsize_varpart/data/preset.R") #West US + Alaska length 
source("~/Desktop/fishsize_varpart/data/NorthSea/IBTSsurvey/preset.R") # get length frequency for North Sea
source("~/Desktop/fishsize_varpart/data/NorthSea/stockassessment/stockassessment.R") # get stockassessment result table

# Variation partitioning (done)
setwd("~/Desktop/fishsize_varpart/")
source("~/Desktop/fishsize_varpart/script/variation_partitioning_revise.r")
source("~/Desktop/fishsize_varpart/script/variation_partitioning_revise1yr.r")
source("~/Desktop/fishsize_varpart/script/variation_partitioning_revise3yr.r")

# Combine all time lag & LM on explained fraction of fishing/temperature v.s. Lifehist
source("~/Desktop/fishsize_varpart/script/meta_variation_partitioning_revise.R")

# Univariate SBIs
source("~/Desktop/fishsize_varpart/script/univariateSBI_revise.R")
source("~/Desktop/fishsize_varpart/script/univariateSBI_revise1yr.R")
source("~/Desktop/fishsize_varpart/script/univariateSBI_revise3yr.R")
source("~/Desktop/fishsize_varpart/script/meta_uniSBIs_revision.R")

# Summary plot + ANOVA test
source("~/Desktop/fishsize_varpart/script/summary_revise.R")

# Comparison plot + binomial test
source("~/Desktop/fishsize_varpart/script/comparison_revise.R")

# Supplementary: bottom/surface temperature plot and correlation analysis
source("~/Desktop/fishsize_varpart/script/envi_temp.R")