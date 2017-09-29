# Run the analysis
# Crystal Tu 20170929

# Variation partitioning
source("~/Desktop/fishsize_varpart/script/variation_partitioning.r")
# Combine all time lag & LM on explained fraction of fishing/temperature v.s. Lifehist
source("~/Desktop/fishsize_varpart/script/meta_variation_partitioning.R")

# Univariate SBIs
source("~/Desktop/fishsize_varpart/script/univariateSBI.R")
source("~/Desktop/fishsize_varpart/script/meta_uniSBIs.R")

# Summary plot + ANOVA test
source("~/Desktop/fishsize_varpart/script/summary.R")

# Comparison plot + binomial test
source("~/Desktop/fishsize_varpart/script/comparison.R")

