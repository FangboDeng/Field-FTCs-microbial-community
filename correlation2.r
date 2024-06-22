##bacterial communites analysis as an example

# install.packages("dplyr")
library(dplyr)

# read module data
module_data <- read.csv("ModuleEigengene.csv", row.names = 1)

# Soil vairbles data file names
prop_files <- c("Propertities.csv", "Alpha.csv", "Potential functions.csv")

results <- data.frame()

for (file in prop_files) {
  prop_data <- read.csv(file, row.names = 1)

  # ensure sample name match
  if (!all(rownames(prop_data) == rownames(module_data))) {
    stop("sample name not match: ", file)
  }

  # calculate correlation for each module
  for (module in colnames(module_data)[1:4]) {  # only calculate main modules #0-3
    for (prop in colnames(prop_data)) {
      corr_test <- cor.test(module_data[,module], prop_data[,prop], method = "pearson")

      results <- bind_rows(results, data.frame(File = file, Module = module, Prop = prop,
                                               Correlation = corr_test$estimate, 
                                               PValue = corr_test$p.value))
    }
  }
}

print(results)


# read keystone OTUs data
OTU_data <- read.csv("keystone_abundance.csv", row.names = 1)

# Soil propertities data file
prop_files <- c("Propertities.csv")

results2 <- data.frame()

for (file in prop_files) {
  prop_data <- read.csv(file, row.names = 1)

  # ensure sample name match
  if (!all(rownames(prop_data) == rownames(OTU_data))) {
    stop("sample name not match: ", file)
  }

  # calculate correlation for each OTU
  for (module in colnames(OTU_data)) {
    for (prop in colnames(prop_data)) {
      corr_test <- cor.test(OTU_data[,module], prop_data[,prop], method = "pearson")

      results2 <- bind_rows(results2, data.frame(File = file, Module = module, Prop = prop,
                                               Correlation = corr_test$estimate, 
                                               PValue = corr_test$p.value))
    }
  }
}

print(results2)

