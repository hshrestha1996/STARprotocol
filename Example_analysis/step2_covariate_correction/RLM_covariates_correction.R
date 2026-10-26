# Step 2 - Covariate correction (age, sex, PMI, batch) with robust linear models.
# Run from inside the step2_covariate_correction folder.

# Load necessary libraries
library(WGCNA)      # corPvalueStudent()
library(ggplot2)
library(reshape2)

# Read CSV file (protein quant)
data <- read.csv("../step1_QC_impute_normalize/proteindf_raw_withNaN_imputed_normalized.csv", row.names = 1)  # Assuming Protein.Accession is in the first column


# Transpose the data and log transform and order them so it align later with metadata order
data <- t(data)
data <- log2(data)
data <- data[order(rownames(data)), ]


# Read CSV file for metadata
metadata <- read.csv("meta.csv", row.names = 'BatchName')  # Assuming TMT.channel is in the first column


### Drop missing values in metadata (remove the IRs) but before that store the info removed
exclude_in_metadata <- metadata[!complete.cases(metadata), ]
excluded_data <- data[rownames(exclude_in_metadata), ]

metadata <- metadata[complete.cases(metadata), ]
metadata <- metadata[order(rownames(metadata)), ]


# Only keep the common channel (removes IRs)
common_channel <- intersect(rownames(data), rownames(metadata))
data <- data[common_channel, ]
metadata <- metadata[common_channel, ]


hist(metadata$PMI, breaks = 20)

# bin the PMI
custom_breaks_PMI <- c(-Inf, 3, 5, Inf)
labels_PMI <- c(3, 4, 5)
metadata$PMI <- cut(metadata$PMI, breaks = custom_breaks_PMI, labels = labels_PMI, right = FALSE)
table(metadata$PMI)


# bin the Age
custom_breaks_age <- c(-Inf, 70, 90, Inf)
labels_age <- c(70, 80, 90)
metadata$Age <- cut(metadata$Age, breaks = custom_breaks_age, labels = labels_age, right = FALSE)
table(metadata$Age)

# Make data frame for corrected data to save
corrected_data <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
colnames(corrected_data) <- colnames(data)
rownames(corrected_data) <- rownames(data)


# Fit multivariate linear model for each protein and obtain residuals
for (protein in colnames(data)) {
  # Extract variables for the model
  y <- data[, protein]
  age <- metadata$Age
  sex <- metadata$Sex
  pmi <- metadata$PMI
  batch <- metadata$Batch
  
  # Create a data frame for the model
  model_data <- data.frame(response_variable = y, age = age, sex = sex, pmi = pmi, batch = batch)
  
  # Fit multivariate linear model
  model <- MASS::rlm(response_variable ~ age + sex + pmi + batch, data = model_data)

    # Obtain the residuals and store in corrected_data
  corrected_data[, protein] <- residuals(model) + coef(model)[1]
}


before_correction <-data
after_correction <- corrected_data


# Calculate correlation and correlation p value as well as create density plot
create_density_plot <- function(covar, variable_title) {
  # Calculate correlation and p-value
  ProExp_Quality_cor_bf <- cor(before_correction, covar, use = "p")
  ProExp_Quality_P_bf <- corPvalueStudent(ProExp_Quality_cor_bf, dim(before_correction)[1])
  
  ProExp_Quality_cor_af <- cor(after_correction, covar, use = "p")
  ProExp_Quality_P_af <- corPvalueStudent(ProExp_Quality_cor_af, dim(after_correction)[1])
  
  # KS test
  ks <- ks.test(ProExp_Quality_P_af, ProExp_Quality_P_bf)
  
  # Plot data
  ProExp_Quality_P_af <- as.data.frame(ProExp_Quality_P_af)
  ProExp_Quality_P_bf <- as.data.frame(ProExp_Quality_P_bf)
  
  dta <- cbind(ProExp_Quality_P_bf, ProExp_Quality_P_af)
  colnames(dta) <- c("before", "after")
  
  # Melting the data
  data_melted <- melt(dta, id.vars = NULL, variable.name = "variable", value.name = "p.value")
  
  # Creating the density plot
  ggplot(data_melted, aes(x = p.value, color = variable)) +
    geom_density(alpha = 0.5, size = 1.5) +  # Adjust size as needed
    theme_minimal() +
    theme(
      axis.line = element_line(size = 1.0),  # Adjust thickness of axes lines
      axis.ticks = element_line(size = 1.0),  # Adjust thickness of tick marks
      axis.text.y = element_text(size = 12),  # Adjust label size of y-axis
      axis.title.y = element_text(size = 14),  # Adjust size of y-axis title
      axis.text.x = element_text(size = 12),  # Adjust label size of y-axis
      axis.title.x = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 15)  # Adjust size of title and center it
    ) +
    ggtitle(paste("Before and After Correction:", variable_title))+
    labs(subtitle = paste("KS Test: D =", round(ks$statistic, 2), "p =", format.pval(ks$p.value, digits = 3)))
}

# Call the function for each variable
# (sex, age, pmi and batch are the covariate vectors defined in the loop above)
create_density_plot(sex, "Sex")
ggsave("Before_and_After_Correction_Sexv2.svg", height = 4.5, width = 4.5)

create_density_plot(age, "Age")
ggsave("Before_and_After_Correction_Agev2.svg", height = 4.5, width = 4.5)

create_density_plot(pmi, "PMI")
ggsave("Before_and_After_Correction_PMIv2.svg", height = 4.5, width = 4.5)

create_density_plot(batch, "Batch")
ggsave("Before_and_After_Correction_batchv2.svg", height = 4.5, width = 4.5)


# Merge corrected_data with excluded_data
complete_correcteddata <- rbind(corrected_data, excluded_data)
complete_correcteddata <- complete_correcteddata[order(rownames(complete_correcteddata)), ]
write.csv(complete_correcteddata,"./covariate_sex_age_pmi_combined_corrected.csv")

# Optional: merge metadata with corrected data (used to build the step 3 input table)
# metadata_complete <- read.csv("meta.csv", row.names = 'BatchName')
# metadata_complete <- metadata_complete[order(rownames(metadata_complete)), ]
# merged_data_meta_completedata <- cbind(metadata_complete, complete_correcteddata)
# 
# merged_data_meta_completedata$Paninfo <- paste("BS13", row.names(merged_data_meta_completedata), sep = "")
# write.csv(merged_data_meta_completedata,"./Metainfo_with_covariate_sex_age_pmi_combined_correcteddata.csv")
# 
# 
# row.names(merged_data_meta_completedata) <- merged_data_meta_completedata$Paninfo
# panData <- merged_data_meta_completedata[, -(1:10)]
# panData <- t(panData)
# write.csv(panData,"./BS13panData.csv")
