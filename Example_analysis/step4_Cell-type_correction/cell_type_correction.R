# Step 4.0 - Cell-type correction
# For each protein, fit a robust linear model (MASS::rlm) on the five
# cell-type marker scores (astrocyte, endothelial, oligodendrocyte, neuron,
# microglia), cap each cell-type coefficient at |beta| <= 1, and subtract the
# cell-type effect while keeping the intercept.
#
# Inputs:
#   ../step2_covariate_correction/covariate_sex_age_pmi_combined_corrected.csv
#   meta_cell_type_correction.csv   (made by cell_type_marker_before.ipynb)
# Output:
#   Celltype_corrected_rlm_beta_less1_BannerSunVaDcohort.csv
#
# Run from inside the step4.0_Cell-type_correction folder.

library(MASS)

# Protein quant (rows = TMT channels, columns = proteins)
data <- read.csv("../step2_covariate_correction/covariate_sex_age_pmi_combined_corrected.csv",
                 row.names = 1)
data <- data[order(rownames(data)), ]

# Metadata with cell-type marker scores (rows = TMT channels)
metadata <- read.csv("meta_cell_type_correction.csv", row.names = 1)

# Drop internal-reference (IR) channels, which have missing metadata
metadata <- metadata[complete.cases(metadata), ]
metadata <- metadata[order(rownames(metadata)), ]

# Keep only the cell-type score columns and mean-centre them
metadata <- metadata[, -(1:10)]
metadata <- sweep(metadata, 2, colMeans(metadata), "-")
metadata <- as.data.frame(round(metadata, 2))

# Optional QC: distribution of each centred cell-type score
hist(metadata$Average_ast_markerabundance)
hist(metadata$Average_end_markerabundance)
hist(metadata$Average_mic_markerabundance)
hist(metadata$Average_neu_markerabundance)
hist(metadata$Average_oil_markerabundance)

# Keep channels present in both tables (removes IRs from the protein data)
common_channel <- intersect(rownames(data), rownames(metadata))
data <- data[common_channel, ]
metadata <- metadata[common_channel, ]

corrected_data <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data),
                                       dimnames = list(rownames(data), colnames(data))))

threshold <- 1  # maximum absolute cell-type coefficient

for (protein in colnames(data)) {
  y  <- data[, protein]
  a1 <- metadata$Average_ast_markerabundance
  e2 <- metadata$Average_end_markerabundance
  o3 <- metadata$Average_oil_markerabundance
  n4 <- metadata$Average_neu_markerabundance
  m5 <- metadata$Average_mic_markerabundance

  a <- summary(rlm(y ~ a1 + e2 + o3 + n4 + m5, maxit = 30))

  # Cap each cell-type coefficient (rows 2-6) at +/- threshold
  for (j in 2:nrow(a$coefficients)) {
    if (abs(a$coefficients[j, 1]) > threshold) {
      a$coefficients[j, 1] <- sign(a$coefficients[j, 1]) * threshold
    }
  }

  pass <- which(!is.na(y))
  mod_residuals <- y[pass] - (a$coefficients[1, 1] +
                                a1[pass] * a$coefficients[2, 1] +
                                e2[pass] * a$coefficients[3, 1] +
                                o3[pass] * a$coefficients[4, 1] +
                                n4[pass] * a$coefficients[5, 1] +
                                m5[pass] * a$coefficients[6, 1])
  corrected_data[, protein] <- mod_residuals + a$coefficients[1, 1]
}

write.csv(corrected_data, "Celltype_corrected_rlm_beta_less1_BannerSunVaDcohort.csv")
