library(tidyverse)
library(Biobase)
library(limma)

##----------Load quant data ---------#########################
df <- readxl::read_excel("Supplement_table_VADv0.0.0.xlsx", sheet = "Sheet1",range = "B17:CP13672",col_names = FALSE)

##---------generate meta data ---------#########################
# Transpose the first two rows for sample metadata
sample_meta <- as.data.frame(t(df[1:2, ]))
# Assign the second row (now first column in transposed data) as column names and  Remove the now redundant first row
colnames(sample_meta) <- sample_meta[1, ]
sample_meta <- sample_meta[-1, ]
rownames(sample_meta) <- sample_meta[,2]
sample_meta <- sample_meta[complete.cases(sample_meta), ]

##---------generate protein info data and abundance matrix---------#########################
colnames(df) <- df[2, ]
df <- df[-c(1,2), ]
df <- as.data.frame(df)
row.names(df) <- df[,1]
protein_info <- df[,c(1:9)]
df <- df[,-c(1:9)]

# Convert all columns in df to numeric
df[] <- lapply(df, as.numeric)
# Apply log2 transformation to the numeric data
dta <- log2(df)
dta <- as.matrix(dta)

# ############################## limma DE analysis ###############################
# # create ExpressionSet object used for limma analysis
df_quan_eSet <- ExpressionSet(assayData = dta)

# create design matrix (allocate samples into groups)
design <- model.matrix(~ 0 + factor(sample_meta$Diagnosis, levels = unique(sample_meta$Diagnosis)))
#colnames(design) <- str_c('group_', unique(sample_meta$Group))
colnames(design) <- unique(sample_meta$Diagnosis)

# fit linear model for each protein
fit <- lmFit(df_quan_eSet, design)

unique(sample_meta$Diagnosis)

cont.matrix <- makeContrasts(contrasts = c("AD-Ctl", "VAD-Ctl"), 
                             levels=design)

# compute contrasts from linear model fit
fit2 <- contrasts.fit(fit, cont.matrix)

# empirical Bayes statistics for differential expression
fit2 <- eBayes(fit2)

# extract differential expression results
df_DE_results <- data.frame(Identifier = row.names(df))
for (i in 1:ncol(cont.matrix)) {
  sub_DE_results <- topTable(fit2, coef = i, sort.by="P", number = Inf) %>% 
    select(logFC, P.Value, adj.P.Val) %>% rownames_to_column(.,"Identifier")
  colnames(sub_DE_results)[2:4] <- str_c(colnames(cont.matrix)[i], colnames(sub_DE_results)[2:4], sep = "_")
  df_DE_results <- merge(df_DE_results, sub_DE_results, by='Identifier')
}

# add protein information to the DE results table
rownames(df_DE_results)<- df_DE_results$Identifier
df_DE_results<- df_DE_results[,-1]

df_DE_results <- df_DE_results[match(rownames(protein_info), rownames(df_DE_results)), ]
all(rownames(df_DE_results) == rownames(protein_info))
DE_results <- cbind(protein_info, df_DE_results)
# output the DE results
writexl::write_xlsx(DE_results, path = paste0("DE_result_bannersunVaD.xlsx"))


##################################################################
######-------Check 2SD to use---------###########
# Define the fold change column
fc_column <- "VAD-Ctl_logFC"
thresholds <- seq(0, 1, by = 0.01)
DE <- df_DE_results[order(abs(df_DE_results[[fc_column]]), decreasing = TRUE), ]
SD <- sapply(thresholds, function(thresh) {
  sd(DE[ceiling(thresh * nrow(DE)):length(DE[[fc_column]]),][[fc_column]])
})

results_df <- data.frame(Threshold = thresholds, SD = SD)
results_df$removed <- results_df$Threshold*100
results_df[is.na(results_df)] <- 0
results_df <- round(results_df, digits = 4)


ggplot(results_df, aes(x = removed, y = SD)) +
  geom_line() +
  geom_point() +
  #geom_hline(yintercept = 0.99, color = "red", linetype = "dashed", size = 0.5) + 
  labs(title = "SD of VADvs.Ctl log2FC",
       x = "% of the data removed",
       y = "SD calculation") +
  theme_minimal()+  # Minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around the panel
    axis.ticks = element_line(color = "black"),  # Color of ticks
    axis.title = element_text(size = 12, face = "bold"),  # Bold axis titles
    axis.text = element_text(size = 10),  # Size of axis text
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # Center title with size and bold
  )

select_removal_perc = 10
sd <- results_df[results_df$removed==select_removal_perc,2]

# Set up the plotting area for two rows
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))  # Adjust margins if necessary

hist(DE[[fc_column]], breaks = 100, main = "Histogram of VAD vs Ctl log2FC (in complete)",xlab = "VAD vs Ctl log2FC",
     col = "grey",border = "black")

filtered_dta <- DE[ceiling((select_removal_perc/100) * nrow(DE)):length(DE[[fc_column]]),] #DE[abs(DE[[fc_column]]) < 2*sd,]
hist(filtered_dta[[fc_column]], breaks = 100, main = "Histogram of VAD vs Ctl log2FC (10% data filtered)",xlab = "VAD vs Ctl log2FC",
     col = "lightblue",border = "black")

dev.off()

##################################################################
#####adding log2FC-z
df_DE_results$`VAD-Ctl_logFC-z` <- (df_DE_results$`VAD-Ctl_logFC`)/sd
##################################################################
######-------Check 2SD to use---------###########
# Define the fold change column
fc_column <- "AD-Ctl_logFC"
thresholds <- seq(0, 1, by = 0.01)
DE <- df_DE_results[order(abs(df_DE_results[[fc_column]]), decreasing = TRUE), ]
SD <- sapply(thresholds, function(thresh) {
  sd(DE[ceiling(thresh * nrow(DE)):length(DE[[fc_column]]),][[fc_column]])
})

results_df <- data.frame(Threshold = thresholds, SD = SD)
results_df$removed <- results_df$Threshold*100
results_df[is.na(results_df)] <- 0
results_df <- round(results_df, digits = 4)


ggplot(results_df, aes(x = removed, y = SD)) +
  geom_line() +
  geom_point() +
  #geom_hline(yintercept = 0.99, color = "red", linetype = "dashed", size = 0.5) + 
  labs(title = "SD of ADvs.Ctl log2FC",
       x = "% of the data removed",
       y = "SD calculation") +
  theme_minimal()+  # Minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around the panel
    axis.ticks = element_line(color = "black"),  # Color of ticks
    axis.title = element_text(size = 12, face = "bold"),  # Bold axis titles
    axis.text = element_text(size = 10),  # Size of axis text
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # Center title with size and bold
  )

select_removal_perc = 10
sd <- results_df[results_df$removed==select_removal_perc,2]

# Set up the plotting area for two rows
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))  # Adjust margins if necessary

hist(DE[[fc_column]], breaks = 100, main = "Histogram of AD vs Ctl log2FC (in complete)",xlab = "AD vs Ctl log2FC",
     col = "grey",border = "black")

filtered_dta <- DE[ceiling((select_removal_perc/100) * nrow(DE)):length(DE[[fc_column]]),] #DE[abs(DE[[fc_column]]) < 2*sd,]
hist(filtered_dta[[fc_column]], breaks = 100, main = "Histogram of AD vs Ctl log2FC (10% data filtered)",xlab = "AD vs Ctl log2FC",
     col = "lightblue",border = "black")

dev.off()

##################################################################
#####adding log2FC-z
df_DE_results$`AD-Ctl_logFC-z` <- (df_DE_results$`AD-Ctl_logFC`)/sd
write.csv(df_DE_results, "DE_result_bannersun_with_logFCzscore_AD_VAD.csv")
