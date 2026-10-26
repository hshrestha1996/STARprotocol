library(Biobase)

###################Discovery cohort#####################################################
### FDR###########
sig <- read.csv("DE_BannerSun_VaD_P0.01_FCz2p5.csv", row.names = 1)
filtered_z <- 2.5
sig <- sig[abs(sig$VAD.Ctl_logFC.z)>=filtered_z,]

##################################################################

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

sig_dta <- dta[intersect(rownames(sig), rownames(dta)),]
VAD <- sample_meta[sample_meta$Diagnosis =="VAD",]
sig_dta_VAD <- sig_dta[,intersect(rownames(VAD), colnames(sig_dta))]
#write.csv(sig_dta_VAD, "DE_VADdata_P001_logFCz2p5.csv")

######################Work on subtyping####################
###----------load significant result--------------------######
library(WGCNA)
yan <- t(sig_dta_VAD)

##-------z-score---------------#
rnames<- rownames(yan)
yan <- apply(t(yan), 1, scale)
rownames(yan) <- rnames

####-----------------Check for outliers and remove them-------------#
sampleTree = hclust(dist(yan), method = "average") ##filter with average than complete

#pdf(file = "./sampleClustering.pdf", width = 200, height = 9); 
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 26, col = "red") #choose best value for you
#dev.off()
#------------keep the cluster below the red line---------------------
clust = cutreeStatic(sampleTree, cutHeight = 26, minSize = 10)
table(clust)
keepSamples = (clust== 1) #Cluster 1 contains the samples we want.
#keepSamples <- clust %in% c(1, 2)
yan = yan[keepSamples, ]
dim(yan)
##--------------------------------------------------------#
##--------------------------------------------------------#


###----------run ConsensusClusterPlus--------------------######
yan <- yan[complete.cases(yan), ]
yan <- t(yan)
boxplot(yan)


library(ConsensusClusterPlus)
library(MSnSet.utils)
library(ggplot2)
library(svglite)
set.seed(123)

Linkages <- c('ward.D2') #'average', 'complete', 
distances <- c('pearson') #, 'spearman', 'euclidean'
clusterAlgs <- c('pam') #,'km'
for (distance in distances) {
  for (Linkage in Linkages) {
    for (clusterAlg in clusterAlgs) {
      # Create a unique title for each combination
      title <- paste0("VaD_P0012p5z_", distance, "_", Linkage, "_", clusterAlg)
      
      # Run the ConsensusClusterPlus function with the current combination
      results <- ConsensusClusterPlus(
        yan,
        maxK = 10,
        reps = 1000,
        pItem = 0.9,
        pFeature = 1,
        title = title,
        distance = distance,
        writeTable = TRUE,
        seed = 123,
        plot = 'pdf',
        verbose = TRUE,
        clusterAlg = clusterAlg,
        innerLinkage = Linkage,
        finalLinkage = Linkage
      )
      
      icl = calcICL(results,title=title,plot="pdf", writeTable=TRUE)
      
      data <- yan
      data <- t(data)
      #data <- log2(data)
      data <- data[order(rownames(data)), ]
      
      
      for (i in 2:10) {
        # Create the meta data frame
        meta <- as.data.frame(results[[i]]$consensusClass)
        meta["subtype"] <- paste0("Subtype", meta$`results[[i]]$consensusClass`)
        metadata <- meta
        
        # Clean the metadata
        metadata <- metadata[complete.cases(metadata), ]
        metadata <- metadata[order(rownames(metadata)), ]
        
        # Keep only the common channel
        common_channel <- intersect(rownames(data), rownames(metadata))
        data_filtered <- data[common_channel, ]
        metadata <- metadata[common_channel, ]
        
        # Correct the data
        corrected_data <- t(data_filtered)
        corrected_data <- as.matrix(corrected_data)
        
        # Create MSnSet object
        msnset <- MSnSet(exprs = corrected_data)
        annotated_data <- AnnotatedDataFrame(metadata)
        phenoData(msnset) <- annotated_data
        
        # Plot PCA and save the plot
        plot_pca(msnset, phenotype = "subtype", legend_title = "Subtype")
        ggsave(paste0(title,"/km_", i, "subtypes.svg"),width = 4,height = 4,dpi = 300)
        #plot_pca(msnset, phenotype = "subtype", legend_title = "Subtype", components = c(1,3))
        #ggsave(paste0(title,"/km_", i, "subtypes_c13.png"),width = 4,height = 4,dpi = 300)
      }
      
    }}
  
}
