setwd("C:/Users/Bacha Zada/Desktop/SysBio_internship/Ranalysis/")

data <- read.delim("ACC_peakCalls.txt", header = TRUE)

###############

# Check how many different genomic annotations exist
unique_annotations <- unique(data$annotation)
head(unique_annotations)

####################


# How many fragments for each genomic annotations
# Count the number of fragments for each genomic annotation
annotation_counts <- table(data$annotation)
print(annotation_counts)



# Visualize the results
# Bar plot of the counts
barplot(annotation_counts, main="Fragment Counts per Genomic Annotation",
        xlab="Genomic Annotation", ylab="Count", col="skyblue", las=3)
#####################

# Check the length distribution of ATAC-seq fragments
data$fragment_length <- data$end - data$start

# Summary of fragment lengths
summary(data$fragment_length)

# Visualize the length distribution
hist(data$fragment_length, breaks=50, main="Length Distribution of ATAC-seq Fragments",
     xlab="Fragment Length", col="lightgreen")
#################

# Extract sequences for promoters and/or distal regions

#if (!requireNamespace("BiocManager", force = TRUE))
# install.packages("BiocManager", lib = Sys.getenv("R_LIBS_USER"))

#BiocManager::install("GenomicRanges", lib = Sys.getenv("R_LIBS_USER"), force = TRUE)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", lib = Sys.getenv("R_LIBS_USER"), force = TRUE)


library(GenomicRanges)
library(XVector)
library(Biostrings)
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

# Filter promoters and distal regions and other regions
promoters <- subset(data, annotation == "Promoter")
distal_regions <- subset(data, annotation == "Distal")
Exons <- subset(data,annotation == "Exon")
UTR3 <- subset(data,annotation == "3' UTR" )
UTR5 <- subset(data, annotation == "5' UTR")



# Create GRanges objects for each annotation group
gr_promoters <- GRanges(seqnames = promoters$seqnames, ranges = IRanges(start = promoters$start, end = promoters$end))
gr_distal_regions <- GRanges(seqnames = distal_regions$seqnames, ranges = IRanges(start = distal_regions$start, end = distal_regions$end))
gr_Exons <- GRanges(seqnames = Exons$seqnames, ranges = IRanges(start = Exons$start, end = Exons$end))
gr_UTR3 <- GRanges(seqnames = UTR3$seqnames, ranges = IRanges(start = UTR3$start, end = UTR3$end))
gr_UTR5 <- GRanges(seqnames = UTR5$seqnames, ranges = IRanges(start = UTR5$start, end = UTR5$end))

# Load the human genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# Extract sequences for the regions of interest
promoter_sequences <- getSeq(genome, gr_promoters)
distal_sequences <- getSeq(genome, gr_distal_regions)
Exon_sequences <- getSeq(genome, gr_Exons)
UTR3_sequences <- getSeq(genome, gr_UTR3)
UTR5_sequences <- getSeq(genome, gr_UTR5)
print(promoters)

###########
# Function to calculate nucleotide frequency and GC/AT content
calculate_content <- function(sequences) {
  seq_string <- as.character(sequences)
  seq_concat <- paste(seq_string, collapse = "")
  seq_length <- nchar(seq_concat)
  
  A_count <- sum(charToRaw(seq_concat) == charToRaw("A"))
  T_count <- sum(charToRaw(seq_concat) == charToRaw("T"))
  G_count <- sum(charToRaw(seq_concat) == charToRaw("G"))
  C_count <- sum(charToRaw(seq_concat) == charToRaw("C"))
  
  GC_content <- (G_count + C_count) / seq_length * 100
  AT_content <- (A_count + T_count) / seq_length * 100
  
  list(A = A_count, T = T_count, G = G_count, C = C_count, GC_content = GC_content, AT_content = AT_content)
}

########
# Calculate content for each annotation group
promoter_content <- calculate_content(promoter_sequences)
distal_content <- calculate_content(distal_sequences)
Exon_content <- calculate_content(Exon_sequences)
UTR3_content <- calculate_content(UTR3_sequences)
UTR5_content <- calculate_content(UTR5_sequences)
############

# Combine results into a data frame
results <- data.frame(
  Annotation = rep(c("Promoter", "Distal", "Exon", "UTR3", "UTR5"), each = 6),
  Nucleotide = rep(c("A", "T", "G", "C", "GC_content", "AT_content"), times = 5),
  Frequency = c(
    promoter_content$A, promoter_content$T, promoter_content$G, promoter_content$C, promoter_content$GC_content, promoter_content$AT_content,
    distal_content$A, distal_content$T, distal_content$G, distal_content$C, distal_content$GC_content, distal_content$AT_content,
    Exon_content$A, Exon_content$T, Exon_content$G, Exon_content$C, Exon_content$GC_content, Exon_content$AT_content,
    UTR3_content$A, UTR3_content$T, UTR3_content$G, UTR3_content$C, UTR3_content$GC_content, UTR3_content$AT_content,
    UTR5_content$A, UTR5_content$T, UTR5_content$G, UTR5_content$C, UTR5_content$GC_content, UTR5_content$AT_content
  )
)
library(ggplot2)
# Plot nucleotide frequencies
ggplot(subset(results, Nucleotide %in% c("A", "T", "G", "C")), aes(x = Annotation, y = Frequency, fill = Nucleotide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Nucleotide Frequencies by Annotation Group", x = "Annotation Group", y = "Frequency")

# Plot GC and AT content percentages
ggplot(subset(results, Nucleotide %in% c("GC_content", "AT_content")), aes(x = Annotation, y = Frequency, fill = Nucleotide)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "GC and AT Content by Annotation Group", x = "Annotation Group", y = "Content Percentage")
####################


# week 2

# Load the new dataset (GBM_Masterfile.txt)
gbm_data <- read.delim("GBM_peakCalls.txt", header = TRUE)
head(gbm_data)

#Convert both your current dataset and the GBM dataset into GRanges objects for comparison.
gr_current <- GRanges(seqnames = data$seqnames, ranges = IRanges(start = data$start, end = data$end))

# Convert GBM dataset to GRanges
gr_GBM <- GRanges(seqnames = gbm_data$seqnames, ranges = IRanges(start = gbm_data$start, end = gbm_data$end))
################
# Calculate unique and overlapping regions using GRanges methods
overlapping_regions <- intersect(gr_current, gr_GBM)
unique_current <- setdiff(gr_current, gr_GBM)
unique_GBM <- setdiff(gr_GBM, gr_current)

# Count the number of overlapping and unique regions
num_overlapping_fragments <- length(overlapping_regions)
num_unique_current <- length(unique_current)
num_unique_GBM <- length(unique_GBM)

# Print out the correct numbers
cat("Number of unique fragments in your current dataset: ", num_unique_current, "\n")
cat("Number of unique fragments in the GBM dataset: ", num_unique_GBM, "\n")
cat("Number of overlapping fragments: ", num_overlapping_fragments, "\n")

library(VennDiagram)

# Define the areas
area1 <- length(unique_current)  # Unique Current data
area2 <- length(unique_GBM)  # Unique GBM data
cross.area <- length(overlapping_regions)  # Overlapping regions

# Create a Venn diagram using `draw.pairwise.venn`
venn.plot <- draw.pairwise.venn(
  area1 = area1,
  area2 = area2,
  cross.area = cross.area,
  category = c("Current_data", "GBM_data"),
  fill = c("blue", "red"),  # Colors for each set
  lty = "blank",  # No outline around the circles
  cex = 2,  # Size of the counts
  cat.cex = 2,  # Size of the category labels
  cat.pos = c(180, 0),  # Positions of the category labels
  cat.dist = c(0.03, 0.04),  # Distance of the labels from the circles
  cat.col = c("blue", "red"),  # Colors for the category labels
  scaled = TRUE  # Scales the circles based on the area size
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Add a label for the overlapping region
grid.text("Overlapping_regions", x = 0.53, y = 0.62, gp = gpar(fontsize = 17, col = "black"))
###########################################

#week 3

# Load the count vector data
acc_count_vector <- readRDS("ACC_countVector.RDS")
gbm_count_vector <- readRDS("GBM_countVector.RDS")

# Load the genomic annotation data (e.g., from a BED file or another format)
acc_data <- read.delim("ACC_peakCalls.txt", header = TRUE)
gbm_data <- read.delim("GBM_peakCalls.txt", header = TRUE)
head(acc_data)
# Filter ACC data by genomic annotations
acc_promoters <- acc_count_vector[acc_data$annotation == "Promoter", ]
acc_distal <- acc_count_vector[acc_data$annotation == "Distal", ]
acc_exons <- acc_count_vector[acc_data$annotation == "Exon", ]
acc_utr3 <- acc_count_vector[acc_data$annotation == "3' UTR", ]
acc_utr5 <- acc_count_vector[acc_data$annotation == "5' UTR", ]

# Similarly for GBM data
gbm_promoters <- gbm_count_vector[gbm_data$annotation == "Promoter", ]
gbm_distal <- gbm_count_vector[gbm_data$annotation == "Distal", ]
gbm_exons <- gbm_count_vector[gbm_data$annotation == "Exon", ]
gbm_utr3 <- gbm_count_vector[gbm_data$annotation == "3' UTR", ]
gbm_utr5 <- gbm_count_vector[gbm_data$annotation == "5' UTR", ]

# Sum transcription factor counts for each genomic annotation in ACC data
acc_promoter_sums <- colSums(acc_promoters)
acc_distal_sums <- colSums(acc_distal)
acc_exon_sums <- colSums(acc_exons)
acc_utr3_sums <- colSums(acc_utr3)
acc_utr5_sums <- colSums(acc_utr5)

# Combine the results into a single data frame
acc_tf_frequencies <- data.frame(
  TF_Class = names(acc_promoter_sums),
  Promoter = acc_promoter_sums,
  Distal = acc_distal_sums,
  Exon = acc_exon_sums,
  UTR3 = acc_utr3_sums,
  UTR5 = acc_utr5_sums
)

# Similarly for GBM data
gbm_promoter_sums <- colSums(gbm_promoters)
gbm_distal_sums <- colSums(gbm_distal)
gbm_exon_sums <- colSums(gbm_exons)
gbm_utr3_sums <- colSums(gbm_utr3)
gbm_utr5_sums <- colSums(gbm_utr5)

gbm_tf_frequencies <- data.frame(
  TF_Class = names(gbm_promoter_sums),
  Promoter = gbm_promoter_sums,
  Distal = gbm_distal_sums,
  Exon = gbm_exon_sums,
  UTR3 = gbm_utr3_sums,
  UTR5 = gbm_utr5_sums
)
library(reshape2)

# Melt the ACC data for easier plotting
acc_tf_frequencies_melted <- melt(acc_tf_frequencies, id.vars = "TF_Class", variable.name = "Annotation", value.name = "Frequency")

# Plot for ACC data
ggplot(acc_tf_frequencies_melted, aes(x = TF_Class, y = Frequency, fill = Annotation)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Transcription Factor Frequencies for ACC Data", x = "Transcription Factor Class", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Melt the GBM data for easier plotting
gbm_tf_frequencies_melted <- melt(gbm_tf_frequencies, id.vars = "TF_Class", variable.name = "Annotation", value.name = "Frequency")

# Plot for GBM data
ggplot(gbm_tf_frequencies_melted, aes(x = TF_Class, y = Frequency, fill = Annotation)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Transcription Factor Frequencies for GBM Data", x = "Transcription Factor Class", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine ACC and GBM data
combined_tf_frequencies <- rbind(
  cbind(acc_tf_frequencies_melted, Dataset = "ACC"),
  cbind(gbm_tf_frequencies_melted, Dataset = "GBM")
)

# Plot comparing ACC and GBM data
ggplot(combined_tf_frequencies, aes(x = TF_Class, y = Frequency, fill = Annotation)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Dataset) +
  theme_minimal() +
  labs(title = "Transcription Factor Frequencies Across Genomic Annotations", x = "Transcription Factor Class", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######
library(corrplot)
library(pheatmap)

# Compute correlation matrix
cor_matrix <- cor(acc_count_vector, use = "pairwise.complete.obs")

# Plot correlation heatmap using pheatmap
pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Correlation Heatmap of Chromatin Accessibility Profiles")


# Assuming that acc_promoters and gbm_promoters have been created and contain promoter-specific data
# Combine ACC and GBM promoter data into one matrix
combined_promoters <- rbind(acc_promoters, gbm_promoters)

# Calculate correlation matrix
promoter_cor_matrix <- cor(combined_promoters, use = "pairwise.complete.obs", method = "pearson")

# Load heatmap library
library(pheatmap)

# Plot the correlation heatmap
pheatmap(promoter_cor_matrix, 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         main = "Correlation Heatmap of Promoter Regions (ACC and GBM)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE)  # Adjust the color palette as per your preference

###################
# Combine ACC and GBM data
combined_data <- rbind(acc_tf_frequencies[, -1], gbm_tf_frequencies[, -1])

# Compute the correlation matrix
tf_cor_matrix <- cor(combined_data, use = "pairwise.complete.obs", method = "pearson")

# View the correlation matrix
head(tf_cor_matrix)
# Perform hierarchical clustering
tf_hclust <- hclust(as.dist(1 - tf_cor_matrix), method = "complete")

# Plot dendrogram
plot(tf_hclust, labels = rownames(tf_cor_matrix), main = "Hierarchical Clustering of Transcription Factors")

library(pheatmap)

# Visualize the correlation matrix with a heatmap
pheatmap(tf_cor_matrix, 
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Correlation Heatmap of Transcription Factors",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE)

# Perform PCA on the combined transcription factor data
pca_result <- prcomp(combined_data, scale. = TRUE)

# Extract the first two principal components
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  TF_Class = rownames(pca_result$x)
)

# Plot the PCA results
ggplot(pca_data, aes(x = PC1, y = PC2, label = TF_Class)) +
  geom_point(aes(color = TF_Class)) +
  # geom_text(aes(label = TF_Class), vjust = 3.5, size = 3) +
  theme_minimal() +
  labs(title = "PCA of Transcription Factors",
       x = "Principal Component 1",
       y = "Principal Component 2")

library(Rtsne)

# Perform t-SNE
set.seed(42)
tsne_result <- Rtsne(as.matrix(combined_data), perplexity = 2, check_duplicates = FALSE)

# Create a data frame for the t-SNE results
tsne_data <- data.frame(
  TSNE1 = tsne_result$Y[, 1],
  TSNE2 = tsne_result$Y[, 2],
  TF_Class = rownames(combined_data)
)

# Plot the t-SNE results
ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, label = TF_Class)) +
  geom_point(aes(color = TF_Class)) +
  # geom_text(aes(label = TF_Class), vjust = 1.5, size = 3) +
  theme_minimal() +
  labs(title = "t-SNE of Transcription Factors",
       x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2")

