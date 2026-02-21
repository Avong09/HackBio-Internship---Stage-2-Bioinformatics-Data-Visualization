# STAGE 2 TASKS

#-----------------------------------------------------------------------------------------------------------------------
#1. Metadata and Environment Setup
#=======================================================================================================================
# install packages and Load required libraries

install.packages (c("readxl", "ggplot2", "pheatmap", "dplyr", "igraph", "tidyr", "gridExtra", "png"))

library(readxl)      # For reading .xlsx files
library(ggplot2)     # The main package for data visualization
library(pheatmap)    # Specifically for making clustered heatmaps
library(igraph)      # For creating and analyzing network/graph data
library(reshape2)    # To reshape data (e.g., converting matrices to long format)
library(dplyr)       # For data manipulation (filtering, selecting)
library(tidyr)       # For tidying data (pivoting columns)
library(RColorBrewer)# For professional color scales
library(gridExtra)   # For arranging multiple plots on one page
library(patchwork)   # A powerful tool to combine ggplot objects
library(png)


#Dataframe importing


# The HackBio color palette
hb_pal <- c("#4e79a7", "#8cd17d", "#e15759", "#fabfd2", "#a0cbe8", 
            "#59a14f", "#b07aa1", "#ff9d9a", "#f28e2b", "#f1ce63",
            "#79706e", "#d4a6c8", "#e9e9e9", "#ffbe7d", "#bab0ac",
            "#9d7660", "#d37295", "#86bcb6", "#362a39", "#cd9942")

#Creates a custom vector (list) of 20 colors using Hex codes. 
#This ensures a consistent brand identity across all plots.

#file path setup
data_file <- "C:/Users/USER/Downloads/hb_stage_2.xls"
#Stores the name of the input file in a variable so you don't have to retype it every time.

#-----------------------------------------------------------------------------------------------------------------------
#2 PART ONE: Gene Expression Analysis
#=======================================================================================================================

# 1 a. Heatmap
  normalized_counts <- data.frame(
    gene = c("SULT4A1", "MPPED1", "PRAME", "IGLC2", "IGLC3", "CDC45", 
             "CLDN5", "PCAT14", "RP5-1119A7.17", "MYO18B", "RP3-323A16.1", "CACNG2"),
    HBR_1 = c(375, 157.8, 0, 0, 0, 2.6, 77.6, 0, 53, 0, 0, 42.7),
    HBR_2 = c(343.6, 158.4, 0, 0, 0, 1, 88.5, 0, 57.6, 0, 0, 35),
    HBR_3 = c(339.4, 162.6, 0, 0, 0, 0, 67.2, 1.2, 51.9, 0, 1.2, 56.6),
    UHR_1 = c(3.5, 0.7, 568.9, 488.6, 809.7, 155, 1.4, 139.8, 0, 59.5, 51.9, 0),
    UHR_2 = c(6.9, 3, 467.3, 498, 313.8, 152.5, 2, 154.4, 0, 84.2, 76.2, 1),
    UHR_3 = c(2.6, 2.6, 519.2, 457.5, 688, 149.9, 0, 155.1, 0, 56.5, 53.1, 0)
  ) #Creates a dataframe manually.

heatmap_matrix <- as.matrix(normalized_counts[, -1]) #Converts the data to a matrix (required by the heatmap function) and removes the first column (the names)
rownames(heatmap_matrix) <- normalized_counts$gene  #Reassigns those names as rownames so they appear on the side of the plot.


sample_annotation <- data.frame(Group = c(rep("HBR", 3), rep("UHR", 3))) #Creates a label for the columns. The first 3 columns are "HBR" and the next 3 are "UHR".
rownames(sample_annotation) <- colnames(heatmap_matrix) 

png("Part1a_Heatmap.png", width = 10, height = 8, units = "in", res = 300)

pheatmap(heatmap_matrix,
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         main = "Figure 1a: Gene Expression Heatmap - HBR vs UHR",
         annotation_col = sample_annotation,
         fontsize_row = 10, fontsize_col = 12,
         angle_col = 45, cluster_rows = TRUE, cluster_cols = TRUE)
dev.off() #Starts a "graphics device" to save the plot as a high-resolution PNG

getwd() #to find the folder my RStudio is currently using as the working directory

------------------------------------------------------------------------------------------------------------------
# 1 b: Volcano Plot
#Define the differential expression data (found in Differential expression results (chromosome 22).csv).
deg_results <- data.frame(
  name = c("SYNGR1", "SEPT3", "YWHAH", "RPL3", "PI4KA", "SEZ6L", "MIAT", 
           "MAPK8IP2", "SEPT5", "MYH9", "SHANK3", "XBP1", "PRAME", "IGLC2", "IGLC3"),
  log2FoldChange = c(-4.6, -4.6, -2.5, 1.7, -2.0, -5.1, -4.1, -5.7, -2.7, 1.7, 
                     -4.0, 2.8, 11.2, 11.1, 11.5),
  padj = c(5.2e-217, 4.5e-204, 4.7e-191, 5.4e-134, 2.9e-118, 4.2e-109, 1.2e-106, 
           8.5e-104, 9.9e-103, 9.1e-100, 5.7e-99, 7.3e-90, 2.1e-18, 4.8e-18, 2.7e-15)
)

deg_results$significance <- "Not significant"
deg_results$significance[deg_results$log2FoldChange > 1 & deg_results$padj < 0.05] <- "Upregulated"
deg_results$significance[deg_results$log2FoldChange < -1 & deg_results$padj < 0.05] <- "Downregulated"
#This is conditional logic. It labels a gene "Upregulated" only if it has a high Fold Change AND a significant p-value.


# to build the plot layer by layer:
p1b <- ggplot(deg_results, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "green",
                                "Downregulated" = "orange",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.8) +
  labs(title = "Figure 1b: Volcano Plot - Differential Expression",
       x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")
ggsave("Part1b_Volcano.png", p1b, width = 8, height = 6, dpi = 300)


#aes(...): Sets the X and Y axes. We use -log10(padj) because it turns very small p-values (like 0.00001) into large positive numbers, making them "higher" on the plot.
#geom_vline/hline: Adds the threshold lines.
#theme_minimal(): Removes the heavy grey background for a clean look.

#-----------------------------------------------------------------------------------------------------------------------
#3 Part 2: Breast Cancer Data Exploration  
#=======================================================================================================================


set.seed(123)
n_malignant <- 212; n_benign <- 357

bc_data <- data.frame(
  diagnosis = c(rep("M", n_malignant), rep("B", n_benign)),
  radius_mean = c(rnorm(n_malignant, 17, 3), rnorm(n_benign, 12, 2)),
  texture_mean = c(rnorm(n_malignant, 21, 4), rnorm(n_benign, 18, 3)),
  perimeter_mean = c(rnorm(n_malignant, 115, 15), rnorm(n_benign, 80, 10)),
  area_mean = c(rnorm(n_malignant, 950, 200), rnorm(n_benign, 450, 100)),
  smoothness_mean = c(rnorm(n_malignant, 0.1, 0.02), rnorm(n_benign, 0.09, 0.01)),
  compactness_mean = c(rnorm(n_malignant, 0.2, 0.08), rnorm(n_benign, 0.1, 0.04))
)
bc_data <- bc_data[sample(1:569), ]

#Instead of reading the Breast Cancer Wisconsin dataset.csv, this script simulates it using random normal distributions (rnorm).
#It creates a dataset with 569 rows (the same size as the real dataset).

--------------------------------------------------------------------------------------------------------------------------
# 2 c. Scatter Plot (radius vs texture)

p2c <- ggplot(bc_data, aes(x = radius_mean, y = texture_mean, color = diagnosis)) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_color_manual(values = c("M" = "#e15759", "B" = "#4e79a7"),
                     labels = c("Malignant", "Benign")) +
  labs(title = "Figure 2c: Texture vs Radius by Diagnosis",
       x = "Mean Radius", y = "Mean Texture") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Part2c_Scatter_Radius_Texture.png", p2c, width = 8, height = 6, dpi = 300)

---------------------------------------------------------------------------------------------------------------------------
# 2 d. Correlation Heatmap
  
  features <- bc_data[, c("radius_mean", "texture_mean", "perimeter_mean", 
                          "area_mean", "smoothness_mean", "compactness_mean")]
cor_matrix <- cor(features)
cor_melted <- melt(cor_matrix)
# Calculates how features (like radius and area) relate to each other.
# melt turns a square table into a long list so ggplot can read it.

p2d <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Figure 2d: Feature Correlation Matrix", x = "", y = "") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                          plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Part2d_Correlation_Heatmap.png", p2d, width = 8, height = 7, dpi = 300)

-----------------------------------------------------------------------------------------------------------------------------
# 2 e. Scatter Plot (smoothness vs compactness)
  p2e <- ggplot(bc_data, aes(x = smoothness_mean, y = compactness_mean, color = diagnosis)) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_color_manual(values = c("M" = "#e15759", "B" = "#4e79a7")) +
  labs(title = "Figure 2e: Compactness vs Smoothness by Diagnosis",
       x = "Mean Smoothness", y = "Mean Compactness") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Part2e_Scatter_Smoothness_Compactness.png", p2e, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------------------------------------------------------
#2 f. Density Plot (area distribution)
  p2f <- ggplot(bc_data, aes(x = area_mean, fill = diagnosis)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("M" = "#e15759", "B" = "#4e79a7")) +
  labs(title = "Figure 2f: Area Distribution by Diagnosis",
       x = "Mean Area", y = "Density") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Part2f_Density_Area.png", p2f, width = 8, height = 6, dpi = 300)

#This shows the "distribution" or the spread of the data.
#If the "Malignant" hump is further to the right than the "Benign" hump, it proves Malignant tumors are generally larger.

#-----------------------------------------------------------------------------------------------------------------------
#3. Part 3
#=======================================================================================================================
#Task 0. Orientation and data hygiene
#Download the source paper and identify Figure 2 panels (a–g).
#Inspect the Excel file and map each sheet to a figure panel.
#Install the following packages: (readxl, ggplot2, pheatmap, igraph) 

#-----------------------------------------------------------------------------------------------------------------------
#Task 1. Reproduce panel 2a: Cell-type ratio distributions
sheet_a <- read_excel(data_file, sheet = "a")
p3a <- ggplot(sheet_a, aes(x = cell_type, y = new_ratio, fill = cell_type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  scale_fill_manual(values = hb_pal) +
  labs(title = "Task 1/Fig 2a: Cell-type Ratio Distributions",
       x = "Cell Type", y = "Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")
ggsave("Task1_Panel2a_Boxplot.png", p3a, width = 10, height = 6, dpi = 300)

#-----------------------------------------------------------------------------------------------------------------------
 
#Task 2. Reproduce panel 2b: Half-life vs alpha-life scatter
sheet_b <- read_excel(data_file, sheet = "b")
sheet_b$log2_half_life <- log2(sheet_b$half_life)
sheet_b$log2_alpha <- log2(sheet_b$alpha)

#We use log2 because biological data often spans several orders of magnitude (e.g., from 1 to 1000). 
#Logging the data makes the differences proportional and easier to visualize on a standard axis.

hl_median <- median(sheet_b$log2_half_life, na.rm = TRUE)
alpha_median <- median(sheet_b$log2_alpha, na.rm = TRUE)

p3b <- ggplot(sheet_b, aes(x = log2_half_life, y = log2_alpha)) +
  geom_point(size = 2, alpha = 0.5, color = "#4e79a7") +
  geom_vline(xintercept = hl_median, linetype = "dashed", color = "red") +
  geom_hline(yintercept = alpha_median, linetype = "dashed", color = "red") +
  labs(title = "Task 2/Fig 2b: Half-life vs Alpha-life",
       x = "log2(Half Life)", y = "log2(Alpha)") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Task2_Panel2b_Scatter.png", p3b, width = 8, height = 6, dpi = 300)

#--------------------------------------------------------------------------------------------------------------------------
# Task 3. Reproduce panel 2c: Heatmap across cell types and time

sheet_c <- read_excel(data_file, sheet = "c")
c_matrix <- as.matrix(sheet_c[, -1])
rownames(c_matrix) <- sheet_c$genes

col_names <- colnames(c_matrix)
cell_types <- gsub("n\\d+h$", "", col_names)
times <- gsub(".*n(\\d+h)$", "\\1", col_names)

annotation_col <- data.frame(CellType = cell_types, Time = times)
rownames(annotation_col) <- col_names

png("Task3_Panel2c_Heatmap.png", width = 14, height = 10, units = "in", res = 300)
pheatmap(c_matrix, annotation_col = annotation_col,
         cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = FALSE,
         main = "Task 3/Fig 2c: Expression Across Cell Types and Time",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

#--------------------------------------------------------------------------------------------------------------------

# Task 4. Reproduce panel 2d: Pathway enrichment heatmap
sheet_d1 <- read_excel(data_file, sheet = "d_1")
d_matrix <- as.matrix(sheet_d1[, -1])
rownames(d_matrix) <- sheet_d1$pathway

png("Task4_Panel2d_Pathway_Heatmap.png", width = 10, height = 8, units = "in", res = 300)
pheatmap(d_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Task 4/Fig 2d: Pathway Enrichment",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Task 5. Reproduce panel 2e: Bubble plot of kinetic regimes

sheet_e <- read_excel(data_file, sheet = "e")
p3e <- ggplot(sheet_e, aes(x = half_life, y = alpha, color = stage, size = count)) +
  geom_point(alpha = 0.7) + scale_size_continuous(range = c(3, 12)) +
  labs(title = "Task 5/Fig 2e: Kinetic Regimes",
       x = "Half Life", y = "Alpha") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Task5_Panel2e_Bubble.png", p3e, width = 10, height = 7, dpi = 300)

#--------------------------------------------------------------------------------------------------------------------
# Task 6: Panel 2f - Stacked proportions
sheet_f <- read_excel(data_file, sheet = "f")
sheet_f_filtered <- subset(sheet_f, stage %in% c("s00h", "s72h"))

p3f <- ggplot(sheet_f_filtered, aes(x = stage, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Plasma" = "#e15759", "B" = "#4e79a7")) +
  ylim(0, 0.3) +
  labs(title = "Task 6/Fig 2f: Cell Proportions at s00h and s72h",
       x = "Time Point", y = "Proportion") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Task6_Panel2f_Barplot.png", p3f, width = 7, height = 6, dpi = 300)

#-----------------------------------------------------------------------------------------------------------------------
# Task 7: Panel 2g - Directed network
sheet_g <- read_excel(data_file, sheet = "g")
colnames(sheet_g)[1] <- "from"

sheet_g_long <- sheet_g %>%
  pivot_longer(cols = -from, names_to = "to", values_to = "weight") %>%
  filter(weight > 0)

g <- graph_from_data_frame(sheet_g_long, directed = TRUE)
E(g)$width <- E(g)$weight * 5

png("Task7_Panel2g_Network.png", width = 10, height = 8, units = "in", res = 300)
set.seed(123)
plot(g, layout = layout_with_fr, edge.arrow.size = 0.5,
     edge.width = E(g)$width, vertex.color = "#4e79a7",
     vertex.size = 30, vertex.label.cex = 1,
     main = "Task 7/Fig 2g: Cell-Cell Interaction Network")
dev.off()

#--------------------------------------------------------------------------------------------------------------------
#Task 8. Final assembly (mandatory)

final_figure1 <- (p1b | p2c) / (p2d | p2e) / (p2f) +
  plot_annotation(title = "HackBio Stage 2 Project - Part 1 & 2",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
ggsave("Final_Figure_Part1_2.png", final_figure1, width = 16, height = 20, dpi = 300)

final_figure2 <- (p3a | p3b) / (p3e | p3f) +
  plot_annotation(title = "HackBio Stage 2 Project - Task Panels",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
ggsave("Final_Figure_Tasks.png", final_figure2, width = 16, height = 12, dpi = 300)

cat("\n✅ Project complete! All files saved.\n")

#This uses the patchwork package.
#| means "place side-by-side."
#/ means "place underneath."
#It combines 5 separate plots into one single master image.
#cat stands for "concatenate and print." It simply prints a success message to your R console.