# Load packages ----------------------------------------
setwd("~/OneDrive - University of Otago/GenomicsAotearoa/Workshops_GA/testing-new-content/RNAseq-new-content/rnaseq-spotty/")

library(ggplot2)
library(tidyr)
library(dplyr)
library(DESeq2)
library(here)

#  Data prep ------------------------------------------

counts.table <- read.table(here("data-files/readspergene_v2.matrix"),
                           header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, row.names = "GeneID")

coldata.table <- read.table(here("data-files/coldata_v2.3.2.txt"),
                            header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, row.names = "sample")


coldata.table$subsetname <- as.factor(coldata.table$subsetname)
coldata.table$histov22 <- as.factor(coldata.table$histov22)


dds <- DESeqDataSetFromMatrix(countData = counts.table,
                              colData   = coldata.table,
                              design    = ~ histov22 + subsetname)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

vsd_normalized_counts <- as.data.frame(assay(vsd))
vsd_normalized_counts$GeneID <- rownames(vsd_normalized_counts)
rownames(vsd_normalized_counts) <- NULL
vsd_normalized_counts <- vsd_normalized_counts[, c(97, 1:96)]

data.long <- pivot_longer(vsd_normalized_counts,
                          cols = -GeneID,
                          names_to  = "Sample",
                          values_to = "expression")

data.long.coldata <- merge(data.long, coldata.table,
                           by.x = "Sample", by.y = "samplenames")

custom_order  <- c("F", "ET", "MT", "LT", "TPM", "IPM")
custom_colors <- c("#F564E3", "#984EA3", "#00BFC4",
                   "#4DAF4A", "#619CFF",  "#CD9600")

data.long.coldata$histov22 <- factor(data.long.coldata$histov22,
                                     levels = custom_order)


# Plotting ------------------------------

pca_data <- plotPCA(vsd, intgroup = c("subsetname","histov22"), returnData = T)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

pca_data$histov22 <- factor(pca_data$histov22,
                            levels = custom_order)


ggplot(pca_data, aes(x = PC1, y = PC2, colour = histov22, shape = subsetname)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  scale_colour_manual(values = custom_colors) +
  scale_shape_manual(values = c( 16, 17, 15, 18)) +
  theme_bw()


table(coldata.table$subsetname, coldata.table$histov22)


# sample choice ----

subset_F_MT_TPM_check <-  coldata.table[coldata.table$histov22 %in% c("F", "MT", "TPM"),]
 
table(subset_F_MT_TPM_check$histov22)
table(coldata.table$histov22)
subset_F_MT_TPM_check$histov22 <- droplevels(subset_F_MT_TPM_check$histov22)
table(subset_F_MT_TPM_check$histov22)

subset_F_MT_TPM_check |> 
  filter(histov22 == "TPM") |> 
  select(samplenames, histov22)


# F
# SI16_21G
# SI16_23G
# SI18_30G
# SI18_48G
# SI18_52G

# MT
# SI16_22G 
# SI18_15G
# SI18_16G
# SI18_18G
# SI18_22G

# TPM
# SI16_1G
# SI16_2G
# SI16_3G
# SI18_35G
# SI18_36G

subset_F_MT_TPM <- subset_F_MT_TPM_check |> 
  filter(samplenames %in% c(
    # F
    "SI16_21G", "SI16_23G", "SI18_30G", "SI18_48G", "SI18_52G",
    # MT
    "SI16_22G", "SI18_15G", "SI18_16G", "SI18_18G", "SI18_22G",
    # TPM
    "SI16_1G", "SI16_2G", "SI16_3G", "SI18_35G", "SI18_36G"
  ))

table(subset_F_MT_TPM$histov22)





