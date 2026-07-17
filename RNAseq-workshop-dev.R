# Load packages ----------------------------------------
setwd("~/OneDrive - University of Otago/GenomicsAotearoa/Workshops_GA/testing-new-content/RNAseq-new-content/rnaseq-spotty/")

library(ggplot2)
library(tidyr)
library(dplyr)
library(DESeq2)
library(here)

#  Data prep ------------------------------------------

counts.table <- read.table(here("data-files/counts_OG.txt"),
                           header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, row.names = "GeneID")

coldata.table <- read.table(here("data-files/coldata_OG.txt"),
                            header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)


coldata.table$batchname <- as.factor(coldata.table$batchname)
coldata.table$histov3 <- as.factor(coldata.table$histov3)


dds <- DESeqDataSetFromMatrix(countData = counts.table,
                              colData   = coldata.table,
                              design    = ~ histov3 + batchname)
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
                           by.x = "Sample", by.y = "sample")

custom_order  <- c("F", "ET", "MT", "LT", "TPM", "IPM")
custom_colors <- c("#F564E3", "#984EA3", "#00BFC4",
                   "#4DAF4A", "#619CFF",  "#CD9600")

data.long.coldata$histov3 <- factor(data.long.coldata$histov3,
                                     levels = custom_order)


# Plotting ------------------------------

pca_data <- plotPCA(vsd, intgroup = c("batchname","histov3"), returnData = T)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

pca_data$histov3 <- factor(pca_data$histov3,
                            levels = custom_order)


ggplot(pca_data, aes(x = PC1, y = PC2, colour = histov3, shape = batchname)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  scale_colour_manual(values = custom_colors) +
  scale_shape_manual(values = c( 16, 17, 15, 18)) +
  theme_bw()


table(coldata.table$batchname, coldata.table$histov3)


# sample choice ----

subset_F_MT_TPM_check <-  coldata.table[coldata.table$histov3 %in% c("F", "MT", "TPM"),]
 
table(subset_F_MT_TPM_check$histov3)
table(coldata.table$histov3)
subset_F_MT_TPM_check$histov3 <- droplevels(subset_F_MT_TPM_check$histov3)
table(subset_F_MT_TPM_check$histov3)

subset_F_MT_TPM_check |> 
  filter(histov3 == "TPM") |> 
  select(sample, histov3)


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

# 5 sample subsets

coldata.subset_F_MT_TPM <- subset_F_MT_TPM_check |> 
  filter(sample %in% c(
    # F
    "SI16_21G", "SI16_23G", "SI18_30G", "SI18_48G", "SI18_52G",
    # MT
    "SI16_22G", "SI18_15G", "SI18_16G", "SI18_18G", "SI18_22G",
    # TPM
    "SI16_1G", "SI16_2G", "SI16_3G", "SI18_35G", "SI18_36G"
  ))

table(coldata.subset_F_MT_TPM$histov3)



counts.subset_F_MT_TPM <- counts.table |> 
  select(#F
    "SI16_21G", "SI16_23G", "SI18_30G", "SI18_48G", "SI18_52G",
    # MT
    "SI16_22G", "SI18_15G", "SI18_16G", "SI18_18G", "SI18_22G",
    # TPM
    "SI16_1G", "SI16_2G", "SI16_3G", "SI18_35G", "SI18_36G" )

# tests
nrow(counts.subset_F_MT_TPM) == nrow(counts.table)
ncol(counts.subset_F_MT_TPM) == nrow(coldata.subset_F_MT_TPM)

head(counts.subset_F_MT_TPM)


dir()

write.table(coldata.subset_F_MT_TPM, "data-files/coldata.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(counts.subset_F_MT_TPM, "data-files/counts.txt", sep="\t", quote=FALSE, row.names=TRUE)

#wrote back out originals to unbreak script at start
#write.table(coldata.table, "data-files/coldata_OG.txt", sep="\t", quote=FALSE, row.names=TRUE)
#write.table(counts.table, "data-files/counts_OG.txt", sep="\t", quote=FALSE, row.names=TRUE)



# test read back in

coldatatest <- read.table("data-files/coldata.txt", sep="\t", header=TRUE)
countstest  <- read.table("data-files/counts.txt",  sep="\t", header=TRUE, row.names=1)

# seems ok as above
rm(coldatatest)
rm(countstest)




library(ggplot2)
library(dplyr)

# Convert counts data to long format for ggplot
counts_long <- counts.subset_F_MT_TPM %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(cols = -Gene, names_to = "sample", values_to = "Counts")

# Create boxplot with ggplot2
ggplot(counts_long, aes(x = sample, y = log(Counts + 1))) +
  geom_boxplot() +
  labs(x = "Samples", y = "Counts") +
  theme_minimal()


counts_long_metadata <- counts_long |> 
  left_join(coldata.subset_F_MT_TPM, by = "sample")




ggplot(counts_long_metadata, aes(x = sample, y = log(Counts + 1), fill = histov3)) +
  geom_boxplot() +
  labs(x = "Samples", y = "Counts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 


# need to reload coldata subset in as have chnage colnames now

coldata.subset_F_MT_TPM <- read.table(here("data-files/coldata.txt"),
                                     header = TRUE, sep = "\t",
                                     stringsAsFactors = FALSE)

library(tidyverse)
colSums(counts.subset_F_MT_TPM) %>%
  as.data.frame() %>%
  setNames("reads") |> 
  tibble::rownames_to_column("sample") %>%
  left_join(coldata.subset_F_MT_TPM, by = "sample") %>%
  mutate(sample = fct_reorder(sample, histology)) %>% #fctreorder comes from forcats
  ggplot(aes(x = sample, y = reads, fill = histology)) +
  geom_col(colour = "black") +
  labs(
    y = "Reads mapped per sample",
    x = "Sample",
    title = "Library size") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# could be useful for data manip w tidyverse
# testing if there are differences in average library size by histology group. also exploring variation for each sample in pct of the mean. 
colSums(counts.subset_F_MT_TPM) %>%
  as.data.frame() %>%
  setNames("reads") |>
  tibble::rownames_to_column("sample") |> 
  left_join(coldata.subset_F_MT_TPM, by = "sample") |> 
  mutate(pct = (reads / mean(reads)) * 100) |> 
  group_by(histology) |> 
  mutate(mean.group = mean(reads)) |> 
  select(sample, reads, histology, pct, mean.group)

# count boxplots

#create new df that includes counts in long form and all metadata in the coldata
counts_long_metadata <- counts_long |> 
  left_join(coldata.subset_F_MT_TPM, by = "sample")  |>
  mutate(sample = fct_reorder(sample, histology))

# Create log counts boxplot with ggplot2
ggplot(counts_long_metadata, aes(
  x = sample,
  y = log2(Counts + 1),
  fill = histology)) +
  geom_boxplot() +
  labs(x = "Samples", y = "Log Counts +1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# heatmap

ggplot(counts_long_metadata, aes(
  x = sample,
  y = Gene,
  fill = log2(Counts + 1))) +
  geom_tile() +
  labs(x = "Samples", y = "Log Counts +1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# top 1000 genes HMP
#top1000 <- counts_long_metadata %>%
#  group_by(Gene) %>%
#  summarise(mean_counts = mean(Counts)) %>%
#  slice_max(mean_counts, n = 1000) %>%
#  pull(Gene)
#
#counts_long_metadata %>%
#  filter(Gene %in% top1000) %>%
#  ggplot(aes(x = sample, y = Gene, fill = log(Counts + 1))) +
#  geom_tile() +
#  labs(x = "Samples", y = "Gene") +
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1),
#        axis.text.y = element_blank())
#


counts_filtered <- counts.subset_F_MT_TPM[rowSums(counts.subset_F_MT_TPM) > 0, ]

pca <- prcomp(t(log2(counts_filtered + 1)))

pca_data <- as.data.frame(pca$x) %>%
  tibble::rownames_to_column("sample") %>%
  left_join(coldata.subset_F_MT_TPM, by = "sample")

pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2),1)

ggplot(pca_data, aes(x = PC1, y = PC2, colour = histology, shape = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample), nudge_y = 15, size = 3) +
  labs(x = paste0("PC1: ", pct_var[1], "% variance"),
       y = paste0("PC2: ", pct_var[2], "% variance")) +
  theme_classic()


# test out using different samples instead ----

coldata.table |> 
  filter(batchname == "SI18") |> 
  filter(histov3 == "TPM")

#SI18_40G
#SI18_41G
#SI18_47G
#SI18_49G

counts.test <- counts.subset_F_MT_TPM
counts.test <- cbind(counts.test, counts.table[, c("SI18_47G", "SI18_49G")])

ncol(counts.test) #17
counts.test <- counts.test[, !colnames(counts.test) %in% c("SI18_35G", "SI18_36G")]
ncol(counts.test) #15


counts_filtered_test <- counts.test[rowSums(counts.test) > 0, ]

pca_test <- prcomp(t(log2(counts_filtered_test + 1)))

pca_data_test <- as.data.frame(pca_test$x) %>%
  tibble::rownames_to_column("sample") %>%
  left_join(coldata.table, by = "sample")

pct_var_test <- round(100 * pca_test$sdev^2 / sum(pca_test$sdev^2),1)

ggplot(pca_data_test, aes(x = PC1, y = PC2, colour = histov3, shape = batchname)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample), nudge_y = 15, size = 3) +
  labs(x = paste0("PC1: ", pct_var_test[1], "% variance"),
       y = paste0("PC2: ", pct_var_test[2], "% variance")) +
  theme_classic()

colnames(counts.test)

#[1] "SI16_21G" "SI16_23G" "SI18_30G" "SI18_48G" "SI18_52G" "SI16_22G" "SI18_15G" "SI18_16G"
#[9] "SI18_18G" "SI18_22G" "SI16_1G"  "SI16_2G"  "SI16_3G"  "SI18_47G" "SI18_49G"

# output these as new subsetted files for workshop

write.table(counts.test, "data-files/counts.txt", sep="\t", quote=FALSE, row.names=TRUE)
coldata_subset_test <- coldata.table |> 
  filter(sample %in% c("SI16_21G", "SI16_23G", "SI18_30G", "SI18_48G", "SI18_52G", "SI16_22G", "SI18_15G", "SI18_16G", "SI18_18G", "SI18_22G", "SI16_1G", "SI16_2G", "SI16_3G", "SI18_47G", "SI18_49G")) 


write.table(coldata_subset_test, "data-files/coldata.txt", sep="\t", quote=FALSE, row.names=TRUE)

# limma ----


coldata <- read.table("data-files/coldata.txt", sep="\t", header=TRUE)
counts <- read.table("data-files/counts.txt",  sep="\t", header=TRUE, row.names=1)

library(limma)
library(edgeR)
library(dplyr)

dge <- DGEList(counts=counts)

dge <- normLibSizes(dge)



logCPM <- cpm(dge, log=TRUE, prior.count=3)


head(logCPM, 3)

coldata$histology <- factor(coldata$histology, levels = c("F", "MT", "TPM"))
coldata$batchname <- factor(coldata$batchname)


# design and custom contrasts ----

coldata <- coldata[match(colnames(counts), coldata$sample), ]
# coldata sample order MUST match the same order in counts, or the dge object and deisng matrix wont match each other!

design <- model.matrix(~0 + histology + batchname, data = coldata)
design

v <- voom(dge, design, plot = TRUE)

# 1. Fit the linear model for each gene
fit <- lmFit(v, design)


contrasts.matrix <- makeContrasts(
  MT_vs_F   = histologyMT - histologyF,
  TPM_vs_F  = histologyTPM - histologyF,
  TPM_vs_MT = histologyTPM - histologyMT,
  levels    = colnames(design)
)

# 3. Fit the contrasts to the model
fitC <- contrasts.fit(fit, contrasts = contrasts.matrix)

fitC <- contrasts.fit(fit, contrasts = contrasts.matrix)
fitC <- eBayes(fitC)



tt_TPM_vs_F <- topTable(fitC, coef="TPM_vs_F", n=nrow(counts))

head(tt_TPM_vs_F)
write.table(tt_TPM_vs_F, file = "tt_TPM_vs_F.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(tt_MT_vs_F, file = "tt_MT_vs_F.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(tt_TPM_vs_MT, file = "tt_TPM_vs_MT.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


sigGenesLimma_TPM_vs_F <- which(tt_TPM_vs_F$adj.P.Val <= 0.05 & tt_TPM_vs_F$logFC > 1)

sigGenesLimma_TPM_vs_F <- tt[sigGenesLimma_TPM_vs_F, ]

nrow(sigGenesLimma_TPM_vs_F )


# absolute vs only up and only down with test 
sigGenesLimma_allDE <- which(tt_TPM_vs_F$adj.P.Val <= 0.05 & abs(tt_TPM_vs_F$logFC) > 1)
sigGenesLimma_allDE <- tt_TPM_vs_F[sigGenesLimma_allDE, ]


sigGenesLimma_up <- which(tt_TPM_vs_F$adj.P.Val <= 0.05 & tt_TPM_vs_F$logFC > 1)
sigGenesLimma_up <- tt_TPM_vs_F[sigGenesLimma_up, ]

sigGenesLimma_down <- which(tt_TPM_vs_F$adj.P.Val <= 0.05 & tt_TPM_vs_F$logFC < -1)
sigGenesLimma_down <- tt_TPM_vs_F[sigGenesLimma_down, ]

# test
nrow(sigGenesLimma_allDE) == nrow(sigGenesLimma_up) + nrow(sigGenesLimma_down)


# deseq2 ----

library(DESeq2)

coldata_nofactor <- coldata %>% mutate(across(where(is.factor), as.character))

coldata$histology <- factor(coldata$histology, levels = c("F", "MT", "TPM"))
coldata$batchname <- factor(coldata$batchname)

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = coldata, 
                              design = ~histology + batchname)



dds <- DESeq(dds)

dds |> str()
resultsNames(dds)


res_MT_vs_F <- results(dds, name = "histology_MT_vs_F") |> 
  na.omit()

res_MT_vs_F %>% head()


res_MT_vs_TPM <- results(dds, contrast = c("histology", "MT", "TPM")) |> na.omit()
res_MT_vs_TPM %>% head()



# deseq2 tpm up vs F ----


res_F_vs_TPM <- results(dds, contrast = c("histology", "F", "TPM")) 

res_F_vs_TPM |> head()

res_F_vs_TPM <- res_F_vs_TPM |> na.omit()

# Keep all rows in the res object if the adjusted p-value < 0.05
resPadj <- res_F_vs_TPM[res_F_vs_TPM$padj <= 0.05 , ]

# Keep all rows in the res object if the adjusted p-value < 0.05 AND the log2 fold change is less than -1. 
resPadjLogFC_TPM_vs_F <- res_F_vs_TPM[res_F_vs_TPM$padj <= 0.05 & res_F_vs_TPM$log2FoldChange < -1,]

resPadjLogFC_TPM_vs_F  |> dim() 

resPadjLogFC_TPM_vs_F |>  head()


# venn ----


sigGenesLimma_TPM_vs_F |> nrow()
sigGenesDESeq_TPM_vs_F |> nrow()

#BiocManager::install("gplots")
library(gplots)

setlist <- list(Limma = rownames(sigGenesLimma_TPM_vs_F), 
                DESeq2 = rownames(sigGenesDESeq_TPM_vs_F))

venn(setlist)

str(setlist)

# Convert venn data to list object
intersect_list <- attr(venn(setlist), "intersections")
intersect_list  |> str()

# Create vectors of the three gene lists 
Limma_only_genes <- intersect_list$Limma
DEseq2_only_genes <- intersect_list$DESeq2
Both_genes <- intersect_list$`Limma:DESeq2`




getwd()
save(res_F_vs_TPM, file = "deseq2-FvsTPM.RData") 
dim(res_F_vs_TPM)

# goseq ----

hypergeoDistMatrix <- matrix(c(10,490,90,9410),2,2)
hypergeoDistMatrix

ftest <- fisher.test(hypergeoDistMatrix)
round(fisher.test(hypergeoDistMatrix)$p.value, 3)


res_F_vs_TPM
dds 
counts |> dim()

# checking effect of na.omit 
res_test<- results(dds, name = "histology_TPM_vs_F") 
res_test


