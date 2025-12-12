# Loading libraries and count data ---------------------------------------------
library(dplyr) 
library(Seurat) 
library(ggplot2)
library(ggrepel)
library(tidyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scran) 
library(scuttle) 
library(tibble)
library(edgeR) 
library(pheatmap) 
library(fgsea)
library(viridis)
library(limma)
library(enrichplot)

counts <- readRDS("objects/merged_counts.rds")
metadata <- readRDS("objects/metadata_updated.rds")
counts$gene_name <- make.unique(counts$gene_name)
counts_mat <- as.matrix(counts[, -(1:2)])
rownames(counts_mat) <- counts$gene_name 
seurat_obj <- readRDS("objects/seurat_normalized.rds")

##Seurat object creation, normalization, removal of batch effects:
# seurat_obj <- CreateSeuratObject(counts=counts_mat, meta.data=metadata)
# meta_all <- seurat_obj@meta.data
# meta_all$nFeature <- seurat_obj$nFeature_RNA
# meta_all$nCount <- seurat_obj$nCount_RNA
# 
# keep_cells <- which(
#   meta_all$nFeature > median(meta_all$nFeature) - 3*mad(meta_all$nFeature) &
#     meta_all$nFeature < median(meta_all$nFeature) + 3*mad(meta_all$nFeature) &
#     meta_all$nCount > median(meta_all$nCount) - 3*mad(meta_all$nCount) &
#     meta_all$nCount < median(meta_all$nCount) + 3*mad(meta_all$nCount)
# )
# 
# meta_all$filter_status <- ifelse(seq_len(nrow(meta_all)) %in% keep_cells, "kept", "removed")
# y_limits <- range(meta_all$nFeature)
# meta_filtered <- subset(meta_all, filter_status == "kept")
# seurat_obj <- subset(seurat_obj, cells = meta_filtered$fiber_id)
# sce <- as.SingleCellExperiment(seurat_obj)
# clusters <- quickCluster(sce)
# sce <- logNormCounts(sce)
# seurat_obj <- SetAssayData(seurat_obj,
#                            slot = "data",
#                            new.data = assay(sce, "logcounts"))
# seurat_obj <- SetAssayData(seurat_obj,
#                            slot = "counts",
#                            new.data = assay(sce, "counts"))
# seurat_obj$size_factors <- sizeFactors(sce)
# seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst")
# DefaultAssay(seurat_obj) <- "RNA"
# seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
# seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
# seurat_obj <- RunHarmony(
#   object = seurat_obj,
#   group.by.vars = c("extra_pcr", "sequencing_batch"),
#   reduction.use = "pca"
# )

# DE condition -----------------------------------------------------------------

pseudo_bulk_CS <- readRDS("objects/pseudo_bulk_CS.rds")
sample_info_CS <- readRDS("objects/sample_info_CS.rds")

##Pseudobulk object creation
# counts <- GetAssayData(seurat_obj, slot="counts")
# meta <- seurat_obj@meta.data
# sep_char <- "__"
# donor_fiber <- paste(meta$sample, meta$fiber_type, sep=sep_char)
# 
# pseudo_bulk_all <- sapply(unique(donor_fiber), function(df){
#   cells <- rownames(meta)[donor_fiber == df]
#   rowSums(counts[, cells, drop=FALSE])
# })
# pseudo_bulk_all <- as.data.frame(pseudo_bulk_all)
# pseudo_bulk_all <- round(pseudo_bulk_all)
# 
# pseudo_bulk_log_all <- log1p(pseudo_bulk_all)
# pseudo_bulk_log_all_t <- t(pseudo_bulk_log_all)
# pseudo_bulk_log_all_t <- pseudo_bulk_log_all_t[, apply(pseudo_bulk_log_all_t, 2, var) != 0]
# 
# pca_res_all <- prcomp(pseudo_bulk_log_all_t, scale.=TRUE)
# pca_df_all <- data.frame(
#   PC1 = pca_res_all$x[,1],
#   PC2 = pca_res_all$x[,2],
#   donor_fiber = rownames(pca_res_all$x),
#   stringsAsFactors = FALSE
# )
# pca_df_all <- pca_df_all %>%
#   mutate(
#     donor = sub(paste0(sep_char, ".*$"), "", donor_fiber),
#     fiber_type = sub(paste0("^.*", sep_char), "", donor_fiber)
#   )
# donor_meta <- unique(meta[, c("sample","group")])
# pca_df_all <- left_join(pca_df_all, donor_meta, by = c("donor" = "sample"))

dge <- DGEList(counts = pseudo_bulk_CS)
dge <- calcNormFactors(dge)

sample_info_CS$group <- factor(sample_info_CS$group, levels = c("C", "S"))
design <- model.matrix(~ group, data = sample_info_CS)

v <- voom(dge, design = design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

res_voom <- topTable(fit, coef = "groupS", number = Inf, sort.by = "P")
res_voom$FDR <- p.adjust(res_voom$P.Value, method = "fdr")

expr_mat <- v$E 
pca_res <- prcomp(t(expr_mat), scale. = TRUE)
var_expl <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

volcano_df <- res_voom %>%
  mutate(
    gene = rownames(res_voom),
    negLogFDR = -log10(adj.P.Val),
    sig = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "sig", "ns")
  )

ggplot(volcano_df, aes(x = logFC, y = negLogFDR)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("significance" = "red", "ns" = "grey")) +
  geom_text_repel(
    data = subset(volcano_df, sig == "sig"),
    aes(label = gene),
    size = 3,
    max.overlaps = 20
  ) +
  labs(
    title = "Volcano — ICU-AW vs Control",
    x = "log2 Fold Change (ICU / C)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()

# Pathway analysis CvS ---------------------------------------------------------

ranked_genes_CS <- res_voom$logFC
names(ranked_genes_CS) <- rownames(res_voom)
ranked_genes_CS <- ranked_genes_CS[!is.na(ranked_genes_CS)]
ranked_genes_CS <- sort(ranked_genes_CS, decreasing = TRUE)

msig <- readRDS("objects/msig_hallmark.rds")
pathways <- split(msig$gene_symbol, msig$gs_name)

fgsea_CS <- fgsea(pathways = pathways, stats = ranked_genes_CS, nperm = 10000)
fgsea_CS <- fgsea_CS[order(fgsea_CS$padj), ]
fgsea_CS$leading_edge_count <- sapply(fgsea_CS$leadingEdge, length)

fgsea_CS %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = NES, y = reorder(pathway, NES), fill = padj)) +
  geom_col() +
  geom_text(
    aes(label = leading_edge_count),
    hjust = -0.1,
    size = 3.5
  ) +
  scale_fill_viridis_c(direction = -1, option = "C") +
  theme_minimal() +
  labs(
    title = "Top Enriched Hallmark Pathways (ICU-AW vs Control)",
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    fill = "FDR"
  ) +
  coord_cartesian(xlim = c(min(fgsea_CS$NES), max(fgsea_CS$NES) * 1.15))

pathway_names <- names(pathways)

jaccard_matrix <- matrix(0, nrow = length(pathways), ncol = length(pathways),
                         dimnames = list(pathway_names, pathway_names))

for (i in seq_along(pathways)) {
  for (j in seq_along(pathways)) {
    A <- pathways[[i]]
    B <- pathways[[j]]
    jaccard_matrix[i, j] <- length(intersect(A, B)) / length(union(A, B))
  }
}

intersect_count <- matrix(0, nrow = length(pathways), ncol = length(pathways),
                          dimnames = list(pathway_names, pathway_names))

for (i in seq_along(pathways)) {
  for (j in seq_along(pathways)) {
    A <- pathways[[i]]
    B <- pathways[[j]]
    intersect_count[i, j] <- length(intersect(A, B))
  }
}

pheatmap(
  jaccard_matrix,
  fontsize_row = 6,
  fontsize_col = 6,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Pathway Jaccard Overlap"
)

# Top nominal gene exploration -------------------------------------------------

top_n <- 50  
top_genes_CS <- rownames(res_voom[order(res_voom$adj.P.Val), ])[1:top_n]

logcounts_mat <- as.matrix(GetAssayData(seurat_obj, slot = "data"))
logcounts_top <- logcounts_mat[top_genes_CS, , drop = FALSE]

logcounts_top_scaled <- t(scale(t(logcounts_top)))
colnames(logcounts_top_scaled) <- colnames(logcounts_top)  

cell_meta <- seurat_obj@meta.data
annotation_col <- cell_meta %>%
  dplyr::select(group, sample)
rownames(annotation_col) <- colnames(logcounts_top_scaled)
cell_meta <- seurat_obj@meta.data
annotation_col <- cell_meta %>%
  dplyr::select(group, sample) %>%
  dplyr::rename(donor = sample) %>%  
  dplyr::rename(condition = group) %>%        
  dplyr::mutate(condition = recode(condition,    
                                   "C" = "Control",
                                   "S" = "ICU"))
rownames(annotation_col) <- colnames(logcounts_top_scaled)

heat <- pheatmap(
  logcounts_top_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  scale = "row",
  fontsize = 10,
  main = "Per-fiber expression of top 50 DE genes (lowest FDR)",
  silent = TRUE
)

# Select cluster ---------------------------------------------------------------

##Locate cluster
# col_clusters <- cutree(heat$tree_col, k = 10)
# annotation_col_with_cluster <- annotation_col %>%
#   dplyr::mutate(sample = rownames(.),
#                 cluster = factor(col_clusters[sample])) %>%
#   dplyr::select(-sample)
# 
# cluster_colors <- rainbow(n_col_clusters)
# annotation_colors <- list(cluster = setNames(cluster_colors, 1:n_col_clusters))
# pheatmap(
#   logcounts_top_scaled,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   annotation_col = annotation_col_with_cluster,
#   annotation_colors = annotation_colors,
#   show_rownames = TRUE,
#   show_colnames = FALSE,
#   scale = "row",
#   fontsize = 10,
#   main = "Per-fiber expression of top 50 DE genes (lowest FDR)"
# )

rownames(seurat_obj@meta.data) <- colnames(logcounts_top_scaled)
branches <- cutree(heat$tree_col, k = 10)

cells_branch <- names(branches)[branches %in% c(3, 10)]

rownames(seurat_obj@meta.data) <- seurat_obj$fiber_id

# Add metadata column for cluster status
seurat_obj$in_cluster <- ifelse(rownames(seurat_obj@meta.data) %in% cells_branch,
                                "cluster", "other")
Idents(seurat_obj) <- "in_cluster"

meta <- seurat_obj@meta.data %>%
  dplyr::mutate(
    in_cluster = factor(in_cluster, levels = c("other", "cluster")),
    condition = recode(group, "C" = "Control", "S" = "ICU-AW") %>%
      factor(levels = c("Control", "ICU-AW"))
  )

plot_data <- meta %>%
  group_by(sample, in_cluster, condition) %>%
  summarise(n = n(), .groups = "drop")

ggplot(plot_data, aes(x = sample, y = n, fill = in_cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("other" = "gray70", "cluster" = "red")) +
  facet_wrap(~condition, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Fibers per donor by cluster status",
    x = "Donor",
    y = "Number of fibers",
    fill = "Cluster status"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# DE cluster v other -----------------------------------------------------------

cluster_cells <- rownames(seurat_obj@meta.data)[seurat_obj$in_cluster == "cluster"]
Idents(seurat_obj) <- factor(ifelse(colnames(seurat_obj) %in% cluster_cells, "cluster", "other"))

counts <- GetAssayData(seurat_obj, slot = "counts")
meta <- seurat_obj@meta.data
all_donors <- unique(meta$sample)

pseudo_bulk_cluster <- sapply(all_donors, function(d) {
  cells <- rownames(meta)[meta$sample == d & meta$in_cluster == "cluster"]
  if(length(cells) == 0) {
    rep(0, nrow(counts))   # placeholder for donors without cluster cells
  } else {
    rowSums(counts[, cells, drop = FALSE])
  }
})

pseudo_bulk_other <- sapply(all_donors, function(d) {
  cells <- rownames(meta)[meta$sample == d & meta$in_cluster == "other"]
  if(length(cells) == 0) {
    rep(0, nrow(counts))
  } else {
    rowSums(counts[, cells, drop = FALSE])
  }
})

pseudo_bulk <- cbind(pseudo_bulk_cluster, pseudo_bulk_other)
colnames(pseudo_bulk) <- paste0(
  rep(all_donors, 2), "_", rep(c("cluster","other"), each = length(all_donors))
)

sample_info <- data.frame(
  sample = colnames(pseudo_bulk),
  donor = rep(all_donors, 2),
  cluster_status = rep(c("cluster","other"), each = length(all_donors))
)
sample_info$cluster_status <- factor(sample_info$cluster_status, levels = c("other","cluster"))
rownames(sample_info) <- sample_info$sample

lib_sizes <- colSums(pseudo_bulk)
zero_lib <- lib_sizes == 0
pseudo_bulk_filtered <- pseudo_bulk[, !zero_lib]
sample_info_filtered <- sample_info[colnames(pseudo_bulk_filtered), ]

cat("Remaining samples for DGEList:", ncol(pseudo_bulk_filtered), "\n")

dge <- DGEList(counts = pseudo_bulk_filtered)
dge <- calcNormFactors(dge)

design <- model.matrix(~ cluster_status, data = sample_info_filtered)
v <- voom(dge, design, plot = TRUE)

corfit <- duplicateCorrelation(v, design, block = sample_info_filtered$donor)
cat("Consensus correlation:", corfit$consensus.correlation, "\n")

fit <- lmFit(v, design, block = sample_info_filtered$donor, correlation = corfit$consensus.correlation)
fit <- eBayes(fit)

top_genes <- topTable(fit, coef = "cluster_statuscluster", number = Inf, sort.by = "P")
top_genes$FDR <- p.adjust(top_genes$P.Value, method = "fdr")
top_genes$negLogFDR <- -log10(top_genes$FDR)
top_genes$sig <- ifelse(top_genes$FDR < 0.05 & abs(top_genes$logFC) > 0.5, "sig", "ns")

top_genes$gene <- rownames(top_genes)

top_label_genes <- top_genes$gene[order(top_genes$FDR)][1:50]

ggplot(top_genes, aes(x = logFC, y = negLogFDR, color = sig)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = subset(top_genes, gene %in% top_label_genes),
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("sig" = "red", "ns" = "grey")) +
  labs(title = "Volcano Plot — Cluster vs Other", x = "log2 Fold Change", y = "-log10(FDR)")

# Cluster DEG annotation -------------------------------------------------------

# GO term 
top_sig <- rownames(top_genes)[1:108]

ego <- enrichGO(gene = top_sig,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "fdr",
                readable = TRUE)
ego_clusters <- pairwise_termsim(ego)

treeplot(ego_clusters, showCategory = 30) +
  ggtitle("GO Term Treeplot for Significant Cluster DEGs") +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  )

# Hallmark

ranked_genes <- top_genes$logFC
names(ranked_genes) <- top_genes$gene
ranked_genes <- ranked_genes[!is.na(ranked_genes)]
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

fgsea_res <- fgsea(pathways = pathways, stats = ranked_genes, nperm = 10000)
fgsea_res <- fgsea_res[order(fgsea_res$padj), ]
fgsea_res$leading_edge_count <- sapply(fgsea_res$leadingEdge, length)
fgsea_sig <- fgsea_res %>%
  filter(padj < 0.05) %>%
  arrange(padj)
fgsea_sig %>%
  slice_head(n = 20) %>%  # optional: top 20 by padj
  ggplot(aes(x = NES, y = reorder(pathway, NES), fill = padj)) +
  geom_col() +
  geom_text(aes(label = leading_edge_count), hjust = -0.1, size = 3.5) +
  scale_fill_viridis_c(direction = -1, option = "C") +
  theme_minimal() +
  labs(
    title = "Top Enriched Hallmark Pathways: Cluster DEGs",
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    fill = "FDR"
  ) +
  coord_cartesian(xlim = c(min(fgsea_sig$NES), max(fgsea_sig$NES) * 1.15))
