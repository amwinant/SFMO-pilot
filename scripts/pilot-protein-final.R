# Load packages and data -------------------------------------------------------
library(vroom)
library(here)
library(PhosR)
library(msigdbr)
library(readr)
library(gt)
library(edgeR)

data_raw <- vroom(here("objects/sepsis_pilot_diann_report.unique_genes_matrix.tsv"))
metadata <- read.csv(here("objects/metadata_protein.csv")) |>
  mutate(fiber_id = Subject.ID, subject = gsub(pattern = "f.*", replacement = "", Subject.ID)
  ) |>
  dplyr::select(!Subject.ID) |>
  mutate(condition = case_when(stringr::str_starts(subject, "s") ~ "s", stringr::str_starts(subject, "c") ~ "c", TRUE ~ "error")
  )
# Data normalizing and filtering -----------------------------------------------

data <- data_raw |>
  tibble::column_to_rownames("Genes") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Analytical.Sample.ID") |>
  dplyr::mutate(
    Analytical.Sample.ID = gsub("_S2-.*", "", Analytical.Sample.ID),
    Analytical.Sample.ID = gsub(".*_", "", Analytical.Sample.ID)
  ) |>
  dplyr::inner_join(
    metadata |> dplyr::select(Analytical.Sample.ID, fiber_id),
    by = "Analytical.Sample.ID"
  ) |>
  dplyr::select(!Analytical.Sample.ID) |>
  tibble::column_to_rownames("fiber_id") |>
  t() |>
  as.data.frame()

transformed_data <- log2(data)

condition_vec <- metadata$condition[match(colnames(transformed_data), metadata$fiber_id)]

stopifnot(!any(is.na(condition_vec)))  

transformed_data |>
  pivot_longer(cols = everything(), names_to = "sample_id", values_to = "intensities") |>
  ggplot(aes(x = sample_id, y = intensities)) +
  geom_boxplot(outlier.size = 0.25) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(face = "bold", colour = "black", size = 6)
  )

filtered_data <- PhosR::selectGrps(
  mat = transformed_data,
  percent = 0.7,
  grps = condition_vec
)

normalized_data <- limma::normalizeBetweenArrays(filtered_data, method = "scale") |> as.data.frame()

normalized_data |>
  pivot_longer(cols = everything(), names_to = "sample_id", values_to = "intensities") |>
  ggplot(aes(x = sample_id, y = intensities)) +
  geom_boxplot(outlier.size = 0.25) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(face = "bold", colour = "black", size = 6)
  )

filtering_columns_Na <- function(.data, percentage_accepted_missing) {
  keep_vector <- .data |>
    is.na() |>
    colSums()
  filtering_vector <- keep_vector / nrow(.data)
  filtering_vector <- filtering_vector <= percentage_accepted_missing
  vector_id <- data.frame(
    ID = colnames(.data),
    keep = filtering_vector) |>
    dplyr::filter(keep == TRUE)
  data_filtered <- .data |>
    dplyr::select(vector_id$ID)
  return(data_filtered)
}

filtered_data <- transformed_data |>
  filtering_columns_Na(percentage_accepted_missing = 0.5)
filtered_metadata <- metadata %>%
  filter(fiber_id %in% colnames(filtered_data)) %>%
  arrange(match(fiber_id, colnames(filtered_data)))

stopifnot(identical(filtered_metadata$fiber_id, colnames(filtered_data)))

filtered_data <- filtered_data |>
  PhosR::selectGrps(grps = filtered_metadata$condition, percent = 0.7)

normalized_data <- filtered_data |>
  limma::normalizeBetweenArrays(method = "scale") |>
  as.data.frame()

# Seurat object ----------------------------------------------------------------
seurat_object <- normalized_data |>
  tImpute(m = 1.8, s = 0.3) |>
  as.data.frame()

seurat_object <- Seurat::CreateSeuratObject(counts = seurat_object,
                                            project = "sepsis",
                                            meta.data = filtered_metadata)
# Seurat processing of filtered data --------------------------------------

seurat_object <- Seurat::FindVariableFeatures(seurat_object,
                                              selection.method = "vst")

seurat_object[["RNA"]]$data <- seurat_object[["RNA"]]$counts

seurat_object <- Seurat::ScaleData(seurat_object)

seurat_object <- Seurat::RunPCA(seurat_object, features = Seurat::VariableFeatures(object = seurat_object))

seurat_object <- Seurat::FindNeighbors(seurat_object, dims = 1:10)

pca_df <- data.frame(
  PC1 = Embeddings(seurat_object[["pca"]])[,1],
  PC2 = Embeddings(seurat_object[["pca"]])[,2],
  seurat_object@meta.data
)

var_expl <- (seurat_object[["pca"]]@stdev^2) / sum(seurat_object[["pca"]]@stdev^2) * 100

# DE CvS -----------------------------------------------------------------------

mat_counts <- as.matrix(GetAssayData(seurat_object, slot = "counts")) 

lin_mat <- 2^as.matrix(GetAssayData(seurat_object, slot = "counts")) - 1 

donor_vec <- seurat_object$subject[colnames(lin_mat)] 
donor_levels <- unique(donor_vec) 
pseudobulk_donor <- sapply(donor_levels, function(d) { rowSums(lin_mat[, donor_vec == d, drop = FALSE], na.rm = TRUE) }) 
rownames(pseudobulk_donor) <- rownames(lin_mat) 
colnames(pseudobulk_donor) <- donor_levels 

metadata$donor_id <- sub("f.*", "", metadata$fiber_id) 

donor_meta <- metadata |> dplyr::select(donor_id, condition) |> distinct() 
donor_meta$condition <- factor(donor_meta$condition, levels = c("c", "s")) 
donor_meta <- donor_meta[match(colnames(pseudobulk_donor), donor_meta$donor_id), ] 
stopifnot(all(colnames(pseudobulk_donor) == donor_meta$donor_id)) 

dge <- DGEList(counts = pseudobulk_donor) 
dge <- calcNormFactors(dge) 
design <- model.matrix(~ 0 + donor_meta$condition) 
colnames(design) <- levels(factor(donor_meta$condition)) 
v <- voom(dge, design) 
fit <- lmFit(v, design) 
contr <- makeContrasts(SvsC = s - c, levels = design) 
fit2 <- contrasts.fit(fit, contr) 
fit2 <- eBayes(fit2) 
DE_results <- topTable(fit2, number = Inf) 
DE_results <- DE_results %>% mutate( Signif_FDR = adj.P.Val < 0.05 & abs(logFC) > 1, Signif_P = P.Value < 0.05 & abs(logFC) > 1 ) 
ggplot(DE_results, aes(x = logFC, y = -log10(adj.P.Val), color = Signif_FDR)) + 
  geom_point(size = 1.5) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + scale_color_manual(values = c("grey", "red")) + 
  labs(title = "Volcano Plot (FDR)", y = "-log10(FDR)") + 
  theme_minimal() 

# DE cluster v other -----------------------------------------------------------

cluster_fibers <- c(
  "c22f17", "s14_2f10", "s14_2f8", "s15_2f12", "s15_2f13",
  "s15_2f2", "s15_2f3", "s15_2f7", "s15_2f9", "s20f1",
  "s20f3", "s21f15", "s23f1", "s23f10", "s23f12", "s23f2",
  "s23f3", "s23f4", "s23f5", "s23f6", "s23f7", "s23f8",
  "s23f9", "s26f1", "s26f12", "s26f3", "s5_2f7", "s5_2f8",
  "s8_2f16", "s8_2f2", "s8_2f3", "s8_2f4"
)

filtered_metadata <- filtered_metadata %>%
  mutate(cluster_status = ifelse(fiber_id %in% cluster_fibers, "cluster", "other"))

mat <- as.matrix(normalized_data)
stopifnot(all(colnames(mat) == filtered_metadata$fiber_id))

donor_cluster <- paste(filtered_metadata$subject, filtered_metadata$cluster_status, sep = "_")
split_idx <- split(seq_along(donor_cluster), donor_cluster)

pseudobulk_cluster <- sapply(names(split_idx), function(d) {
  Matrix::rowSums(mat[, split_idx[[d]], drop = FALSE])
})
pseudobulk_cluster[is.na(pseudobulk_cluster)] <- 0
keep <- rowSums(edgeR::cpm(pseudobulk_cluster) > 1) >= 2  
pseudobulk_cluster <- pseudobulk_cluster[keep, ]
colnames(pseudobulk_cluster) <- names(split_idx)

donor_meta_donor <- donor_meta

donor_meta <- filtered_metadata %>%
  dplyr::select(subject, cluster_status) %>%
  distinct() %>%
  dplyr::slice(match(colnames(pseudobulk_cluster), paste(subject, cluster_status, sep = "_")))

stopifnot(identical(paste(donor_meta$subject, donor_meta$cluster_status, sep = "_"),
                    colnames(pseudobulk_cluster)))

design <- model.matrix(~ factor(cluster_status, levels = c("other", "cluster")),
                       data = donor_meta)

dge <- DGEList(counts = pseudobulk_cluster)
dge <- calcNormFactors(dge, method = "TMM")

v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

coef_name <- colnames(design)[2]
results_pseudobulk <- topTable(fit, coef = coef_name, number = Inf)

volcano_df <- results_pseudobulk %>%
  rownames_to_column("protein") %>%
  mutate(
    sig = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "sig", "ns"),
    negLogFDR = -log10(adj.P.Val)
  )

ggplot(volcano_df, aes(x = logFC, y = negLogFDR)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("sig" = "red", "ns" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = subset(volcano_df, sig == "sig"),
    aes(label = protein),
    size = 3,
    max.overlaps = 30
  ) +
  labs(
    title = "Volcano — Cluster vs Other (pseudobulk by donor)",
    x = "log2 Fold Change (cluster / other)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()

# DEP annotation ---------------------------------------------------------------

sig_df <- results_pseudobulk %>%
  rownames_to_column("protein") %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

up_proteins <- sig_df %>% filter(logFC > 0) %>% pull(protein)
down_proteins <- sig_df %>% filter(logFC < 0) %>% pull(protein)

ego_up <- enrichGO(
  gene = up_proteins,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "fdr",
  readable = TRUE
)

ego_up_sim <- pairwise_termsim(ego_up)

ego_down <- enrichGO(
  gene = down_proteins,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "fdr",
  readable = TRUE
)

ego_down_sim <- pairwise_termsim(ego_down)
tree_up <- treeplot(ego_up_sim, showCategory = 30) +
  ggtitle("GO Treeplot — Upregulated Proteins")

ego_down_df <- as.data.frame(ego_down)

down <- ggplot(ego_down_df, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "C", direction = -1) +
  labs(
    title = "GO Terms — Downregulated Proteins",
    x = NULL,
    y = "Gene Count",
    fill = "FDR"
  ) +
  theme_minimal()

tree_up / down


# Compare sig proteins to genes ------------------------------------------------

rna <- read.csv("objects/significant_genes_cluster.csv")
prot <- read.csv("objects/significant_proteins_cluster.csv")
colnames(rna)[1] <- "Gene"
colnames(prot)[1] <- "Protein"

overlap_genes <- intersect(rna$Gene, prot$Protein)
length(overlap_genes) 
overlap_genes

overlap_df <- merge(
  rna %>% filter(Gene %in% overlap_genes) %>% dplyr::select(Gene, logFC_RNA = logFC),
  prot %>% filter(Protein %in% overlap_genes) %>% dplyr::select(Protein, logFC_prot = logFC),
  by.x = "Gene", by.y = "Protein"
)

overlap_df$same_direction <- sign(overlap_df$logFC_RNA) == sign(overlap_df$logFC_prot)

table(overlap_df$same_direction)
xlim <- range(overlap_df$logFC_RNA) * 1.1
ylim <- range(overlap_df$logFC_prot) * 1.1

quadrants <- data.frame(
  xmin = c(0, -Inf, -Inf, 0),
  xmax = c(Inf, 0, 0, Inf),
  ymin = c(0, 0, -Inf, -Inf),
  ymax = c(Inf, Inf, 0, 0),
  fill = c("lightgreen", "orange", "lightgreen", "red")
)

ggplot(overlap_df, aes(x = logFC_RNA, y = logFC_prot, label = Gene)) +

  geom_rect(data = quadrants, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE, alpha = 0.2) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_text_repel(max.overlaps = 20) +
  scale_fill_manual(values = c("lightgreen" = "lightgreen", "orange" = "orange", "red" = "red")) +
  theme_minimal() +
  labs(x = "log2FC RNA", y = "log2FC Protein", title = "RNA vs Protein Fold Change for Overlapping Genes") +
  theme(legend.position = "none")

# Manual annotation 

annotations <- read_csv("objects/coord_annotation.csv")

annotations %>%
  gt() %>%
  tab_header(
    title = "Overlapping feature annotation",
    subtitle = "Manually annotated from NCBI gene summaries"
  ) %>%
  tab_options(
    table.font.size = 13,
    table.border.top.color = "white",
    table.border.bottom.color = "white"
  ) %>%
  cols_width(
    Gene ~ px(120),          
    everything() ~ px(400)  
  ) %>%
  tab_style(
    style = list(
      cell_borders(sides = "bottom", color = "grey80"),
      cell_text(size = 12)
    ),
    locations = cells_body()
  )

<<<<<<< HEAD
# Count data -------------------------------------------------------------------
=======
# Feature counts
>>>>>>> 4db7856bfb54e97078c84aff1a6bd30ef206424a
rna_counts <- GetAssayData(seurat_obj, slot = "counts")
genes_per_cell <- colSums(rna_counts > 0)

rna_df <- data.frame(
  cell = names(genes_per_cell),
  genes = as.numeric(genes_per_cell),
  sample = seurat_obj$sample,
  condition = recode(seurat_obj$group, "C" = "Control", "S" = "ICU-AW")
)

prot_counts <- colSums(!is.na(normalized_data))

prot_df <- data.frame(
  cell = names(prot_counts),
  proteins = prot_counts
) %>%
  left_join(filtered_metadata, by = c("cell" = "fiber_id")) %>%  # metadata has sample/subject info
  mutate(
    condition = recode(condition, "c" = "Control", "s" = "ICU-AW")
  )

prot_df$sample <- gsub("^[csCS]", "", prot_df$subject)
rna_long <- rna_df %>%
  dplyr::select(sample, condition, count = genes) %>%
  mutate(omic = "RNA: Detected Genes")

prot_long <- prot_df %>%
  dplyr::select(sample, condition, count = proteins) %>%
  mutate(omic = "Proteomics: Detected Proteins")

combined_long <- bind_rows(rna_long, prot_long)

rna <- ggplot(rna_df, aes(x = sample, y = genes, fill = sample)) +
  geom_violin(alpha = 0.65) +
  geom_point(size = 1.5, alpha = 0.75) +
  facet_grid(~ condition, scales = "free_x") +
  scale_fill_viridis_d(option = "plasma") +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Detected Genes per Fiber (RNA-seq)",
    x = "Donor",
    y = "Detected Genes"
  )
prot <- ggplot(prot_df, aes(x = sample, y = proteins, fill = sample)) +
  geom_violin(alpha = 0.65) +
  geom_point(size = 1.5, alpha = 0.75) +
  facet_grid(~ condition, scales = "free_x") +
  scale_fill_viridis_d(option = "plasma") +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Detected Proteins per Fiber (Proteomics)",
    x = "Donor",
    y = "Detected Proteins"
  )
rna / prot


# Scoring and integration with functional data ---------------------------------
cluster_fibers <- c(
  "c22f17", "s14_2f10", "s14_2f8", "s15_2f12", "s15_2f13",
  "s15_2f2", "s15_2f3", "s15_2f7", "s15_2f9", "s20f1",
  "s20f3", "s21f15", "s23f1", "s23f10", "s23f12", "s23f2",
  "s23f3", "s23f4", "s23f5", "s23f6", "s23f7", "s23f8",
  "s23f9", "s26f1", "s26f12", "s26f3", "s5_2f7", "s5_2f8",
  "s8_2f16", "s8_2f2", "s8_2f3", "s8_2f4"
)
rna_norm <- readRDS("../pilot-rna/objects/seurat_normalized.rds")
DefaultAssay(rna_norm) <- "RNA"
rna_matrix <- GetAssayData(rna_norm, slot = "data")
prot_matrix <- normalized_data

rna_overlap <- rna_matrix[rownames(rna_matrix) %in% overlap_genes, , drop = FALSE]
prot_overlap <- prot_matrix[rownames(prot_matrix) %in% overlap_genes, , drop = FALSE]

common_cells <- intersect(colnames(rna_overlap), colnames(prot_overlap))
rna_overlap_sub  <- rna_overlap[, common_cells, drop = FALSE]
prot_overlap_sub <- prot_overlap[, common_cells, drop = FALSE]

rna_score  <- colMeans(rna_overlap_sub,  na.rm = TRUE)
prot_score <- colMeans(prot_overlap_sub, na.rm = TRUE)

rna_score_z  <- as.numeric(scale(rna_score))
prot_score_z <- as.numeric(scale(prot_score))

combined_score <- (rna_score_z + prot_score_z) / 2
names(combined_score) <- common_cells

rna_norm@meta.data <- rna_norm@meta.data[1:length(colnames(rna_overlap)), ]
rownames(rna_norm@meta.data) <- colnames(rna_overlap)
rna_norm@meta.data$in_cluster <- ifelse(
  rownames(rna_norm@meta.data) %in% cluster_fibers,
  "cluster",
  "other"
)

meta_data_subset <- rna_norm@meta.data[common_cells, ]


meta_data_subset$Combined_Score <- combined_score
meta_data_subset$in_cluster <- factor(
  ifelse(rownames(meta_data_subset) %in% cluster_fibers, "cluster", "other"),
  levels = c("other", "cluster")
)

func_metrics <- c("srx","drx","t1","t2")
meta_data_subset[, func_metrics] <- lapply(meta_data_subset[, func_metrics], as.numeric)
meta_data_subset <- meta_data_subset[complete.cases(meta_data_subset[, func_metrics]), ]

p_srx <- wilcox.test(meta_data_subset$srx ~ meta_data_subset$in_cluster)$p.value
p_drx <- wilcox.test(meta_data_subset$drx ~ meta_data_subset$in_cluster)$p.value
p_t1  <- wilcox.test(meta_data_subset$t1  ~ meta_data_subset$in_cluster)$p.value
p_t2  <- wilcox.test(meta_data_subset$t2  ~ meta_data_subset$in_cluster)$p.value
atp  <- wilcox.test(meta_data_subset$ATP_turnover  ~ meta_data_subset$in_cluster)$p.value

plot_srx <- ggplot(meta_data_subset, aes(x = srx, y = Combined_Score, color = in_cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("cluster"="red", "other"="gray70")) +
  theme_bw() +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("p = ", signif(p_srx, 3))),
    hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 5
  ) +
  labs(
    x = "SRX",
    y = "Combined RNA+Protein Score",
    color = "Cluster",
    title = "SRX vs Combined Score"
  )

plot_drx <- ggplot(meta_data_subset, aes(x = drx, y = Combined_Score, color = in_cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("cluster"="red", "other"="gray70")) +
  theme_bw() +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("p = ", signif(p_drx, 3))),
    hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 5
  ) +
  labs(
    x = "DRX",
    y = "Combined RNA+Protein Score",
    color = "Cluster",
    title = "DRX vs Combined Score"
  )

plot_t1 <- ggplot(meta_data_subset, aes(x = t1, y = Combined_Score, color = in_cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("cluster"="red", "other"="gray70")) +
  theme_bw() +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("p = ", signif(p_t1, 3))),
    hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 5
  ) +
  labs(
    x = "T1",
    y = "Combined RNA+Protein Score",
    color = "Cluster",
    title = "T1 vs Combined Score"
  )

plot_t2 <- ggplot(meta_data_subset, aes(x = t2, y = Combined_Score, color = in_cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("cluster"="red", "other"="gray70")) +
  theme_bw() +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("p = ", signif(p_t2, 3))),
    hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 5
  ) +
  labs(
    x = "T2",
    y = "Combined RNA+Protein Score",
    color = "Cluster",
    title = "T2 vs Combined Score"
  )

combined_plot <- (plot_srx + plot_drx) / (plot_t1 + plot_t2)
combined_plot

meta_data_subset$ATP_turnover <- 220 * (
  (meta_data_subset$drx * 60) / (100 * meta_data_subset$t1) +
    (meta_data_subset$srx * 60) / (100 * meta_data_subset$t2)
)
ggplot(meta_data_subset, aes(x = ATP_turnover, y = Combined_Score, color = in_cluster)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("cluster"="red", "other"="gray70")) +
  theme_bw() +
  geom_text(
    aes(x = Inf, y = Inf, label = paste0("p = ", signif(atp, 3))),
    hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 5
  ) +
  labs(
    x = "Time (seconds)",
    y = "Combined RNA+Protein Score",
    color = "Cluster",
    title = "Theoretical ATP turnover time"
  )

# Citations
c("vroom", "here", "PhosR", "msigdbr", "readr", "gt", "dplyr", "Seurat", "ggplot2", "ggrepel", "tidyr", "stringr", "clusterProfiler", "org.Hs.eg.db", "scran", "scuttle", "tibble", "edgeR", "pheatmap", "fgsea", "viridis", "limma", "enrichplot") %>%
  map(citation) %>%
  print(style = "text")
