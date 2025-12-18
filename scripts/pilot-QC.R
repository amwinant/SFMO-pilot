#RNA
meta_all <- seurat_all@meta.data
meta_all$nFeature <- seurat_all$nFeature_RNA
meta_all$nCount   <- seurat_all$nCount_RNA
feat_med <- median(meta_all$nFeature)
feat_mad <- mad(meta_all$nFeature)
count_med <- median(meta_all$nCount)
count_mad <- mad(meta_all$nCount)
meta_all$keep <- with(meta_all,
                      nFeature > feat_med - 3 * feat_mad &
                        nFeature < feat_med + 3 * feat_mad &
                        nCount   > count_med - 3 * count_mad &
                        nCount   < count_med + 3 * count_mad
)
ggplot(meta_all, aes(x = nCount, y = nFeature)) +
  geom_point(data = subset(meta_all, !keep),color = "gray",alpha = 0.6) +
  geom_point(data = subset(meta_all, keep),color = "black",alpha = 0.8) +
  geom_smooth(method = "lm",se = FALSE,color = "lightgray",linewidth = 1) +
  geom_vline(xintercept = c(count_med - 3 * count_mad,count_med + 3 * count_mad),
    linetype = "dashed",
    color = "red"
    ) +
  geom_hline(
    yintercept = c(feat_med - 3 * feat_mad,feat_med + 3 * feat_mad),
    linetype = "dashed",
    color = "red") +
  theme_bw() +
  labs(x = "Number of reads",y = "Number of unique genes", title = "RNA-seq QC filtering")


#Protein 
lib_size_all <- colSums(transformed_data, na.rm = TRUE)
n_proteins_all <- colSums(!is.na(transformed_data))

n_total_proteins <- nrow(transformed_data)
threshold_proteins <- 0.5 * n_total_proteins
qc_df <- data.frame(
  fiber_id = colnames(transformed_data),
  library_size = lib_size_all,
  n_proteins = n_proteins_all,
  kept = colnames(transformed_data) %in% colnames(normalized_data)
  )

ggplot(qc_df, aes(x = library_size, y = n_proteins)) +
  geom_point(aes(color = kept), alpha = 0.7) +
  scale_color_manual(values = c("TRUE"="black", "FALSE"="lightgray")) +
  geom_smooth(method = "lm", se = FALSE, color = "gray") +
  geom_hline(yintercept = threshold_proteins, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = "Total intensity per fiber",
    y = "Number of proteins",
    title = "Proteomics QC filtering"
    )
      
      
      