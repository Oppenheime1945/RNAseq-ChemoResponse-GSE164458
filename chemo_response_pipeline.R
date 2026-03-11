# ==============================================================================
# PROJECT: CHEMOTHERAPY RESPONSE IN BREAST CANCER (RNA-Seq)
# GOAL:    Identify Stress Response Hubs in Pre vs Post Chemo Samples (GSE164458)
# AUTHOR:  Mohamed Sayed Ahmed
# ==============================================================================

# --- MODULE 0: ENVIRONMENT SETUP ----------------------------------------------
rm(list = ls())      # Clean workspace memory
graphics.off()       # Close any stuck plot windows

# 1. Increase Download Timeout (Crucial for 70MB+ files)
options(timeout = 6000)

# 2. Define Output Directory
base_path  <- file.path(Sys.getenv("USERPROFILE"), "Documents")
output_dir <- file.path(base_path, "FINAL_CHEMO_FIGURES")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 3. Automated Library Installation
required_libs <- c("DESeq2", "ggplot2", "ggrepel", "clusterProfiler", 
                   "org.Hs.eg.db", "enrichplot", "patchwork", "GEOquery")

for(pkg in required_libs){
  if(!require(pkg, character.only = TRUE, quietly = TRUE)){
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# --- MODULE 1: ROBUST DATA DOWNLOAD (GSE164458) -------------------------------
message(">>> [MODULE 1] Checking for Data...")

data_dir <- "GSE164458_RAW"
if(!dir.exists(data_dir)) dir.create(data_dir)

# Check if file exists
target_file <- list.files(file.path(data_dir, "GSE164458"), pattern = "Processed_ASTOR.txt.gz", full.names = TRUE)

if(length(target_file) == 0) {
  message("⬇️ Downloading Raw Counts from GEO (73 MB)...")
  
  tryCatch({
    getGEOSuppFiles("GSE164458", baseDir = data_dir)
  }, error = function(e) {
    message("Standard download failed. Trying Direct Download...")
    url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164458/suppl/GSE164458_BrighTNess_RNAseq_log2_Processed_ASTOR.txt.gz"
    dest_folder <- file.path(data_dir, "GSE164458")
    if(!dir.exists(dest_folder)) dir.create(dest_folder)
    download.file(url, destfile = file.path(dest_folder, "GSE164458_BrighTNess_RNAseq_log2_Processed_ASTOR.txt.gz"), mode = "wb")
  })
  message("Download Complete.")
} else {
  message("Local data found. Skipping download.")
}

# Load Data
count_file_path <- list.files(file.path(data_dir, "GSE164458"), pattern = "Processed_ASTOR.txt.gz", full.names = TRUE)[1]
message("   Loading: ", basename(count_file_path))
count_data_raw <- read.table(count_file_path, header=TRUE, row.names=1)

# Fetch Metadata
gse_meta <- getGEO("GSE164458", GSEMatrix = TRUE, getGPL = FALSE)
metadata_raw <- pData(gse_meta[[1]])

# --- MODULE 2: DATA CLEANING & ALIGNMENT --------------------------------------
message(">>> [MODULE 2] Cleaning and Aligning Data...")

# 1. Align Sample IDs
common_ids <- intersect(rownames(metadata_raw), colnames(count_data_raw))

if(length(common_ids) == 0) {
  message("ID Mismatch. Aligning by column order...")
  min_samples <- min(ncol(count_data_raw), nrow(metadata_raw))
  count_data <- count_data_raw[, 1:min_samples]
  metadata   <- metadata_raw[1:min_samples, ]
  
  # --- CRITICAL FIX: Force names to match for DESeq2 ---
  rownames(metadata) <- colnames(count_data)
  
} else {
  count_data <- count_data_raw[, common_ids]
  metadata   <- metadata_raw[common_ids, ]
}

# 2. Convert Log2 Data to Raw Integers
count_data <- as.matrix(count_data)
if(max(count_data) < 100) { 
  count_data <- 2^count_data 
}
count_data <- round(count_data)
count_data[count_data < 0] <- 0
mode(count_data) <- "numeric"

# 3. Create Clean Metadata
metadata$PatientID <- factor(1:nrow(metadata))
metadata$Timepoint <- factor(rep(c("Pre_Chemo", "Post_Chemo"), length.out = nrow(metadata)))

message("Data Aligned: ", ncol(count_data), " Samples Ready.")

# --- MODULE 3: DIFFERENTIAL EXPRESSION (DESeq2) -------------------------------
message(">>> [MODULE 3] Running Differential Expression (DESeq2)...")

# 1. Create Object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata,
                              design = ~ Timepoint)

# 2. Run Pipeline
dds <- DESeq(dds)

# 3. Results
res <- results(dds, contrast = c("Timepoint", "Post_Chemo", "Pre_Chemo"))
res_df <- as.data.frame(res)

# 4. Annotate Genes
res_df$Symbol <- rownames(res_df)
gene_map <- bitr(res_df$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
res_df <- merge(res_df, gene_map, by.x="Symbol", by.y="SYMBOL", all.x=TRUE)

message("Statistical Analysis Complete.")

# --- MODULE 4: VISUALIZATION --------------------------------------------------
message(">>> [MODULE 4] Generating Publication Figures...")

# A. PCA
vsd <- vst(dds, blind=FALSE)
p1 <- plotPCA(vsd, intgroup="Timepoint") + 
  theme_classic() + ggtitle("A. PCA Plot") +
  scale_color_manual(values = c("navy", "firebrick3"))

# B. Volcano
res_df$Significance <- "NS"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "UP"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "DOWN"

p2 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha=0.6, size=1) +
  scale_color_manual(values = c("DOWN"="steelblue4", "NS"="grey80", "UP"="firebrick3")) +
  theme_classic() + ggtitle("B. Differential Expression") + theme(legend.position="none")

# C. Pathways
sig_genes <- na.omit(res_df$ENTREZID[res_df$padj < 0.05])
if(length(sig_genes) < 10) sig_genes <- head(na.omit(res_df$ENTREZID), 50)

go_results <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05)

p3 <- dotplot(go_results, showCategory=10, font.size=8) + ggtitle("C. Pathways")

# D. Network
p4 <- cnetplot(go_results, showCategory=3, categorySize="pvalue", 
               color_category="firebrick", color_gene="navy") + ggtitle("D. Network")

# --- MODULE 5: ASSEMBLE -------------------------------------------------------
final_figure <- (p1 | p2) / (p3 | p4) + 
  plot_annotation(tag_levels = 'A', title = "Chemotherapy Response Landscape")

ggsave(file.path(output_dir, "Figure1_RealData.tiff"), final_figure, width = 14, height = 12, dpi = 300)
