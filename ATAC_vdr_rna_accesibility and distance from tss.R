# =========================
# FROM SCRATCH:
# 1) ATAC peaks common with VDR peaks (overlap or within 500 bp)
# 2) RNA genes whose TSS is within ±500 kb of those common peaks
# 3) Distance from each gene TSS to nearest common peak (bp)
# =========================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(GenomicRanges)
  library(openxlsx)
})

# --------- PATHS ----------
atac_file <- "/Users/tanyatripathi/Dropbox/Lab stuff Carlberg group/MS manuscript 1/Final/Table S3.xlsx"
atac_sheet <- 1

vdr_file <- "/Users/tanyatripathi/Desktop/VDR_MS overlap_D84/priman_vdr_binding_sites.xlsx"
vdr_sheet <- 1

rna_file <- "/Users/tanyatripathi/Dropbox/Lab stuff Carlberg group/MS manuscript 1/Final/Table S2.xlsx"
rna_sheet <- 1

out_file <- "/Users/tanyatripathi/Desktop/VDR_MS overlap_D84/ATAC_VDR_RNA_annotated_NEW.xlsx"

# --------- PARAMETERS ----------
maxgap_atac_vdr_bp <- 500      # ATAC-VDR "common" if overlap OR within 500 bp
tss_window_bp      <- 500000   # RNA gene TSS within ±500 kb of any common peak

# --------- HELPERS ----------
standardize_chr <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  has_chr <- grepl("^chr", x, ignore.case = TRUE)
  x[!has_chr] <- paste0("chr", x[!has_chr])
  x <- sub("^CHR", "chr", x)
  x <- sub("^Chr", "chr", x)
  x
}

to_int <- function(x) suppressWarnings(as.integer(as.character(x)))

must_have_cols <- function(df, cols, df_name) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) stop(df_name, " missing columns: ", paste(miss, collapse = ", "))
}

sanity_intervals <- function(df, name) {
  cat("\n--- Sanity:", name, "---\n")
  cat("Rows:", nrow(df), "\n")
  cat("NA chr/start/end:", sum(is.na(df$Chromosome)), "/", sum(is.na(df$Start)), "/", sum(is.na(df$End)), "\n")
  bad <- which(df$Start > df$End)
  cat("Start > End:", length(bad), "\n")
  if (length(bad) > 0) cat("Example bad rows:", paste(head(bad, 5), collapse = ", "), "\n")
  w <- df$End - df$Start + 1
  cat("Width summary:\n"); print(summary(w))
}

swap_if_needed <- function(df) {
  idx <- which(df$Start > df$End)
  if (length(idx) > 0) {
    tmp <- df$Start[idx]
    df$Start[idx] <- df$End[idx]
    df$End[idx] <- tmp
  }
  df
}

# =========================
# 1) READ DATA
# =========================
atac <- read_excel(atac_file, sheet = atac_sheet) %>% as.data.frame()
vdr  <- read_excel(vdr_file,  sheet = vdr_sheet)  %>% as.data.frame()
rna  <- read_excel(rna_file,  sheet = rna_sheet)  %>% as.data.frame()

# Required columns
must_have_cols(atac, c("Chromosome", "Start", "End"), "ATAC")
must_have_cols(vdr,  c("Chromosome", "Start", "End"), "VDR")
must_have_cols(rna,  c("Chromosome", "Start", "End", "Strand"), "RNA")

# Standardize chr + convert coords
atac$Chromosome <- standardize_chr(atac$Chromosome)
vdr$Chromosome  <- standardize_chr(vdr$Chromosome)
rna$Chromosome  <- standardize_chr(rna$Chromosome)

atac$Start <- to_int(atac$Start); atac$End <- to_int(atac$End)
vdr$Start  <- to_int(vdr$Start);  vdr$End  <- to_int(vdr$End)
rna$Start  <- to_int(rna$Start);  rna$End  <- to_int(rna$End)

# Drop missing coords (strict)
atac <- atac %>% filter(!is.na(Chromosome), !is.na(Start), !is.na(End))
vdr  <- vdr  %>% filter(!is.na(Chromosome), !is.na(Start), !is.na(End))
rna  <- rna  %>% filter(!is.na(Chromosome), !is.na(Start), !is.na(End), !is.na(Strand))

# Fix swapped coords
atac <- swap_if_needed(atac)
vdr  <- swap_if_needed(vdr)
rna  <- swap_if_needed(rna)

# Sanity checks
sanity_intervals(atac, "ATAC")
sanity_intervals(vdr,  "VDR")
sanity_intervals(rna,  "RNA (gene body)")

# Chromosome overlap sanity
cat("\n--- Chromosome sanity ---\n")
cat("ATAC chr count:", length(unique(atac$Chromosome)), "\n")
cat("VDR  chr count:", length(unique(vdr$Chromosome)), "\n")
cat("RNA  chr count:", length(unique(rna$Chromosome)), "\n")
cat("Chr in ATAC not in VDR:", paste(setdiff(unique(atac$Chromosome), unique(vdr$Chromosome)), collapse=", "), "\n")
cat("Chr in VDR not in ATAC:", paste(setdiff(unique(vdr$Chromosome), unique(atac$Chromosome)), collapse=", "), "\n")

# =========================
# 2) GRanges OBJECTS
# =========================
gr_atac <- GRanges(
  seqnames = atac$Chromosome,
  ranges   = IRanges(start = atac$Start, end = atac$End)
)

gr_vdr <- GRanges(
  seqnames = vdr$Chromosome,
  ranges   = IRanges(start = vdr$Start, end = vdr$End)
)

# =========================
# 3) STEP A: ATAC common with VDR (<=500 bp)
# =========================
hits_av <- findOverlaps(gr_atac, gr_vdr, maxgap = maxgap_atac_vdr_bp, type = "any", ignore.strand = TRUE)

cat("\n--- ATAC vs VDR overlap summary ---\n")
cat("ATAC peaks:", length(gr_atac), "\n")
cat("VDR peaks :", length(gr_vdr), "\n")
cat("Hit pairs :", length(hits_av), "\n")
cat("Unique ATAC peaks hit:", length(unique(queryHits(hits_av))), "\n")
cat("Unique VDR peaks hit :", length(unique(subjectHits(hits_av))), "\n")

atac$Common_with_VDR_within500bp <- "No"
atac$Common_with_VDR_within500bp[unique(queryHits(hits_av))] <- "Yes"

# Define "common accessible sites" as the ATAC peaks that match VDR (standard choice)
gr_common <- gr_atac[unique(queryHits(hits_av))]
cat("Common ATAC peaks used as 'accessible sites':", length(gr_common), "\n")

# =========================
# 4) STEP B: RNA TSS within ±500 kb of any common accessible site
# =========================
rna$Strand <- as.character(rna$Strand)
rna$TSS <- ifelse(rna$Strand %in% c("+", "plus", "1"), rna$Start, rna$End)
rna$TSS <- to_int(rna$TSS)

cat("\n--- TSS sanity ---\n")
cat("NA TSS:", sum(is.na(rna$TSS)), "\n")
cat("TSS summary:\n"); print(summary(rna$TSS))

gr_tss <- GRanges(
  seqnames = rna$Chromosome,
  ranges   = IRanges(start = rna$TSS, end = rna$TSS)
)

# Extend common sites by ±500kb
gr_common_ext <- gr_common
start(gr_common_ext) <- pmax(1L, start(gr_common_ext) - tss_window_bp)
end(gr_common_ext)   <- end(gr_common_ext) + tss_window_bp

hits_tss <- findOverlaps(gr_tss, gr_common_ext, type = "any", ignore.strand = TRUE)

cat("\n--- RNA TSS proximity summary ---\n")
cat("RNA genes:", length(gr_tss), "\n")
cat("Hit pairs:", length(hits_tss), "\n")
cat("Unique genes within ±500kb:", length(unique(queryHits(hits_tss))), "\n")

rna$TSS_within_500kb_of_common_ATAC_VDR_site <- "No"
rna$TSS_within_500kb_of_common_ATAC_VDR_site[unique(queryHits(hits_tss))] <- "Yes"

# =========================
# 5) STEP C: Distance (bp) from each gene TSS to nearest common accessible site
# =========================
rna$Distance_bp_to_nearest_common_ATAC_VDR_site <- NA_integer_

if (length(gr_common) == 0) {
  cat("\nNo common ATAC–VDR sites found -> distance column will remain NA.\n")
} else {
  nn <- distanceToNearest(gr_tss, gr_common, ignore.strand = TRUE)
  rna$Distance_bp_to_nearest_common_ATAC_VDR_site[queryHits(nn)] <- mcols(nn)$distance
  
  cat("\n--- Distance sanity ---\n")
  print(summary(rna$Distance_bp_to_nearest_common_ATAC_VDR_site))
  cat("Genes with distance computed:", sum(!is.na(rna$Distance_bp_to_nearest_common_ATAC_VDR_site)), "\n")
  
  # Strong sanity check: all "Yes" genes should have distance <= 500kb
  yes_idx <- which(rna$TSS_within_500kb_of_common_ATAC_VDR_site == "Yes")
  if (length(yes_idx) > 0) {
    max_yes <- max(rna$Distance_bp_to_nearest_common_ATAC_VDR_site[yes_idx], na.rm = TRUE)
    cat("Max distance among 'Yes' genes (should be <= 500000):", max_yes, "\n")
  }
}

# =========================
# 6) WRITE OUTPUT EXCEL
# =========================
wb <- createWorkbook()
addWorksheet(wb, "ATAC_annotated")
addWorksheet(wb, "RNA_annotated")
addWorksheet(wb, "VDR_original")

writeData(wb, "ATAC_annotated", atac)
writeData(wb, "RNA_annotated",  rna)
writeData(wb, "VDR_original",   vdr)

saveWorkbook(wb, out_file, overwrite = TRUE)
cat("\n✅ Done. Output written to:\n", out_file, "\n")

