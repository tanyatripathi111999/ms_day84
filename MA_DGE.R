# ================================
# MA plots (same X for all; Y fixed to ±2 with decimal labels)
# ================================
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(rlang)
  library(scales)   # for label_number()
})

# ---- Set your file paths ----
file_coimbra     <- "/Users/tanyatripathi/Desktop/MS_data/FINAL/new_analysis_raw/DGE_Coimbra_MS_all_genes_edgeR.xlsx"
file_noncoimbra  <- "/Users/tanyatripathi/Desktop/MS_data/FINAL/new_analysis_raw/DGE_NonCoimbra_MS_all_genes_edgeR.xlsx"
file_healthy     <- "/Users/tanyatripathi/Desktop/MS_data/FINAL/new_analysis_raw/DGE_Healthy_all_genes_edgeR.xlsx"

# ---- thresholds ----
cpm_thresh <- 1.5
lfc_thresh <- 0.25
pval_thresh <- 0.05

# y-axis settings
ylim_range <- c(-1.25, 1.25)     # FIXED
y_break_by <- 0.5         # tick every 0.2 (change to 0.1 or 0.5 if you want)

# ---- reader (robust to column name case) ----
read_dge <- function(path) {
  tb <- read_excel(path)
  names(tb) <- tolower(names(tb))
  req <- c("hgnc_symbol","logcpm","logfc","pvalue")
  miss <- setdiff(req, names(tb))
  if (length(miss)) stop("Missing columns in ", path, ": ", paste(miss, collapse=", "))
  tb %>%
    transmute(
      hgnc_symbol = as.character(.data$hgnc_symbol),
      logCPM = as.numeric(.data$logcpm),
      logFC  = as.numeric(.data$logfc),
      PValue = as.numeric(.data$pvalue)
    ) %>%
    filter(is.finite(logCPM), is.finite(logFC), is.finite(PValue))
}

coimbra     <- read_dge(file_coimbra)
noncoimbra  <- read_dge(file_noncoimbra)
healthy     <- read_dge(file_healthy)

# ---- shared X-axis (same across all three) with a little padding ----
all_logCPM <- c(coimbra$logCPM, noncoimbra$logCPM, healthy$logCPM)
x_min <- min(all_logCPM, na.rm=TRUE)
x_max <- max(all_logCPM, na.rm=TRUE)
x_pad <- max((x_max - x_min) * 0.03, 0.25)
xlim_range <- c(floor(x_min) - x_pad, ceiling(x_max) + x_pad)

# ---- plot helper (NS under, significant over) ----
mk_ma_plot <- function(df, title, xlim_range, ylim_range) {
  df <- df %>%
    mutate(
      status = case_when(
        PValue < pval_thresh & logCPM > cpm_thresh &  logFC >  lfc_thresh ~ "Upregulated",
        PValue < pval_thresh & logCPM > cpm_thresh &  logFC < -lfc_thresh ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      status = factor(status, levels = c("Not significant","Downregulated","Upregulated"))
    )
  
  df_ns  <- df %>% filter(status == "Not significant")
  df_sig <- df %>% filter(status != "Not significant")
  
  up_n   <- sum(df$status == "Upregulated",   na.rm = TRUE)
  down_n <- sum(df$status == "Downregulated", na.rm = TRUE)
  ns_n   <- sum(df$status == "Not significant", na.rm = TRUE)
  
  legend_labels <- c(
    paste0("Upregulated (", up_n, ")"),
    paste0("Downregulated (", down_n, ")"),
    paste0("Not significant (", ns_n, ")")
  )
  
  ggplot() +
    # non-significant underneath
    geom_point(data = df_ns, aes(x = logCPM, y = logFC, color = status),
               size = 2.8, alpha = 0.9) +
    # significant on top
    geom_point(data = df_sig, aes(x = logCPM, y = logFC, color = status),
               size = 2.8, alpha = 0.95) +
    # thresholds
    geom_hline(yintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
    geom_vline(xintercept = cpm_thresh, linetype = "dotted") +
    # fixed limits (no expansion) and friendly decimal labels
    scale_x_continuous(limits = xlim_range, expand = expansion(mult = 0)) +
    scale_y_continuous(
      limits = ylim_range,
      breaks = seq(ylim_range[1], ylim_range[2], by = y_break_by),
      labels = label_number(accuracy = 0.01),   # 0.1 precision, no scientific notation
      expand = expansion(mult = 0)
    ) +
    scale_color_manual(
      values = c("Upregulated"="darkred", "Downregulated"="darkblue", "Not significant"="grey70"),
      breaks = c("Upregulated","Downregulated","Not significant"),
      labels = legend_labels
    ) +
    labs(title = title, x = "logCPM", y = "log2 Fold Change", color = NULL) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line  = element_line(color = "black"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text  = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 13),
      plot.title   = element_text(size = 16, face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

# ---- build & show ----
p_coimbra    <- mk_ma_plot(coimbra,    "MA Plot — MS Coimbra (Day 84 vs Day 0)",     xlim_range, ylim_range)
p_noncoimbra <- mk_ma_plot(noncoimbra, "MA Plot — MS Non-Coimbra (Day 84 vs Day 0)", xlim_range, ylim_range)
p_healthy    <- mk_ma_plot(healthy,    "MA Plot — Healthy (Day 84 vs Day 0)",        xlim_range, ylim_range)

graphics.off()
if (Sys.getenv("RSTUDIO") == "1") options(device = "RStudioGD")
print(p_coimbra)
print(p_noncoimbra)
print(p_healthy)

# Optional saves:
# ggsave("MA_MS_Coimbra.png",     p_coimbra,    width=8, height=6, dpi=300)
# ggsave("MA_MS_NonCoimbra.png",  p_noncoimbra, width=8, height=6, dpi=300)
# ggsave("MA_Healthy.png",        p_healthy,    width=8, height=6, dpi=300)




















# =========================
# edgeR unpaired DGE on RAW COUNTS
# d0_MS vs d0_Ind  AND  d84_MS vs d84_Ind
# with stricter filtering (expressed in BOTH groups) + identical MA scaling
# =========================

suppressPackageStartupMessages({
  library(edgeR)
  library(readxl)
  library(dplyr)
  library(writexl)
  library(ggplot2)
  library(qvalue)
  library(tibble)
})

# ==== Paths ====
file_path <- "/Users/tanyatripathi/Desktop/MS_data/FINAL/new_analysis_raw/mixed_raw.xlsx"
outdir    <- "/Users/tanyatripathi/Desktop/MS_data/FINAL/new_analysis_raw"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==== Load RAW COUNT data (NOT CPM/TPM) ====
df <- read_excel(file_path)
stopifnot("hgnc_symbol" %in% names(df))
df <- df %>% filter(!is.na(hgnc_symbol))
rownames(df) <- df$hgnc_symbol

# ---- Column groups ----
ms_d0_cols       <- grep("^d0_MS",   names(df), value = TRUE)
healthy_d0_cols  <- grep("^d0_Ind",  names(df), value = TRUE)
ms_d84_cols      <- grep("^d84_MS",  names(df), value = TRUE)
healthy_d84_cols <- grep("^d84_Ind", names(df), value = TRUE)

# ========================= Helpers =========================

# Stricter filter: keep genes with CPM >= cpm_cut in AT LEAST min_per_grp samples
# in EACH group (prevents one-group-only artifacts)
filter_both_groups <- function(counts, group, cpm_cut = 1, min_per_grp = NULL) {
  y <- DGEList(counts = counts)
  cpm_mat <- cpm(y)
  glev <- levels(group)
  if (is.null(min_per_grp)) {
    # sensible default: at least 2 per group, or half of the smaller group
    min_per_grp <- max(2, floor(min(table(group)) / 2))
  }
  keep <- rowSums(cpm_mat[, group == glev[1], drop = FALSE] >= cpm_cut) >= min_per_grp &
    rowSums(cpm_mat[, group == glev[2], drop = FALSE] >= cpm_cut) >= min_per_grp
  keep
}

# Run one unpaired contrast: group1 (MS) vs group2 (Healthy)
run_unpaired_dge <- function(group1_cols, group2_cols, group1_label, group2_label,
                             both_groups_filter = TRUE, cpm_cut = 1, min_per_grp = NULL) {
  cols <- c(group1_cols, group2_cols)
  counts <- as.matrix(sapply(df[, cols], as.numeric))
  rownames(counts) <- df$hgnc_symbol
  
  # guards: RAW integer counts only
  if (any(is.na(counts))) stop("Counts contain NA.")
  if (any(counts < 0)) stop("Counts contain negatives.")
  if (any(abs(counts - round(counts)) > 1e-6)) stop("Input looks like CPM/TPM. Provide RAW counts.")
  
  group <- factor(c(rep(group1_label, length(group1_cols)),
                    rep(group2_label, length(group2_cols))))
  group <- relevel(group, ref = group2_label)  # Healthy as reference
  
  # Optional diagnostic: how many "one-group-only" genes before filtering?
  y_diag <- DGEList(counts = counts)
  cpm_mat <- cpm(y_diag)
  glev <- levels(group)
  n1 <- sum(group == glev[1]); n2 <- sum(group == glev[2])
  only_g1 <- rowSums(cpm_mat[, group==glev[1], drop=FALSE] >= 1) >= max(2, floor(n1/2)) &
    rowSums(cpm_mat[, group==glev[2], drop=FALSE] >= 1) == 0
  only_g2 <- rowSums(cpm_mat[, group==glev[2], drop=FALSE] >= 1) >= max(2, floor(n2/2)) &
    rowSums(cpm_mat[, group==glev[1], drop=FALSE] >= 1) == 0
  
  # Filtering
  if (both_groups_filter) {
    keep <- filter_both_groups(counts, group, cpm_cut = cpm_cut, min_per_grp = min_per_grp)
  } else {
    y_tmp <- DGEList(counts = counts, group = group)
    keep <- filterByExpr(y_tmp, group = group)
  }
  
  # edgeR pipeline
  y <- DGEList(counts = counts, group = group)
  total_genes <- nrow(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  kept_genes <- nrow(y)
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ group)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = 2)  # group1 vs reference
  
  res <- topTags(qlf, n = Inf)$table %>%
    rownames_to_column("hgnc_symbol") %>%
    mutate(
      FDR_BH  = p.adjust(PValue, method = "BH"),
      FDR_BY  = p.adjust(PValue, method = "BY"),
      q_value = qvalue(PValue)$qvalues
    )
  
  # summary for console
  message(sprintf(
    "[%s] total=%d, kept=%d, dropped=%d, only_%s=%d, only_%s=%d",
    paste0(group1_label," vs ",group2_label),
    total_genes, kept_genes, total_genes - kept_genes,
    glev[1], sum(only_g1), glev[2], sum(only_g2)
  ))
  
  res
}

# MA plot (FDR-based coloring); pass shared x/y ranges for identical scaling
plot_MA <- function(data, title, xlim_range, ylim_range) {
  lfc_thresh <- 1
  cpm_thresh <- 1.5
  fdr_thresh <- 0.05
  
  data <- data %>% filter(is.finite(logCPM)) %>%
    mutate(
      status = case_when(
        FDR_BH < fdr_thresh & logFC >  lfc_thresh & logCPM > cpm_thresh ~ "Upregulated",
        FDR_BH < fdr_thresh & logFC < -lfc_thresh & logCPM > cpm_thresh ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      status = factor(status, levels = c("Not significant","Downregulated","Upregulated"))
    )
  
  up_count <- sum(data$status == "Upregulated")
  down_count <- sum(data$status == "Downregulated")
  legend_labels <- c(
    paste0("Upregulated (", up_count, ")"),
    paste0("Downregulated (", down_count, ")"),
    "Not significant"
  )
  
  ggplot(data, aes(x = logCPM, y = logFC, color = status)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("Not significant"="grey70","Downregulated"="darkblue","Upregulated"="darkred"),
                       breaks = c("Upregulated","Downregulated","Not significant"),
                       labels = legend_labels) +
    geom_hline(yintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
    geom_vline(xintercept = c(cpm_thresh), linetype = "dotted") +
    coord_cartesian(xlim = xlim_range, ylim = ylim_range) +
    labs(title = title, x = "logCPM", y = "log2 Fold Change", color = "Gene Status") +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold")
    )
}

# ========================= Run both contrasts =========================
res_d0  <- run_unpaired_dge(ms_d0_cols,  healthy_d0_cols,  "d0_MS",  "d0_Ind",
                            both_groups_filter = TRUE, cpm_cut = 1, min_per_grp = NULL)
res_d84 <- run_unpaired_dge(ms_d84_cols, healthy_d84_cols, "d84_MS", "d84_Ind",
                            both_groups_filter = TRUE, cpm_cut = 1, min_per_grp = NULL)

# ---- Shared axes (identical scaling) ----
x_min <- floor(min(c(res_d0$logCPM, res_d84$logCPM), na.rm = TRUE))
x_max <- ceiling(max(c(res_d0$logCPM, res_d84$logCPM), na.rm = TRUE))
ycap  <- as.numeric(quantile(abs(c(res_d0$logFC, res_d84$logFC)), 0.995, na.rm = TRUE))
y_lim <- c(-ycap, ycap)

# ---- MA plots (display + save) ----
p_d0  <- plot_MA(res_d0,  "MA Plot: d0_MS vs d0_Ind",   c(x_min, x_max), y_lim)
p_d84 <- plot_MA(res_d84, "MA Plot: d84_MS vs d84_Ind", c(x_min, x_max), y_lim)
print(p_d0); print(p_d84)
ggsave(file.path(outdir, "MA_d0_MS_vs_d0_Ind.png"),  p_d0,  width = 8, height = 6, dpi = 300)
ggsave(file.path(outdir, "MA_d84_MS_vs_d84_Ind.png"), p_d84, width = 8, height = 6, dpi = 300)

# ---- Filtered DEGs (adjust if you prefer P-value) ----
deg_lfc <- 1.0  # |log2FC| > 1 (~2x)
deg_cpm <- 1.5

res_d0_FDR  <- res_d0  %>% filter(FDR_BH < 0.05, abs(logFC) > deg_lfc, logCPM > deg_cpm)
res_d84_FDR <- res_d84 %>% filter(FDR_BH < 0.05, abs(logFC) > deg_lfc, logCPM > deg_cpm)

# Optional: P-value based sets too
res_d0_P  <- res_d0  %>% filter(PValue < 0.05, abs(logFC) > deg_lfc, logCPM > deg_cpm)
res_d84_P <- res_d84 %>% filter(PValue < 0.05, abs(logFC) > deg_lfc, logCPM > deg_cpm)

# ---- Save results ----
final_all_d0  <- df %>% select(hgnc_symbol, all_of(c(ms_d0_cols,  healthy_d0_cols)))  %>% left_join(res_d0,  by = "hgnc_symbol")
final_all_d84 <- df %>% select(hgnc_symbol, all_of(c(ms_d84_cols, healthy_d84_cols))) %>% left_join(res_d84, by = "hgnc_symbol")

write_xlsx(final_all_d0,  file.path(outdir, "DGE_d0_MS_vs_d0_Ind_all_genes_normal.xlsx"))
write_xlsx(final_all_d84, file.path(outdir, "DGE_d84_MS_vs_d84_Ind_all_genes_normal.xlsx"))
write_xlsx(res_d0_FDR,    file.path(outdir, "DGE_d0_MS_vs_d0_Ind_filtered_FDR.xlsx"))
write_xlsx(res_d84_FDR,   file.path(outdir, "DGE_d84_MS_vs_d84_Ind_filtered_FDR.xlsx"))
write_xlsx(res_d0_P,      file.path(outdir, "DGE_d0_MS_vs_d0_Ind_filtered_P.xlsx"))
write_xlsx(res_d84_P,     file.path(outdir, "DGE_d84_MS_vs_d84_Ind_filtered_P.xlsx"))

cat("✅ Finished: d0 & d84 contrasts on RAW counts with both-groups filtering. Plots saved to:", outdir, "\n")


