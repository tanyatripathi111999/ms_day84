# install.packages(c("readxl","dplyr","ggplot2","stringr"))  # if needed
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

# ==== 1) Read your Excel ====
# Change these to your actual file & sheet name
xlsx_path <- "/Users/tanyatripathi/Downloads/125_levels.xlsx"
sheet_name <- 1  # or "Sheet1"

raw <- read_excel(xlsx_path, sheet = sheet_name, col_types = c("text","text"))


# Convert comma decimals to numeric
df <- raw %>%
  rename(OHD_raw = `25OHD (ng/ml)`) %>%
  mutate(OHD = as.numeric(gsub(",", ".", OHD_raw, fixed = TRUE)))

# ==== 2) Assign groups & IDs ====
df <- df %>%
  mutate(
    isMS   = str_starts(Code, "MS"),
    isInd  = str_starts(Code, "Ind"),
    ID_MS  = suppressWarnings(as.integer(str_match(Code, "^MS(\\d+)_")[,2])),
    ID_Ind = suppressWarnings(as.integer(str_match(Code, "^Ind(\\d+)_")[,2])),
    Time   = case_when(
      str_detect(Code, "_1$") ~ "d0",
      str_detect(Code, "_2$") ~ "d84",
      TRUE ~ NA_character_
    ),
    Cohort = case_when(
      isMS & !is.na(ID_MS) & ID_MS <= 11 ~ "Coimbra",
      isMS & !is.na(ID_MS) & ID_MS >  11 ~ "NonCoimbra",
      isInd                              ~ "Healthy",
      TRUE ~ NA_character_
    ),
    Group = paste0(Cohort, "_", Time)
  )

# Order for x-axis
lvl <- c("Coimbra_d0","Coimbra_d84","NonCoimbra_d0","NonCoimbra_d84","Healthy_d0","Healthy_d84")
df <- df %>% filter(!is.na(Group), !is.na(OHD))
df$Group <- factor(df$Group, levels = lvl)

# ==== 3) Legend stats (mean & n) ====
avg_df <- df %>%
  group_by(Group) %>%
  summarise(mean_OHD = mean(OHD, na.rm = TRUE), count = n(), .groups = "drop") %>%
  mutate(label = paste0(Group, " (Avg: ", round(mean_OHD, 2), ", n=", count, ")"))

df <- df %>% left_join(avg_df %>% select(Group, label), by = "Group")

# ==== 4) Compute paired p-values (Wilcoxon) for each cohort ====
# We pair by subject ID within each cohort.
paired_pvals <- bind_rows(
  # Coimbra (MS01–MS11)
  df %>% filter(Cohort == "Coimbra") %>%
    mutate(ID = ID_MS) %>%
    select(Cohort, ID, Time, OHD) %>%
    pivot_wider(names_from = Time, values_from = OHD) %>%
    drop_na(d0, d84) %>%
    summarise(
      p = wilcox.test(d0, d84, paired = TRUE)$p.value,
      group1 = "Coimbra_d0", group2 = "Coimbra_d84"
    ),
  # Non-Coimbra (MS12–MS18)
  df %>% filter(Cohort == "NonCoimbra") %>%
    mutate(ID = ID_MS) %>%
    select(Cohort, ID, Time, OHD) %>%
    pivot_wider(names_from = Time, values_from = OHD) %>%
    drop_na(d0, d84) %>%
    summarise(
      p = wilcox.test(d0, d84, paired = TRUE)$p.value,
      group1 = "NonCoimbra_d0", group2 = "NonCoimbra_d84"
    ),
  # Healthy (Ind*)
  df %>% filter(Cohort == "Healthy") %>%
    mutate(ID = ID_Ind) %>%
    select(Cohort, ID, Time, OHD) %>%
    pivot_wider(names_from = Time, values_from = OHD) %>%
    drop_na(d0, d84) %>%
    summarise(
      p = wilcox.test(d0, d84, paired = TRUE)$p.value,
      group1 = "Healthy_d0", group2 = "Healthy_d84"
    )
) %>%
  # y positions a bit above the max of each pair
  rowwise() %>%
  mutate(
    y.position = max(df$OHD[df$Group %in% c(group1, group2)], na.rm = TRUE) + 5,
    p.label = paste0("p = ", signif(p, 3))
  ) %>%
  ungroup()





# Get unique labels
labels_unique <- unique(avg_df$label)

# Assign colors to each label
bar_cols  <- c("skyblue","blue","orange","red","gray70","gray30")[seq_along(labels_unique)]
point_cols<- c("steelblue4","blue4","orange4","red4","gray50","black")[seq_along(labels_unique)]
names(bar_cols)   <- labels_unique
names(point_cols) <- labels_unique

# Plot
p <- ggplot(df, aes(x = Group, y = OHD, fill = label)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black") +
  geom_jitter(aes(color = label), width = 0.15, size = 2) +
  scale_fill_manual(values = bar_cols) +
  scale_color_manual(values = point_cols) +
  labs(
    title = "25OHD Levels by Group (with paired Wilcoxon p-values)",
    x = NULL, y = "25OHD (ng/ml)",
    fill = "Groups (Avg, n)", color = "Groups (Avg, n)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  stat_pvalue_manual(
    paired_pvals,
    label = "p.label",
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  )

print(p)

p <- ggplot(df, aes(x = Group, y = OHD, fill = label)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black") +
  geom_jitter(aes(color = label), width = 0.15, size = 2.5) +
  scale_fill_manual(values = bar_cols) +
  scale_color_manual(values = point_cols) +
  labs(
    title = "25OHD Levels by Group (with paired Wilcoxon p-values)",
    x = NULL, y = "25OHD (ng/ml)",
    fill = "Groups (Avg, n)", color = "Groups (Avg, n)"
  ) +
  theme_classic(base_size = 18) +  # increase base font
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", size = 16),
    axis.text.y = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.line = element_line(linewidth = 1.2)
  ) +
  stat_pvalue_manual(
    paired_pvals,
    label = "p.label",
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 6,              # 🔥 makes the p-value text bigger
    fontface = "bold"      # 🔥 makes it bold too
  )

print(p)
