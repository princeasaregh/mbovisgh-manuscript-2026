getwd()
setwd("D:/Dropbox/bovineWGSA/BTBGhAnalysis/pairsnp_&_iqtree_&_itol")



# ============================================================
# Pairwise SNP distance and geographic distance analyses
# M. bovis transmission / spatial-genomic analysis
# Outputs: PNG + SVG figures and CSV statistics
# ============================================================

packages <- c(
  "tidyverse",
  "readxl",
  "ggplot2",
  "rstatix",
  "broom",
  "patchwork"
)

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)
library(readxl)
library(ggplot2)
library(rstatix)
library(broom)
library(patchwork)

# -----------------------------
# Load data
# -----------------------------

df <- read_excel(
  "Ghana_Mbovis_pairwise_geographic_SNP_distance_merged_REVISED.xlsx",
  sheet = "Merged_pairwise_data"
)

# -----------------------------
# Clean variables
# -----------------------------

put_mixed_last <- function(x) {
  x <- as.character(x)
  x[x == ""] <- NA
  
  observed <- unique(x[!is.na(x)])
  
  if ("mixed" %in% observed) {
    factor(x, levels = c(sort(observed[observed != "mixed"]), "mixed"))
  } else {
    factor(x, levels = sort(observed))
  }
}

df <- df %>%
  rename(
    snp_dist = `SNP-dist`,
    geo_dist = Geographic_distance_km
  ) %>%
  mutate(
    snp_dist = as.numeric(snp_dist),
    geo_dist = as.numeric(geo_dist),
    
    SB_plot = case_when(
      is.na(SB) | SB == "" ~ NA_character_,
      SB == "SB0944" ~ "SB0944",
      SB == "SB0300" ~ "SB0300",
      TRUE ~ "mixed"
    ),
    
    SB_plot = factor(SB_plot, levels = c("SB0944", "SB0300", "mixed")),
    Lineage = put_mixed_last(Lineage),
    Region = put_mixed_last(Region),
    Cattle_Breed = put_mixed_last(Cattle_Breed),
    District = put_mixed_last(District)
  )

# -----------------------------
# Output folders
# -----------------------------

dir.create("figures_png", showWarnings = FALSE)
dir.create("figures_svg", showWarnings = FALSE)
dir.create("statistics", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

dir.create("figures_png/geo_vs_snp", showWarnings = FALSE, recursive = TRUE)
dir.create("figures_svg/geo_vs_snp", showWarnings = FALSE, recursive = TRUE)

dir.create("figures_png/boxplots", showWarnings = FALSE, recursive = TRUE)
dir.create("figures_svg/boxplots", showWarnings = FALSE, recursive = TRUE)

dir.create("figures_png/density_plots", showWarnings = FALSE, recursive = TRUE)
dir.create("figures_svg/density_plots", showWarnings = FALSE, recursive = TRUE)

dir.create("figures_png/low_snp", showWarnings = FALSE, recursive = TRUE)
dir.create("figures_svg/low_snp", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# QC summary
# -----------------------------

required_cols <- c(
  "snp_dist", "geo_dist", "Lineage", "SB_plot",
  "Cattle_Breed", "Region", "District",
  "Sample1", "Sample2"
)

missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

qc_summary <- data.frame(
  total_pairwise_rows = nrow(df),
  missing_snp_dist = sum(is.na(df$snp_dist)),
  missing_geo_dist = sum(is.na(df$geo_dist)),
  missing_lineage = sum(is.na(df$Lineage) | df$Lineage == ""),
  missing_SB_plot = sum(is.na(df$SB_plot) | df$SB_plot == ""),
  missing_breed = sum(is.na(df$Cattle_Breed) | df$Cattle_Breed == ""),
  missing_region = sum(is.na(df$Region) | df$Region == ""),
  missing_district = sum(is.na(df$District) | df$District == "")
)

write.csv(qc_summary, "tables/data_quality_summary.csv", row.names = FALSE)
print(qc_summary)

# -----------------------------
# Publication theme
# -----------------------------

pub_theme <- theme_bw(base_size = 28) +
  theme(
    axis.title = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 24),
    strip.text = element_text(size = 24, face = "bold"),
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22)
  )

# -----------------------------
# Save helper
# -----------------------------

save_plot_png_svg <- function(plot_object,
                              filename_stub,
                              width = 14,
                              height = 10,
                              dpi = 500,
                              subfolder = NULL) {
  
  png_dir <- ifelse(is.null(subfolder),
                    "figures_png",
                    paste0("figures_png/", subfolder))
  
  svg_dir <- ifelse(is.null(subfolder),
                    "figures_svg",
                    paste0("figures_svg/", subfolder))
  
  dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(svg_dir, showWarnings = FALSE, recursive = TRUE)
  
  ggsave(
    filename = paste0(png_dir, "/", filename_stub, ".png"),
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi
  )
  
  ggsave(
    filename = paste0(svg_dir, "/", filename_stub, ".svg"),
    plot = plot_object,
    width = width,
    height = height
  )
}

# ============================================================
# ANALYSIS 1: Geographic distance vs SNP distance
# ============================================================

geo_snp_stats_plot <- function(var, filename_stub) {
  
  plot_df <- df %>%
    filter(
      !is.na(.data[[var]]),
      .data[[var]] != "",
      !is.na(snp_dist),
      !is.na(geo_dist)
    )
  
  if (nrow(plot_df) < 3 || n_distinct(plot_df[[var]]) < 1) {
    warning(paste("Skipping", var, "- insufficient data or no valid categories."))
    return(NULL)
  }
  
  stats <- plot_df %>%
    group_by(.data[[var]]) %>%
    filter(n() >= 3) %>%
    do({
      fit <- lm(geo_dist ~ snp_dist, data = .)
      
      tibble(
        n = nrow(.),
        slope = coef(fit)[2],
        intercept = coef(fit)[1],
        r_squared = summary(fit)$r.squared,
        adjusted_r_squared = summary(fit)$adj.r.squared,
        slope_pvalue = summary(fit)$coefficients[2, 4],
        spearman_rho = cor(
          .$snp_dist,
          .$geo_dist,
          method = "spearman",
          use = "complete.obs"
        ),
        spearman_pvalue = cor.test(
          .$snp_dist,
          .$geo_dist,
          method = "spearman",
          exact = FALSE
        )$p.value
      )
    }) %>%
    ungroup()
  
  write.csv(
    stats,
    paste0("statistics/slope_stats_geo_vs_snp_by_", var, ".csv"),
    row.names = FALSE
  )
  
  p <- ggplot(plot_df, aes(x = snp_dist, y = geo_dist)) +
    geom_point(alpha = 0.45, size = 3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "blue") +
    facet_wrap(as.formula(paste("~", var)), scales = "free") +
    pub_theme +
    labs(
      x = "Pairwise SNP distance",
      y = "Geographic distance (km)",
      title = paste("Geographic distance vs SNP distance by", var)
    )
  
  save_plot_png_svg(
    p,
    filename_stub,
    width = 14,
    height = 10,
    subfolder = "geo_vs_snp"
  )
  
  return(list(plot = p, stats = stats))
}

geo_lineage <- geo_snp_stats_plot("Lineage", "geo_vs_snp_by_lineage")
geo_sb <- geo_snp_stats_plot("SB_plot", "geo_vs_snp_by_SB")
geo_breed <- geo_snp_stats_plot("Cattle_Breed", "geo_vs_snp_by_cattle_breed")
geo_region <- geo_snp_stats_plot("Region", "geo_vs_snp_by_region")
geo_district <- geo_snp_stats_plot("District", "geo_vs_snp_by_district") # I removed the district level analysis.

geo_plots <- list(
  geo_lineage,
  geo_sb,
  geo_breed,
  geo_region
)

geo_plots <- geo_plots[!sapply(geo_plots, is.null)]

if (length(geo_plots) > 0) {
  
  combined_geo_vs_snp <- wrap_plots(
    lapply(geo_plots, function(x) x$plot),
    ncol = 2
  ) +
    plot_annotation(
      title = "Geographic distance versus pairwise SNP distance"
    )
  
  save_plot_png_svg(
    combined_geo_vs_snp,
    "combined_geo_vs_snp_plots",
    width = 36,
    height = 28,
    subfolder = "geo_vs_snp"
  )
}

# ============================================================
# ANALYSIS 2: Histogram of pairwise SNP distances
# ============================================================

hist_df <- df %>%
  filter(!is.na(snp_dist))

p_hist <- ggplot(hist_df, aes(x = snp_dist)) +
  geom_histogram(binwidth = 10, color = "black", fill = "grey80") +
  pub_theme +
  labs(
    x = "Pairwise SNP distance",
    y = "Frequency",
    title = "Distribution of pairwise SNP distances"
  )

save_plot_png_svg(
  p_hist,
  "histogram_pairwise_snps",
  width = 10,
  height = 8
)

# ============================================================
# ANALYSIS 3: SNP distance by categorical variables
# ============================================================

snp_category_plot <- function(var, filename_stub) {
  
  plot_df <- df %>%
    filter(
      !is.na(.data[[var]]),
      .data[[var]] != "",
      !is.na(snp_dist)
    )
  
  if (n_distinct(plot_df[[var]]) < 2) {
    warning(paste("Skipping", var, "- fewer than two categories available."))
    return(NULL)
  }
  
  kw <- kruskal_test(
    plot_df,
    as.formula(paste("snp_dist ~", var))
  )
  
  pairwise <- pairwise_wilcox_test(
    plot_df,
    as.formula(paste("snp_dist ~", var)),
    p.adjust.method = "BH",
    exact = FALSE
  )
  
  write.csv(
    kw,
    paste0("statistics/kruskal_stats_snp_by_", var, ".csv"),
    row.names = FALSE
  )
  
  write.csv(
    pairwise,
    paste0("statistics/pairwise_wilcox_stats_snp_by_", var, ".csv"),
    row.names = FALSE
  )
  
  summary_stats <- plot_df %>%
    group_by(.data[[var]]) %>%
    summarise(
      n_pairs = n(),
      median_snp = median(snp_dist),
      mean_snp = mean(snp_dist),
      IQR_snp = IQR(snp_dist),
      sd_snp = sd(snp_dist),
      min_snp = min(snp_dist),
      max_snp = max(snp_dist),
      .groups = "drop"
    )
  
  write.csv(
    summary_stats,
    paste0("statistics/summary_stats_snp_by_", var, ".csv"),
    row.names = FALSE
  )
  
  p <- ggplot(plot_df, aes(x = .data[[var]], y = snp_dist)) +
    geom_boxplot(outlier.shape = NA, linewidth = 1.2) +
    geom_jitter(width = 0.2, alpha = 0.25, size = 2.5) +
    pub_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = var,
      y = "Pairwise SNP distance",
      title = paste0(
        "Pairwise SNP distance by ", var,
        "\nKruskal-Wallis p = ", signif(kw$p, 3)
      )
    )
  
  save_plot_png_svg(
    p,
    filename_stub,
    width = 14,
    height = 10,
    subfolder = "boxplots"
  )
  
  return(list(plot = p, kruskal = kw, pairwise = pairwise, summary = summary_stats))
}

box_lineage <- snp_category_plot("Lineage", "pairwise_snps_by_lineage")
box_sb <- snp_category_plot("SB_plot", "pairwise_snps_by_SB")
box_breed <- snp_category_plot("Cattle_Breed", "pairwise_snps_by_cattle_breed")
box_region <- snp_category_plot("Region", "pairwise_snps_by_region")

box_plots <- list(
  box_lineage,
  box_sb,
  box_breed,
  box_region
)

box_plots <- box_plots[!sapply(box_plots, is.null)]

if (length(box_plots) > 0) {
  
  combined_boxplots <- wrap_plots(
    lapply(box_plots, function(x) x$plot),
    ncol = 2
  ) +
    plot_annotation(
      title = "Pairwise SNP distance distributions by metadata category"
    )
  
  save_plot_png_svg(
    combined_boxplots,
    "combined_pairwise_snp_boxplots",
    width = 30,
    height = 28,
    subfolder = "boxplots"
  )
}

# ============================================================
# ANALYSIS 4: Overall isolation-by-distance correlation
# ============================================================

ibd_df <- df %>%
  filter(!is.na(snp_dist), !is.na(geo_dist))

ibd_cor <- cor.test(
  ibd_df$geo_dist,
  ibd_df$snp_dist,
  method = "spearman",
  exact = FALSE
)

ibd_result <- data.frame(
  test = "Spearman correlation",
  rho = unname(ibd_cor$estimate),
  p_value = ibd_cor$p.value,
  n = nrow(ibd_df)
)

write.csv(
  ibd_result,
  "statistics/isolation_by_distance_spearman.csv",
  row.names = FALSE
)

p_ibd <- ggplot(ibd_df, aes(x = snp_dist, y = geo_dist)) +
  geom_point(alpha = 0.35, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "blue") +
  pub_theme +
  labs(
    x = "Pairwise SNP distance",
    y = "Geographic distance (km)",
    title = paste0(
      "Isolation-by-distance pattern\nSpearman rho = ",
      round(unname(ibd_cor$estimate), 3),
      ", p = ",
      signif(ibd_cor$p.value, 3)
    )
  )

save_plot_png_svg(
  p_ibd,
  "isolation_by_distance",
  width = 10,
  height = 8
)

# ============================================================
# ANALYSIS 5: Recent transmission pairs by SNP threshold
# ============================================================

recent_5 <- df %>%
  filter(!is.na(snp_dist), snp_dist <= 5)

recent_10 <- df %>%
  filter(!is.na(snp_dist), snp_dist <= 10)

recent_12 <- df %>%
  filter(!is.na(snp_dist), snp_dist <= 12)

write.csv(
  recent_5,
  "tables/recent_pairs_leq5SNPs.csv",
  row.names = FALSE
)

write.csv(
  recent_10,
  "tables/recent_pairs_leq10SNPs.csv",
  row.names = FALSE
)

write.csv(
  recent_12,
  "tables/recent_pairs_leq12SNPs.csv",
  row.names = FALSE
)

recent_summary_region <- recent_12 %>%
  filter(!is.na(Region), Region != "") %>%
  count(Region, sort = TRUE)

recent_summary_district <- recent_12 %>%
  filter(!is.na(District), District != "") %>%
  count(District, sort = TRUE)

write.csv(
  recent_summary_region,
  "tables/recent_pairs_leq12SNPs_by_region.csv",
  row.names = FALSE
)

write.csv(
  recent_summary_district,
  "tables/recent_pairs_leq12SNPs_by_district.csv",
  row.names = FALSE
)

# ============================================================
# ANALYSIS 6: Density plot helper
# ============================================================

density_plot_with_stats <- function(var, filename_stub) {
  
  plot_df <- df %>%
    filter(
      !is.na(.data[[var]]),
      .data[[var]] != "",
      !is.na(snp_dist)
    )
  
  if (n_distinct(plot_df[[var]]) < 2) {
    warning(paste("Skipping", var, "- fewer than two categories available."))
    return(NULL)
  }
  
  summary_stats <- plot_df %>%
    group_by(.data[[var]]) %>%
    summarise(
      n_pairs = n(),
      median_snp = median(snp_dist),
      mean_snp = mean(snp_dist),
      IQR_snp = IQR(snp_dist),
      sd_snp = sd(snp_dist),
      min_snp = min(snp_dist),
      max_snp = max(snp_dist),
      .groups = "drop"
    )
  
  kw <- kruskal_test(
    plot_df,
    as.formula(paste("snp_dist ~", var))
  )
  
  pairwise <- pairwise_wilcox_test(
    plot_df,
    as.formula(paste("snp_dist ~", var)),
    p.adjust.method = "BH",
    exact = FALSE
  )
  
  write.csv(
    summary_stats,
    paste0("statistics/density_summary_pairwise_snps_by_", var, ".csv"),
    row.names = FALSE
  )
  
  write.csv(
    kw,
    paste0("statistics/kruskal_density_pairwise_snps_by_", var, ".csv"),
    row.names = FALSE
  )
  
  write.csv(
    pairwise,
    paste0("statistics/pairwise_density_pairwise_snps_by_", var, ".csv"),
    row.names = FALSE
  )
  
  p <- ggplot(plot_df, aes(x = snp_dist, fill = .data[[var]])) +
    geom_density(alpha = 0.4) +
    pub_theme +
    labs(
      x = "Pairwise SNP distance",
      y = "Density",
      title = paste0(
        "Density distribution of pairwise SNP distances by ", var,
        "\nKruskal-Wallis p = ", signif(kw$p, 3)
      ),
      fill = var
    )
  
  save_plot_png_svg(
    p,
    filename_stub,
    width = 12,
    height = 8,
    subfolder = "density_plots"
  )
  
  return(list(plot = p, summary = summary_stats, kruskal = kw, pairwise = pairwise))
}

density_region <- density_plot_with_stats("Region", "density_pairwise_snps_by_region")
density_breed <- density_plot_with_stats("Cattle_Breed", "density_pairwise_snps_by_cattle_breed")
density_lineage <- density_plot_with_stats("Lineage", "density_pairwise_snps_by_lineage")
density_sb <- density_plot_with_stats("SB_plot", "density_pairwise_snps_by_SB")

density_plots <- list(
  density_lineage,
  density_sb,
  density_region,
  density_breed
)

density_plots <- density_plots[!sapply(density_plots, is.null)]

if (length(density_plots) > 0) {
  
  combined_density <- wrap_plots(
    lapply(density_plots, function(x) x$plot),
    ncol = 2
  ) +
    plot_annotation(
      title = "Density distribution of pairwise SNP distances"
    )
  
  save_plot_png_svg(
    combined_density,
    "combined_density_pairwise_snps",
    width = 24,
    height = 18,
    subfolder = "density_plots"
  )
}

# ============================================================
# ANALYSIS 7: Low SNP pairs by geographic distance, ≤12 SNPs
# ============================================================

low_snp_geo <- df %>%
  filter(!is.na(snp_dist), !is.na(geo_dist), snp_dist <= 12) %>%
  arrange(snp_dist, geo_dist)

write.csv(
  low_snp_geo,
  "tables/low_snp_pairs_with_geographic_distance_leq12.csv",
  row.names = FALSE
)

if (nrow(low_snp_geo) >= 3) {
  
  low_snp_lm <- lm(snp_dist ~ geo_dist, data = low_snp_geo)
  
  low_snp_lm_stats <- broom::tidy(low_snp_lm)
  low_snp_lm_glance <- broom::glance(low_snp_lm)
  
  low_snp_regression_summary <- data.frame(
    n_pairs = nrow(low_snp_geo),
    slope = coef(low_snp_lm)[2],
    intercept = coef(low_snp_lm)[1],
    r_squared = summary(low_snp_lm)$r.squared,
    adjusted_r_squared = summary(low_snp_lm)$adj.r.squared,
    slope_pvalue = summary(low_snp_lm)$coefficients[2, 4]
  )
  
  write.csv(
    low_snp_regression_summary,
    "statistics/low_snp_pairs_leq12_regression_summary.csv",
    row.names = FALSE
  )
  
  write.csv(
    low_snp_lm_stats,
    "statistics/low_snp_pairs_leq12_lm_coefficients.csv",
    row.names = FALSE
  )
  
  write.csv(
    low_snp_lm_glance,
    "statistics/low_snp_pairs_leq12_lm_model_fit.csv",
    row.names = FALSE
  )
  
  low_snp_spearman <- cor.test(
    low_snp_geo$geo_dist,
    low_snp_geo$snp_dist,
    method = "spearman",
    exact = FALSE
  )
  
  low_snp_spearman_summary <- data.frame(
    n_pairs = nrow(low_snp_geo),
    spearman_rho = unname(low_snp_spearman$estimate),
    p_value = low_snp_spearman$p.value
  )
  
  write.csv(
    low_snp_spearman_summary,
    "statistics/low_snp_pairs_leq12_spearman_summary.csv",
    row.names = FALSE
  )
  
  p_low_snp_geo <- ggplot(low_snp_geo, aes(x = geo_dist, y = snp_dist)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "blue") +
    pub_theme +
    labs(
      x = "Geographic distance (km)",
      y = "Pairwise SNP distance",
      title = paste0(
        "Geographic separation among ≤12 SNP isolate pairs\n",
        "Slope = ", round(coef(low_snp_lm)[2], 4),
        ", R² = ", round(summary(low_snp_lm)$r.squared, 3),
        ", Spearman ρ = ", round(unname(low_snp_spearman$estimate), 3),
        ", p = ", signif(low_snp_spearman$p.value, 3)
      )
    )
  
  save_plot_png_svg(
    p_low_snp_geo,
    "low_snp_pairs_geo_distance_leq12_with_slope",
    width = 10,
    height = 8,
    subfolder = "low_snp"
  )
}

# ============================================================
# ANALYSIS 8: Summary of low SNP pairs
# ============================================================

low_snp_summary <- recent_12 %>%
  summarise(
    total_pairs_leq12 = n(),
    median_geo_distance = median(geo_dist, na.rm = TRUE),
    mean_geo_distance = mean(geo_dist, na.rm = TRUE),
    min_geo_distance = min(geo_dist, na.rm = TRUE),
    max_geo_distance = max(geo_dist, na.rm = TRUE)
  )

write.csv(
  low_snp_summary,
  "tables/summary_low_snp_pairs_leq12.csv",
  row.names = FALSE
)

low_snp_by_lineage <- recent_12 %>%
  filter(!is.na(Lineage), Lineage != "") %>%
  count(Lineage, sort = TRUE)

low_snp_by_SB <- recent_12 %>%
  filter(!is.na(SB_plot), SB_plot != "") %>%
  count(SB_plot, sort = TRUE)

low_snp_by_breed <- recent_12 %>%
  filter(!is.na(Cattle_Breed), Cattle_Breed != "") %>%
  count(Cattle_Breed, sort = TRUE)

low_snp_by_region <- recent_12 %>%
  filter(!is.na(Region), Region != "") %>%
  count(Region, sort = TRUE)

low_snp_by_district <- recent_12 %>%
  filter(!is.na(District), District != "") %>%
  count(District, sort = TRUE)

write.csv(
  low_snp_by_lineage,
  "tables/low_snp_pairs_leq12_by_lineage.csv",
  row.names = FALSE
)

write.csv(
  low_snp_by_SB,
  "tables/low_snp_pairs_leq12_by_SB.csv",
  row.names = FALSE
)

write.csv(
  low_snp_by_breed,
  "tables/low_snp_pairs_leq12_by_cattle_breed.csv",
  row.names = FALSE
)

write.csv(
  low_snp_by_region,
  "tables/low_snp_pairs_leq12_by_region.csv",
  row.names = FALSE
)

write.csv(
  low_snp_by_district,
  "tables/low_snp_pairs_leq12_by_district.csv",
  row.names = FALSE
)

cat("Analysis complete.\n")
cat("PNG figures saved in: figures_png/\n")
cat("SVG figures saved in: figures_svg/\n")
cat("Statistics saved in: statistics/\n")
cat("Tables saved in: tables/\n")


