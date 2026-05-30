getwd()
setwd("D:/Dropbox/bovineWGSA/BTBGhAnalysis/maps")

# ============================================================
# Ghana maps for M. bovis sampling points
# Outputs:
# 1. Map showing sampled districts using District_midpoints sheet
# 2. Jittered isolate point map with district labels
# 3. District point-size map by number of isolates
# 4. District pie maps for SB, Host, Lineage, Cattle_Breed
# ============================================================

packages <- c(
  "tidyverse",
  "readxl",
  "sf",
  "terra",
  "geodata",
  "ggplot2",
  "ggrepel",
  "scatterpie",
  "scales"
)

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)
library(readxl)
library(sf)
library(terra)
library(geodata)
library(ggplot2)
library(ggrepel)
library(scatterpie)
library(scales)

# -----------------------------
# Output folders
# -----------------------------

dir.create("figures_png", showWarnings = FALSE)
dir.create("figures_svg", showWarnings = FALSE)
dir.create("spatial_data", showWarnings = FALSE)

# -----------------------------
# Load sample coordinates
# -----------------------------

df_map <- read_excel(
  "Ghana_Mbovis_pairwise_geographic_SNP_distance_merged_REVISED.xlsx",
  sheet = "Sample_coordinates"
)

df_map <- df_map %>%
  mutate(
    Latitude_plot = Latitude,
    Longitude_plot = Longitude
  ) %>%
  filter(
    !is.na(Latitude_plot),
    !is.na(Longitude_plot),
    !is.na(District),
    District != "NA",
    District != ""
  ) %>%
  mutate(
    Region = ifelse(is.na(Region) | Region == "", "Unknown", Region),
    District = as.character(District),
    SB_number = ifelse(is.na(SB) | SB == "", "Unknown", SB),
    Host = ifelse(is.na(Host) | Host == "", "Unknown", Host),
    Lineage = ifelse(is.na(sub_lineage) | sub_lineage == "", "Unknown", sub_lineage),
    Cattle_Breed = ifelse(is.na(Cattle_Breed) | Cattle_Breed == "", "Unknown", Cattle_Breed)
  )

# -----------------------------
# Load district midpoint sheet
# -----------------------------

district_midpoints <- read_excel(
  "Ghana_Mbovis_pairwise_geographic_SNP_distance_merged_REVISED.xlsx",
  sheet = "District_midpoints"
)

district_midpoints <- district_midpoints %>%
  filter(
    !is.na(Latitude),
    !is.na(Longitude),
    !is.na(District),
    District != "",
    District != "NA"
  ) %>%
  mutate(
    District = as.character(District),
    Region = ifelse(is.na(Region) | Region == "", "Unknown", Region)
  )

# -----------------------------
# Ghana administrative boundaries
# Level 1 = regions
# Level 2 = districts
# -----------------------------

ghana_regions <- geodata::gadm(
  country = "GHA",
  level = 1,
  path = "spatial_data"
)

ghana_districts <- geodata::gadm(
  country = "GHA",
  level = 2,
  path = "spatial_data"
)

ghana_regions_sf <- st_as_sf(ghana_regions)
ghana_districts_sf <- st_as_sf(ghana_districts)

# -----------------------------
# Common map theme
# -----------------------------

map_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.key = element_rect(fill = "white", color = NA)
  )

# -----------------------------
# Colour palettes
# -----------------------------

sb_palette <- c(
  "SB0944" = "#6B8E23",
  "SB0300" = "#FFD92F",
  "SB0134" = "#FBB4B9",
  "SB0878" = "#00E676",
  "SB1025" = "#008000",
  "SB1026" = "#8B0000",
  "SB1027" = "#4DD0E1",
  "SB1410" = "#FF33CC",
  "SB1418" = "#40E0D0",
  "SB1432" = "#F28E2B",
  "SB1439" = "#4B3CFA",
  "SB1472" = "#E91E63",
  "SB1517" = "#FFFF33",
  "SB2286" = "#2ECC40",
  "SB2745" = "#3A33FF",
  "SB2749" = "#66E52B",
  "SB2750" = "#EC2A8C",
  "SB2751" = "#35C9E8",
  "SB2753" = "#F4D03F",
  "SB2756" = "#8A2BE2",
  "SB2757" = "#2EE65A",
  "SB2758" = "#FF2B2B",
  "SB2762" = "#2B6FF6",
  "SB2763" = "#99FF33",
  "SB2764" = "#E91EDB",
  "SB2899" = "#F92672",
  "SB2900" = "#3498DB",
  "SB2901" = "#D6F21A",
  "SB2902" = "#C026E8",
  "Unknown" = "grey80"
)

host_palette <- c(
  "Bovine Ghana" = "#2CA02C",
  "Bovine" = "#2CA02C",
  "Human Ghana" = "#E377C2",
  "Human" = "#E377C2",
  "Unknown" = "grey80"
)

lineage_palette <- c(
  "La1.6" = "#8C564B",
  "La1.8.2" = "#39346E",
  "Unknown" = "grey80"
)

breed_palette <- c(
  "Ndama" = "#1F77B4",
  "Sanga" = "#FF7F0E",
  "WASH" = "#2CA02C",
  "Unknown" = "grey80"
)

region_palette <- c(
  "Greater Accra" = "#1F77B4",
  "North_East" = "#FF7F0E",
  "Northern" = "#2CA02C",
  "Savannah" = "#D62728",
  "Upper_West" = "#9467BD",
  "MD" = "grey60",
  "Unknown" = "grey80"
)

get_palette <- function(variable, values) {
  
  if (variable == "SB_number") {
    pal <- sb_palette
  } else if (variable == "Host") {
    pal <- host_palette
  } else if (variable == "Lineage") {
    pal <- lineage_palette
  } else if (variable == "Cattle_Breed") {
    pal <- breed_palette
  } else if (variable == "Region") {
    pal <- region_palette
  } else {
    pal <- setNames(scales::hue_pal()(length(values)), values)
  }
  
  missing_cols <- setdiff(values, names(pal))
  
  if (length(missing_cols) > 0) {
    extra_cols <- setNames(
      scales::hue_pal()(length(missing_cols)),
      missing_cols
    )
    pal <- c(pal, extra_cols)
  }
  
  pal[values]
}

# ============================================================
# 1. Map showing sampled districts using District_midpoints
# ============================================================

district_midpoints_sf <- st_as_sf(
  district_midpoints,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

district_colour_palette <- setNames(
  scales::hue_pal()(length(sort(unique(district_midpoints$District)))),
  sort(unique(district_midpoints$District))
)

p_sampled_districts <- ggplot() +
  geom_sf(
    data = ghana_districts_sf,
    fill = "grey97",
    color = "grey75",
    linewidth = 0.25
  ) +
  geom_sf(
    data = ghana_regions_sf,
    fill = NA,
    color = "black",
    linewidth = 0.7
  ) +
  geom_sf(
    data = district_midpoints_sf,
    aes(color = District),
    size = 4,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = district_midpoints,
    aes(
      x = Longitude,
      y = Latitude,
      label = District
    ),
    size = 3.2,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey40"
  ) +
  scale_color_manual(values = district_colour_palette) +
  coord_sf(
    xlim = c(-3.6, 1.5),
    ylim = c(4.5, 11.5),
    expand = FALSE
  ) +
  map_theme +
  labs(
    title = "Districts represented in the M. bovis dataset",
    x = "Longitude",
    y = "Latitude",
    color = "Sampled district"
  )

ggsave(
  "figures_png/ghana_sampled_districts_midpoints.png",
  p_sampled_districts,
  width = 9,
  height = 10,
  dpi = 500
)

ggsave(
  "figures_svg/ghana_sampled_districts_midpoints.svg",
  p_sampled_districts,
  width = 9,
  height = 10
)

# ============================================================
# 2. Jittered isolate-level point map with district labels
# ============================================================

set.seed(123)

df_jitter <- df_map %>%
  mutate(
    Longitude_jitter = jitter(Longitude_plot, amount = 0.05),
    Latitude_jitter = jitter(Latitude_plot, amount = 0.05)
  )

points_jitter_sf <- st_as_sf(
  df_jitter,
  coords = c("Longitude_jitter", "Latitude_jitter"),
  crs = 4326,
  remove = FALSE
)

district_label_df <- df_map %>%
  group_by(District, Region, Latitude_plot, Longitude_plot) %>%
  summarise(n_isolates = n(), .groups = "drop")

p_jitter <- ggplot() +
  geom_sf(
    data = ghana_districts_sf,
    fill = "grey97",
    color = "grey75",
    linewidth = 0.25
  ) +
  geom_sf(
    data = ghana_regions_sf,
    fill = NA,
    color = "black",
    linewidth = 0.7
  ) +
  geom_sf(
    data = points_jitter_sf,
    aes(color = Region),
    size = 2,
    alpha = 0.8
  ) +
  geom_text_repel(
    data = district_label_df,
    aes(
      x = Longitude_plot,
      y = Latitude_plot,
      label = District
    ),
    size = 3.2,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey40"
  ) +
  scale_color_manual(values = region_palette) +
  coord_sf(
    xlim = c(-3.6, 1.5),
    ylim = c(4.5, 11.5),
    expand = FALSE
  ) +
  map_theme +
  labs(
    title = "Jittered isolate locations by district midpoint",
    x = "Longitude",
    y = "Latitude",
    color = "Region"
  )

ggsave(
  "figures_png/ghana_jittered_isolate_points.png",
  p_jitter,
  width = 8,
  height = 10,
  dpi = 500
)

ggsave(
  "figures_svg/ghana_jittered_isolate_points.svg",
  p_jitter,
  width = 8,
  height = 10
)

# ============================================================
# 3. District-level point-size map by number of isolates
# ============================================================

district_counts <- df_map %>%
  group_by(District, Region, Latitude_plot, Longitude_plot) %>%
  summarise(
    n_isolates = n(),
    .groups = "drop"
  )

district_counts_sf <- st_as_sf(
  district_counts,
  coords = c("Longitude_plot", "Latitude_plot"),
  crs = 4326,
  remove = FALSE
)

p_counts <- ggplot() +
  geom_sf(
    data = ghana_districts_sf,
    fill = "grey97",
    color = "grey75",
    linewidth = 0.25
  ) +
  geom_sf(
    data = ghana_regions_sf,
    fill = NA,
    color = "black",
    linewidth = 0.7
  ) +
  geom_sf(
    data = district_counts_sf,
    aes(size = n_isolates, color = Region),
    alpha = 0.85
  ) +
  geom_text_repel(
    data = district_counts,
    aes(
      x = Longitude_plot,
      y = Latitude_plot,
      label = District
    ),
    size = 3.2,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey40"
  ) +
  scale_size_continuous(range = c(1.8, 6)) +
  scale_color_manual(values = region_palette) +
  coord_sf(
    xlim = c(-3.6, 1.5),
    ylim = c(4.5, 11.5),
    expand = FALSE
  ) +
  map_theme +
  labs(
    title = "Number of M. bovis isolates by district midpoint",
    x = "Longitude",
    y = "Latitude",
    color = "Region",
    size = "No. isolates"
  )

ggsave(
  "figures_png/ghana_district_isolate_counts.png",
  p_counts,
  width = 8,
  height = 10,
  dpi = 500
)

ggsave(
  "figures_svg/ghana_district_isolate_counts.svg",
  p_counts,
  width = 8,
  height = 10
)

# ============================================================
# Helper function for district pie maps
# Pie radius reflects total isolates, while size legend uses
# point-style legend matching the point-size map.
# ============================================================

make_district_pie_map <- function(data,
                                  variable,
                                  filename_stub,
                                  title_text,
                                  top_n = NULL,
                                  show_size_legend = TRUE) {
  
  var_sym <- sym(variable)
  
  data2 <- data %>%
    filter(
      !is.na(!!var_sym),
      !!var_sym != ""
    )
  
  if (!is.null(top_n)) {
    dominant_categories <- data2 %>%
      count(!!var_sym, sort = TRUE) %>%
      slice_head(n = top_n) %>%
      pull(!!var_sym)
    
    data2 <- data2 %>%
      mutate(
        plot_group = ifelse(
          !!var_sym %in% dominant_categories,
          as.character(!!var_sym),
          "Other"
        )
      )
    
    variable2 <- "plot_group"
    var2_sym <- sym(variable2)
  } else {
    variable2 <- variable
    var2_sym <- var_sym
  }
  
  category_order <- data2 %>%
    count(!!var2_sym, sort = TRUE) %>%
    pull(!!var2_sym) %>%
    as.character()
  
  if (variable == "SB_number") {
    category_order <- c(
      intersect(c("SB0944", "SB0300"), category_order),
      setdiff(category_order, c("SB0944", "SB0300"))
    )
  }
  
  pie_df <- data2 %>%
    mutate(
      plot_var = factor(as.character(!!var2_sym), levels = category_order)
    ) %>%
    group_by(District, Latitude_plot, Longitude_plot, plot_var) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = plot_var,
      values_from = n,
      values_fill = 0
    )
  
  pie_cols <- setdiff(
    names(pie_df),
    c("District", "Latitude_plot", "Longitude_plot")
  )
  
  pie_cols <- category_order[category_order %in% pie_cols]
  
  min_radius <- 0.045
  max_radius <- 0.14
  
  radius_from_count <- function(x, count_range) {
    if (diff(count_range) == 0) {
      rep((min_radius + max_radius) / 2, length(x))
    } else {
      scales::rescale(
        sqrt(x),
        from = sqrt(count_range),
        to = c(min_radius, max_radius)
      )
    }
  }
  
  pie_df <- pie_df %>%
    mutate(
      total_isolates = rowSums(across(all_of(pie_cols)))
    )
  
  count_range <- range(pie_df$total_isolates)
  
  pie_df <- pie_df %>%
    mutate(
      pie_radius = radius_from_count(total_isolates, count_range)
    )
  
  set.seed(123)
  
  pie_df <- pie_df %>%
    mutate(
      Longitude_pie = Longitude_plot + jitter(rep(0, n()), amount = 0.025),
      Latitude_pie = Latitude_plot + jitter(rep(0, n()), amount = 0.025)
    )
  
  pie_palette <- get_palette(variable, pie_cols)
  
  p <- ggplot() +
    geom_sf(
      data = ghana_districts_sf,
      fill = "grey97",
      color = "grey75",
      linewidth = 0.25
    ) +
    geom_sf(
      data = ghana_regions_sf,
      fill = NA,
      color = "black",
      linewidth = 0.7
    ) +
    geom_scatterpie(
      data = pie_df,
      aes(
        x = Longitude_pie,
        y = Latitude_pie,
        r = pie_radius
      ),
      cols = pie_cols,
      alpha = 0.88,
      color = "black",
      linewidth = 0.25
    ) +
    geom_text_repel(
      data = pie_df,
      aes(
        x = Longitude_pie,
        y = Latitude_pie,
        label = District
      ),
      size = 3,
      max.overlaps = 100,
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey40"
    ) +
    scale_fill_manual(
      values = pie_palette,
      breaks = pie_cols
    ) +
    coord_sf(
      xlim = c(-3.6, 1.5),
      ylim = c(4.5, 11.5),
      expand = FALSE
    ) +
    map_theme +
    labs(
      title = title_text,
      x = "Longitude",
      y = "Latitude",
      fill = variable
    )
  
  if (show_size_legend) {
    
    size_values <- pretty(pie_df$total_isolates, n = 4)
    size_values <- size_values[size_values > 0]
    size_values <- unique(round(size_values))
    
    # Keep only observed range
    size_values <- size_values[
      size_values >= min(pie_df$total_isolates) &
        size_values <= max(pie_df$total_isolates)
    ]
    
    legend_radii <- radius_from_count(size_values, count_range)
    
    size_legend_df <- data.frame(
      n_isolates = size_values,
      radius = legend_radii
    )
    
    write.csv(
      size_legend_df,
      paste0(filename_stub, "_size_legend_values.csv"),
      row.names = FALSE
    )
    
    p <- p +
      scatterpie::geom_scatterpie_legend(
        r = legend_radii,
        x = -3.35,
        y = 4.85,
        n = length(legend_radii),
        labeller = function(x) {
          sapply(x, function(rr) {
            closest <- size_legend_df$n_isolates[
              which.min(abs(size_legend_df$radius - rr))
            ]
            paste0(closest, " isolates")
          })
        }
      )
  }
  
  print(p)
  
  ggsave(
    paste0("figures_png/", filename_stub, ".png"),
    p,
    width = 10,
    height = 10,
    dpi = 500
  )
  
  ggsave(
    paste0("figures_svg/", filename_stub, ".svg"),
    p,
    width = 10,
    height = 10
  )
  
  write.csv(
    pie_df,
    paste0(filename_stub, "_district_counts.csv"),
    row.names = FALSE
  )
  
  return(p)
}

# ============================================================
# 4. District pie chart map by SB number: ALL SBs shown
# ============================================================

p_sb_pie <- make_district_pie_map(
  data = df_map,
  variable = "SB_number",
  filename_stub = "ghana_district_pie_SB_number_all_SBs",
  title_text = "District-level SB number composition of M. bovis isolates",
  top_n = NULL,
  show_size_legend = TRUE
)

# ============================================================
# 5. District pie chart map by Host
# ============================================================

p_host_pie <- make_district_pie_map(
  data = df_map,
  variable = "Host",
  filename_stub = "ghana_district_pie_Host",
  title_text = "District-level host composition of M. bovis isolates",
  show_size_legend = TRUE
)

# ============================================================
# 6. District pie chart map by Lineage
# ============================================================

p_lineage_pie <- make_district_pie_map(
  data = df_map,
  variable = "Lineage",
  filename_stub = "ghana_district_pie_Lineage",
  title_text = "District-level lineage composition of M. bovis isolates",
  show_size_legend = TRUE
)

# ============================================================
# 7. District pie chart map by Cattle_Breed
# ============================================================

p_breed_pie <- make_district_pie_map(
  data = df_map,
  variable = "Cattle_Breed",
  filename_stub = "ghana_district_pie_Cattle_Breed",
  title_text = "District-level cattle breed composition of M. bovis isolates",
  show_size_legend = TRUE
)

cat("Mapping complete.\n")
cat("PNG figures saved in figures_png/\n")
cat("SVG figures saved in figures_svg/\n")