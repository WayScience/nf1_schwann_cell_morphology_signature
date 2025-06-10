suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(arrow))

figure_dir <- "../figures/supplementary"
output_supp_figure_8 <- file.path(
    figure_dir, "supp_figure_8_third_top_feature_montage.png"
)

top_feat_path_min = file.path("./cytoplasm_radial_feature_actin_montage_min.png")
top_feat_img_min = png::readPNG(top_feat_path_min)

# Get the dimensions of the image
img_height <- nrow(top_feat_img_min)
img_width <- ncol(top_feat_img_min)

# Calculate the aspect ratio
aspect_ratio <- img_height / img_width

# Plot the image montage to a ggplot object
top_feat_montage_min <- ggplot() +
  annotation_custom(
    rasterGrob(top_feat_img_min, interpolate = TRUE),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  theme_void() +
  coord_fixed(ratio = aspect_ratio, clip = "off") +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))  # Adjust margins as needed

top_feat_montage_min

top_feat_path_max = file.path("./cytoplasm_radial_feature_actin_montage_max.png")
top_feat_img_max = png::readPNG(top_feat_path_max)

# Get the dimensions of the image
img_height <- nrow(top_feat_img_max)
img_width <- ncol(top_feat_img_max)

# Calculate the aspect ratio
aspect_ratio <- img_height / img_width

# Plot the image montage to a ggplot object
top_feat_montage_max <- ggplot() +
  annotation_custom(
    rasterGrob(top_feat_img_max, interpolate = TRUE),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  theme_void() +
  coord_fixed(ratio = aspect_ratio, clip = "off") +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))  # Adjust margins as needed

top_feat_montage_max

image_montage <- (
   free(top_feat_montage_min) +
   top_feat_montage_max
) + plot_layout(heights = c(1,1), guides = "collect")

image_montage


supp_fig_8_gg <- (
  image_montage
) + plot_annotation(tag_levels = list(c("A", "B", "", ""))) & theme(plot.tag = element_text(size = 30))

# Save or display the plot
ggsave(output_supp_figure_8, plot = supp_fig_8_gg, dpi = 500, height = 16, width = 12)
