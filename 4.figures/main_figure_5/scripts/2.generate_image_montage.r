suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(RColorBrewer))

load_image <- function(path){
    img <- png::readPNG(path)
    # Convert the image to a raster object
    g <- grid::rasterGrob(img, interpolate=TRUE)

    # Create a ggplot
    p <- ggplot() +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme_void()
    return(p)
}

# Directory with single-cell crops
sc_crop_dir <- "./sc_crops"

# Path to each composite image (min or max) per top feature
## Top null feature (actin radial distribution)
max_radial_feat1 <- file.path(sc_crop_dir, "max_radial_feature_1", "max_radial_feature_1_composite_cropped.png")
max_radial_feat2 <- file.path(sc_crop_dir, "max_radial_feature_2", "max_radial_feature_2_composite_cropped.png")
max_radial_feat3 <- file.path(sc_crop_dir, "max_radial_feature_3", "max_radial_feature_3_composite_cropped.png")
max_radial_feat4 <- file.path(sc_crop_dir, "max_radial_feature_4", "max_radial_feature_4_composite_cropped.png")
max_radial_feat5 <- file.path(sc_crop_dir, "max_radial_feature_5", "max_radial_feature_5_composite_cropped.png")
max_radial_feat6 <- file.path(sc_crop_dir, "max_radial_feature_6", "max_radial_feature_6_composite_cropped.png")
min_radial_feat1 <- file.path(sc_crop_dir, "min_radial_feature_1", "min_radial_feature_1_composite_cropped.png")
min_radial_feat2 <- file.path(sc_crop_dir, "min_radial_feature_2", "min_radial_feature_2_composite_cropped.png")
min_radial_feat3 <- file.path(sc_crop_dir, "min_radial_feature_3", "min_radial_feature_3_composite_cropped.png")
min_radial_feat4 <- file.path(sc_crop_dir, "min_radial_feature_4", "min_radial_feature_4_composite_cropped.png")
min_radial_feat5 <- file.path(sc_crop_dir, "min_radial_feature_5", "min_radial_feature_5_composite_cropped.png")
min_radial_feat6 <- file.path(sc_crop_dir, "min_radial_feature_6", "min_radial_feature_6_composite_cropped.png")

# Top WT feature (DAPI + ER cell correlation)
max_corr_feat1 <- file.path(sc_crop_dir, "max_corr_feature_1", "max_corr_feature_1_composite_cropped.png")
max_corr_feat2 <- file.path(sc_crop_dir, "max_corr_feature_2", "max_corr_feature_2_composite_cropped.png")
max_corr_feat3 <- file.path(sc_crop_dir, "max_corr_feature_3", "max_corr_feature_3_composite_cropped.png")
max_corr_feat4 <- file.path(sc_crop_dir, "max_corr_feature_4", "max_corr_feature_4_composite_cropped.png")
max_corr_feat5 <- file.path(sc_crop_dir, "max_corr_feature_5", "max_corr_feature_5_composite_cropped.png")
max_corr_feat6 <- file.path(sc_crop_dir, "max_corr_feature_6", "max_corr_feature_6_composite_cropped.png")
min_corr_feat1 <- file.path(sc_crop_dir, "min_corr_feature_1", "min_corr_feature_1_composite_cropped.png")
min_corr_feat2 <- file.path(sc_crop_dir, "min_corr_feature_2", "min_corr_feature_2_composite_cropped.png")
min_corr_feat3 <- file.path(sc_crop_dir, "min_corr_feature_3", "min_corr_feature_3_composite_cropped.png")
min_corr_feat4 <- file.path(sc_crop_dir, "min_corr_feature_4", "min_corr_feature_4_composite_cropped.png")
min_corr_feat5 <- file.path(sc_crop_dir, "min_corr_feature_5", "min_corr_feature_5_composite_cropped.png")
min_corr_feat6 <- file.path(sc_crop_dir, "min_corr_feature_6", "min_corr_feature_6_composite_cropped.png")

# load top radial feat images 
max_radial_feat1_image <- load_image(max_radial_feat1)
max_radial_feat2_image <- load_image(max_radial_feat2)
max_radial_feat3_image <- load_image(max_radial_feat3)
max_radial_feat4_image <- load_image(max_radial_feat4)
max_radial_feat5_image <- load_image(max_radial_feat5)
max_radial_feat6_image <- load_image(max_radial_feat6)
min_radial_feat1_image <- load_image(min_radial_feat1)
min_radial_feat2_image <- load_image(min_radial_feat2)
min_radial_feat3_image <- load_image(min_radial_feat3)
min_radial_feat4_image <- load_image(min_radial_feat4)
min_radial_feat5_image <- load_image(min_radial_feat5)
min_radial_feat6_image <- load_image(min_radial_feat6)

# load top corr feat images 
max_corr_feat1_image <- load_image(max_corr_feat1)
max_corr_feat2_image <- load_image(max_corr_feat2)
max_corr_feat3_image <- load_image(max_corr_feat3)
max_corr_feat4_image <- load_image(max_corr_feat4)
max_corr_feat5_image <- load_image(max_corr_feat5)
max_corr_feat6_image <- load_image(max_corr_feat6)
min_corr_feat1_image <- load_image(min_corr_feat1)
min_corr_feat2_image <- load_image(min_corr_feat2)
min_corr_feat3_image <- load_image(min_corr_feat3)
min_corr_feat4_image <- load_image(min_corr_feat4)
min_corr_feat5_image <- load_image(min_corr_feat5)
min_corr_feat6_image <- load_image(min_corr_feat6)

# Create list of images
list_of_images <- list(
    max_radial_feat1_image,
    max_radial_feat2_image,
    max_radial_feat3_image,
    max_radial_feat4_image,
    max_radial_feat5_image,
    max_radial_feat6_image,
    min_radial_feat1_image,
    min_radial_feat2_image,
    min_radial_feat3_image,
    min_radial_feat4_image,
    min_radial_feat5_image,
    min_radial_feat6_image,

    max_corr_feat1_image,
    max_corr_feat2_image,
    max_corr_feat3_image,
    max_corr_feat4_image,
    max_corr_feat5_image,
    max_corr_feat6_image,
    min_corr_feat1_image,
    min_corr_feat2_image,
    min_corr_feat3_image,
    min_corr_feat4_image,
    min_corr_feat5_image,
    min_corr_feat6_image
)

width <- 2.5
height <- 2.5

text_size <- 13

options(repr.plot.width = width, repr.plot.height = height)

# blank
blank <- (
    ggplot()
    + geom_text(aes(x = 0.5, y = 0.5, label = ""), size = text_size) 
    + theme_void()
)

# ggplot of just text for labelling y axis
WT_text <- (
    ggplot()
    + geom_text(aes(x = 0.5, y = 0.5, label = "WT (NF1 +/+) genotype\nrepresentative single cells"), size = text_size, angle = 90) 
    + theme_void()
)
Null_text <- (
    ggplot()
    + geom_text(aes(x = 0.5, y = 0.5, label = "Null (NF1 -/-) genotype\nrepresentative single cells"), size = text_size, angle = 90) 
    + theme_void()
)

# stich the images together for each top feature
## top features for Null cells
max_radial_feat_images_1 <- (
    list_of_images[[1]]
    + list_of_images[[2]]
    + list_of_images[[3]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
max_radial_feat_images_2 <- (
    list_of_images[[4]]
    + list_of_images[[5]]
    + list_of_images[[6]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_radial_feat_images_1 <- (
    list_of_images[[7]]
    + list_of_images[[8]]
    + list_of_images[[9]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_radial_feat_images_2 <- (
    list_of_images[[10]]
    + list_of_images[[11]]
    + list_of_images[[12]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
## Top WT feature
max_corr_feat_images_1 <- (
    list_of_images[[13]]
    + list_of_images[[14]]
    + list_of_images[[15]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
max_corr_feat_images_2 <- (
    list_of_images[[16]]
    + list_of_images[[17]]
    + list_of_images[[18]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_corr_feat_images_1 <- (
    list_of_images[[19]]
    + list_of_images[[20]]
    + list_of_images[[21]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_corr_feat_images_2 <- (
    list_of_images[[22]]
    + list_of_images[[23]]
    + list_of_images[[24]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)

# Generate labels for each plot with CellProfiler feature
width <- 2.5
height <- 2.5

text_size <- 14

options(repr.plot.width = width, repr.plot.height = height)

# ggplot of just text
radial_feat_text <- (
    ggplot()
    + geom_text(aes(x = 0.5, y = 0.5, label = "Top Null (NF1 -/-) predicting feature:\nIntensity of F-actin at the edge of cytoskeleton"), size = text_size) 
    + theme_void()
)
corr_feat_text <- (
    ggplot()
    + geom_text(aes(x = 0.5, y = 0.5, label = "Top WT (NF1 +/+) predicting feature:\nCorrelation of ER and nucleus stains"), size = text_size) 
    + theme_void()
)

# patch feature texts together
radial_patch_text <- (
    radial_feat_text
    + plot_layout(nrow = 1)
)
corr_patch_text <- (
    corr_feat_text
    + plot_layout(nrow = 1)
)

width <- 17
height <- 2.5

options(repr.plot.width = width, repr.plot.height = height)

radial_patch_text


# Create montage
width <- 15
height <- 11

options(repr.plot.width = width, repr.plot.height = height)

radial_feat_plot_max <- (
  Null_text + 
  (
    wrap_elements(full = radial_patch_text) + 
    wrap_elements(max_radial_feat_images_1) + 
    wrap_elements(max_radial_feat_images_2) + 
    plot_layout(ncol = 1, heights = c(0.35, 1, 1))
  ) +
  plot_layout(widths = c(0.1, 1)) # Adjusts the width of the left area for the text
)

radial_feat_plot_max

# save plot
ggsave(
    file.path(
        paste0(
            "./","radial_feature_montage_max.png"
        )
    ),
    radial_feat_plot_max, width = width, height = height, dpi = 600
)

# Create montage
width <- 15
height <- 11

options(repr.plot.width = width, repr.plot.height = height)

# patch the images together
radial_feat_plot_min <- (
    WT_text + 
    (
        wrap_elements(full = radial_patch_text) + 
        wrap_elements(min_radial_feat_images_1) + 
        wrap_elements(min_radial_feat_images_2) + 
        plot_layout(ncol = 1, heights = c(0.35, 1, 1))
    ) +
    plot_layout(widths = c(0.1, 1)) # Adjusts the width of the left area for the text
)

radial_feat_plot_min

# save plot
ggsave(
    file.path(
        paste0(
            "./","radial_feature_montage_min.png"
        )
    ),
    radial_feat_plot_min, width = width, height = height, dpi = 600
)

# Create montage
width <- 15
height <- 11

options(repr.plot.width = width, repr.plot.height = height)

# patch the images together
corr_feat_plot_max <- (
    WT_text + 
    (
        wrap_elements(full = corr_patch_text) + 
        wrap_elements(max_corr_feat_images_1) + 
        wrap_elements(max_corr_feat_images_2) + 
        plot_layout(ncol = 1, heights = c(0.35, 1, 1))
    ) +
    plot_layout(widths = c(0.1, 1)) # Adjusts the width of the left area for the text
)

corr_feat_plot_max

# save plot
ggsave(
    file.path(
        paste0(
            "./","correlation_feature_montage_max.png"
        )
    ),
    corr_feat_plot_max, width = width, height = height, dpi = 600
)

# Create montage
width <- 15
height <- 11

options(repr.plot.width = width, repr.plot.height = height)

# patch the images together
corr_feat_plot_min <- (
    Null_text + 
    (
        wrap_elements(full = corr_patch_text) + 
        wrap_elements(min_corr_feat_images_1) + 
        wrap_elements(min_corr_feat_images_2) + 
        plot_layout(ncol = 1, heights = c(0.35, 1, 1))
    ) +
    plot_layout(widths = c(0.1, 1)) # Adjusts the width of the left area for the text
)

corr_feat_plot_min

# save plot
ggsave(
    file.path(
        paste0(
            "./","correlation_feature_montage_min.png"
        )
    ),
    corr_feat_plot_min, width = width, height = height, dpi = 600
)
