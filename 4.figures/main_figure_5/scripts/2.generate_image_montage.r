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
## Top coefficient feature
max_top_feat1 <- file.path(sc_crop_dir, "max_top_feature_1", "max_top_feature_1_composite_cropped.png")
max_top_feat2 <- file.path(sc_crop_dir, "max_top_feature_2", "max_top_feature_2_composite_cropped.png")
max_top_feat3 <- file.path(sc_crop_dir, "max_top_feature_3", "max_top_feature_3_composite_cropped.png")
max_top_feat4 <- file.path(sc_crop_dir, "max_top_feature_4", "max_top_feature_4_composite_cropped.png")
max_top_feat5 <- file.path(sc_crop_dir, "max_top_feature_5", "max_top_feature_5_composite_cropped.png")
max_top_feat6 <- file.path(sc_crop_dir, "max_top_feature_6", "max_top_feature_6_composite_cropped.png")
min_top_feat1 <- file.path(sc_crop_dir, "min_top_feature_1", "min_top_feature_1_composite_cropped.png")
min_top_feat2 <- file.path(sc_crop_dir, "min_top_feature_2", "min_top_feature_2_composite_cropped.png")
min_top_feat3 <- file.path(sc_crop_dir, "min_top_feature_3", "min_top_feature_3_composite_cropped.png")
min_top_feat4 <- file.path(sc_crop_dir, "min_top_feature_4", "min_top_feature_4_composite_cropped.png")
min_top_feat5 <- file.path(sc_crop_dir, "min_top_feature_5", "min_top_feature_5_composite_cropped.png")
min_top_feat6 <- file.path(sc_crop_dir, "min_top_feature_6", "min_top_feature_6_composite_cropped.png")

# Second top coefficient feature
max_second_top_feat1 <- file.path(sc_crop_dir, "max_second_top_feature_1", "max_second_top_feature_1_composite_cropped.png")
max_second_top_feat2 <- file.path(sc_crop_dir, "max_second_top_feature_2", "max_second_top_feature_2_composite_cropped.png")
max_second_top_feat3 <- file.path(sc_crop_dir, "max_second_top_feature_3", "max_second_top_feature_3_composite_cropped.png")
max_second_top_feat4 <- file.path(sc_crop_dir, "max_second_top_feature_4", "max_second_top_feature_4_composite_cropped.png")
max_second_top_feat5 <- file.path(sc_crop_dir, "max_second_top_feature_5", "max_second_top_feature_5_composite_cropped.png")
max_second_top_feat6 <- file.path(sc_crop_dir, "max_second_top_feature_6", "max_second_top_feature_6_composite_cropped.png")
min_second_top_feat1 <- file.path(sc_crop_dir, "min_second_top_feature_1", "min_second_top_feature_1_composite_cropped.png")
min_second_top_feat2 <- file.path(sc_crop_dir, "min_second_top_feature_2", "min_second_top_feature_2_composite_cropped.png")
min_second_top_feat3 <- file.path(sc_crop_dir, "min_second_top_feature_3", "min_second_top_feature_3_composite_cropped.png")
min_second_top_feat4 <- file.path(sc_crop_dir, "min_second_top_feature_4", "min_second_top_feature_4_composite_cropped.png")
min_second_top_feat5 <- file.path(sc_crop_dir, "min_second_top_feature_5", "min_second_top_feature_5_composite_cropped.png")
min_second_top_feat6 <- file.path(sc_crop_dir, "min_second_top_feature_6", "min_second_top_feature_6_composite_cropped.png")

# load top coefficient feat images 
max_top_feat1_image <- load_image(max_top_feat1)
max_top_feat2_image <- load_image(max_top_feat2)
max_top_feat3_image <- load_image(max_top_feat3)
max_top_feat4_image <- load_image(max_top_feat4)
max_top_feat5_image <- load_image(max_top_feat5)
max_top_feat6_image <- load_image(max_top_feat6)
min_top_feat1_image <- load_image(min_top_feat1)
min_top_feat2_image <- load_image(min_top_feat2)
min_top_feat3_image <- load_image(min_top_feat3)
min_top_feat4_image <- load_image(min_top_feat4)
min_top_feat5_image <- load_image(min_top_feat5)
min_top_feat6_image <- load_image(min_top_feat6)

# load second top coefficient feat images 
max_second_top_feat1_image <- load_image(max_second_top_feat1)
max_second_top_feat2_image <- load_image(max_second_top_feat2)
max_second_top_feat3_image <- load_image(max_second_top_feat3)
max_second_top_feat4_image <- load_image(max_second_top_feat4)
max_second_top_feat5_image <- load_image(max_second_top_feat5)
max_second_top_feat6_image <- load_image(max_second_top_feat6)
min_second_top_feat1_image <- load_image(min_second_top_feat1)
min_second_top_feat2_image <- load_image(min_second_top_feat2)
min_second_top_feat3_image <- load_image(min_second_top_feat3)
min_second_top_feat4_image <- load_image(min_second_top_feat4)
min_second_top_feat5_image <- load_image(min_second_top_feat5)
min_second_top_feat6_image <- load_image(min_second_top_feat6)

# Create list of images
list_of_images <- list(
    max_top_feat1_image,
    max_top_feat2_image,
    max_top_feat3_image,
    max_top_feat4_image,
    max_top_feat5_image,
    max_top_feat6_image,
    min_top_feat1_image,
    min_top_feat2_image,
    min_top_feat3_image,
    min_top_feat4_image,
    min_top_feat5_image,
    min_top_feat6_image,

    max_second_top_feat1_image,
    max_second_top_feat2_image,
    max_second_top_feat3_image,
    max_second_top_feat4_image,
    max_second_top_feat5_image,
    max_second_top_feat6_image,
    min_second_top_feat1_image,
    min_second_top_feat2_image,
    min_second_top_feat3_image,
    min_second_top_feat4_image,
    min_second_top_feat5_image,
    min_second_top_feat6_image
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
    + geom_text(aes(x = 0.5, y = 0.5, label = "WT A3 (NF1 +/+) genotype\nrepresentative single cells"), size = text_size, angle = 90) 
    + theme_void()
)
Null_text <- (
    ggplot()
    + geom_text(aes(x = 0.5, y = 0.5, label = "Null C04 (NF1 -/-) genotype\nrepresentative single cells"), size = text_size, angle = 90) 
    + theme_void()
)

# stich the images together for each top feature
## Top coefficient feature
max_top_feat_images_1 <- (
    list_of_images[[1]]
    + list_of_images[[2]]
    + list_of_images[[3]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
max_top_feat_images_2 <- (
    list_of_images[[4]]
    + list_of_images[[5]]
    + list_of_images[[6]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_top_feat_images_1 <- (
    list_of_images[[7]]
    + list_of_images[[8]]
    + list_of_images[[9]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_top_feat_images_2 <- (
    list_of_images[[10]]
    + list_of_images[[11]]
    + list_of_images[[12]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
## Second top coefficient feature
max_second_top_feat_images_1 <- (
    list_of_images[[13]]
    + list_of_images[[14]]
    + list_of_images[[15]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
max_second_top_feat_images_2 <- (
    list_of_images[[16]]
    + list_of_images[[17]]
    + list_of_images[[18]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_second_top_feat_images_1 <- (
    list_of_images[[19]]
    + list_of_images[[20]]
    + list_of_images[[21]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_second_top_feat_images_2 <- (
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
top_feat_text <- (
    ggplot()
    + geom_text(aes(x = 0.5, y = 0.5, label = "Highest model coefficient feature:\nDistribution of nuclei stain on edge of nuclei"), size = text_size) 
    + theme_void()
)
second_top_feat_text <- (
    ggplot()
    + geom_text(aes(x = 0.5, y = 0.5, label = "Second highest coefficient feature:\nDistribution of mito stain within nuclei"), size = text_size) 
    + theme_void()
)

# patch feature texts together
top_patch_text <- (
    top_feat_text
    + plot_layout(nrow = 1)
)
second_top_patch_text <- (
    second_top_feat_text
    + plot_layout(nrow = 1)
)

width <- 17
height <- 2.5

options(repr.plot.width = width, repr.plot.height = height)

top_patch_text


# Create montage
width <- 15
height <- 11

options(repr.plot.width = width, repr.plot.height = height)

top_feat_plot_max <- (
  WT_text + 
  (
    wrap_elements(full = top_patch_text) + 
    wrap_elements(max_top_feat_images_1) + 
    wrap_elements(max_top_feat_images_2) + 
    plot_layout(ncol = 1, heights = c(0.35, 1, 1))
  ) +
  plot_layout(widths = c(0.1, 1)) # Adjusts the width of the left area for the text
)

top_feat_plot_max

# save plot
ggsave(
    file.path(
        paste0(
            "./","nuclei_radial_feature_montage_max.png"
        )
    ),
    top_feat_plot_max, width = width, height = height, dpi = 600
)

# Create montage
width <- 15
height <- 11

options(repr.plot.width = width, repr.plot.height = height)

# patch the images together
top_feat_plot_min <- (
    Null_text + 
    (
        wrap_elements(full = blank) + 
        wrap_elements(min_top_feat_images_1) + 
        wrap_elements(min_top_feat_images_2) + 
        plot_layout(ncol = 1, heights = c(0.35, 1, 1))
    ) +
    plot_layout(widths = c(0.1, 1)) # Adjusts the width of the left area for the text
)

top_feat_plot_min

# save plot
ggsave(
    file.path(
        paste0(
            "./","nuclei_radial_feature_montage_min.png"
        )
    ),
    top_feat_plot_min, width = width, height = height, dpi = 600
)

# Create montage
width <- 15
height <- 11

options(repr.plot.width = width, repr.plot.height = height)

# patch the images together
second_top_feat_plot_max <- (
    blank + 
    (
        wrap_elements(full = second_top_patch_text) + 
        wrap_elements(max_second_top_feat_images_1) + 
        wrap_elements(max_second_top_feat_images_2) + 
        plot_layout(ncol = 1, heights = c(0.35, 1, 1))
    ) +
    plot_layout(widths = c(0.1, 1)) # Adjusts the width of the left area for the text
)

second_top_feat_plot_max

# save plot
ggsave(
    file.path(
        paste0(
            "./","mito_radial_feature_montage_max.png"
        )
    ),
    second_top_feat_plot_max, width = width, height = height, dpi = 600
)

# Create montage
width <- 15
height <- 11

options(repr.plot.width = width, repr.plot.height = height)

# patch the images together
second_top_feat_plot_min <- (
    blank + 
    (
        wrap_elements(full = blank) + 
        wrap_elements(min_second_top_feat_images_1) + 
        wrap_elements(min_second_top_feat_images_2) + 
        plot_layout(ncol = 1, heights = c(0.35, 1, 1))
    ) +
    plot_layout(widths = c(0.1, 1)) # Adjusts the width of the left area for the text
)

second_top_feat_plot_min

# save plot
ggsave(
    file.path(
        paste0(
            "./","mito_radial_feature_montage_min.png"
        )
    ),
    second_top_feat_plot_min, width = width, height = height, dpi = 600
)
