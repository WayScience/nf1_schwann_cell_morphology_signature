suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(magick))

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

# Top coefficient feature
max_top_feat_orig_1 <- file.path(sc_crop_dir, "max_top_feature_orig_1", "max_top_feature_orig_1_composite_cropped_annotated.png")
max_top_feat_orig_2 <- file.path(sc_crop_dir, "max_top_feature_orig_2", "max_top_feature_orig_2_composite_cropped_annotated.png")
max_top_feat_orig_3 <- file.path(sc_crop_dir, "max_top_feature_orig_3", "max_top_feature_orig_3_composite_cropped_annotated.png")
max_top_feat_deriv_1 <- file.path(sc_crop_dir, "max_top_feature_deriv_1", "max_top_feature_deriv_1_composite_cropped_annotated.png")
max_top_feat_deriv_2 <- file.path(sc_crop_dir, "max_top_feature_deriv_2", "max_top_feature_deriv_2_composite_cropped_annotated.png")
max_top_feat_deriv_3 <- file.path(sc_crop_dir, "max_top_feature_deriv_3", "max_top_feature_deriv_3_composite_cropped_annotated.png")

min_top_feat_orig_1 <- file.path(sc_crop_dir, "min_top_feature_orig_1", "min_top_feature_orig_1_composite_cropped_annotated.png")
min_top_feat_orig_2 <- file.path(sc_crop_dir, "min_top_feature_orig_2", "min_top_feature_orig_2_composite_cropped_annotated.png")
min_top_feat_orig_3 <- file.path(sc_crop_dir, "min_top_feature_orig_3", "min_top_feature_orig_3_composite_cropped_annotated.png")
min_top_feat_deriv_1 <- file.path(sc_crop_dir, "min_top_feature_deriv_1", "min_top_feature_deriv_1_composite_cropped_annotated.png")
min_top_feat_deriv_2 <- file.path(sc_crop_dir, "min_top_feature_deriv_2", "min_top_feature_deriv_2_composite_cropped_annotated.png")
min_top_feat_deriv_3 <- file.path(sc_crop_dir, "min_top_feature_deriv_3", "min_top_feature_deriv_3_composite_cropped_annotated.png")

# Second top coefficient feature
max_second_top_feat_orig_1 <- file.path(sc_crop_dir, "max_second_top_feature_orig_1", "max_second_top_feature_orig_1_composite_cropped_annotated.png")
max_second_top_feat_orig_2 <- file.path(sc_crop_dir, "max_second_top_feature_orig_2", "max_second_top_feature_orig_2_composite_cropped_annotated.png")
max_second_top_feat_orig_3 <- file.path(sc_crop_dir, "max_second_top_feature_orig_3", "max_second_top_feature_orig_3_composite_cropped_annotated.png")
max_second_top_feat_deriv_1 <- file.path(sc_crop_dir, "max_second_top_feature_deriv_1", "max_second_top_feature_deriv_1_composite_cropped_annotated.png")
max_second_top_feat_deriv_2 <- file.path(sc_crop_dir, "max_second_top_feature_deriv_2", "max_second_top_feature_deriv_2_composite_cropped_annotated.png")
max_second_top_feat_deriv_3 <- file.path(sc_crop_dir, "max_second_top_feature_deriv_3", "max_second_top_feature_deriv_3_composite_cropped_annotated.png")

min_second_top_feat_orig_1 <- file.path(sc_crop_dir, "min_second_top_feature_orig_1", "min_second_top_feature_orig_1_composite_cropped_annotated.png")
min_second_top_feat_orig_2 <- file.path(sc_crop_dir, "min_second_top_feature_orig_2", "min_second_top_feature_orig_2_composite_cropped_annotated.png")
min_second_top_feat_orig_3 <- file.path(sc_crop_dir, "min_second_top_feature_orig_3", "min_second_top_feature_orig_3_composite_cropped_annotated.png")
min_second_top_feat_deriv_1 <- file.path(sc_crop_dir, "min_second_top_feature_deriv_1", "min_second_top_feature_deriv_1_composite_cropped_annotated.png")
min_second_top_feat_deriv_2 <- file.path(sc_crop_dir, "min_second_top_feature_deriv_2", "min_second_top_feature_deriv_2_composite_cropped_annotated.png")
min_second_top_feat_deriv_3 <- file.path(sc_crop_dir, "min_second_top_feature_deriv_3", "min_second_top_feature_deriv_3_composite_cropped_annotated.png")

# Load top coefficient feat images
max_top_feat_orig_1_image <- load_image(max_top_feat_orig_1)
max_top_feat_orig_2_image <- load_image(max_top_feat_orig_2)
max_top_feat_orig_3_image <- load_image(max_top_feat_orig_3)
max_top_feat_deriv_1_image <- load_image(max_top_feat_deriv_1)
max_top_feat_deriv_2_image <- load_image(max_top_feat_deriv_2)
max_top_feat_deriv_3_image <- load_image(max_top_feat_deriv_3)

min_top_feat_orig_1_image <- load_image(min_top_feat_orig_1)
min_top_feat_orig_2_image <- load_image(min_top_feat_orig_2)
min_top_feat_orig_3_image <- load_image(min_top_feat_orig_3)
min_top_feat_deriv_1_image <- load_image(min_top_feat_deriv_1)
min_top_feat_deriv_2_image <- load_image(min_top_feat_deriv_2)
min_top_feat_deriv_3_image <- load_image(min_top_feat_deriv_3)

# Load second top coefficient feat images
max_second_top_feat_orig_1_image <- load_image(max_second_top_feat_orig_1)
max_second_top_feat_orig_2_image <- load_image(max_second_top_feat_orig_2)
max_second_top_feat_orig_3_image <- load_image(max_second_top_feat_orig_3)
max_second_top_feat_deriv_1_image <- load_image(max_second_top_feat_deriv_1)
max_second_top_feat_deriv_2_image <- load_image(max_second_top_feat_deriv_2)
max_second_top_feat_deriv_3_image <- load_image(max_second_top_feat_deriv_3)

min_second_top_feat_orig_1_image <- load_image(min_second_top_feat_orig_1)
min_second_top_feat_orig_2_image <- load_image(min_second_top_feat_orig_2)
min_second_top_feat_orig_3_image <- load_image(min_second_top_feat_orig_3)
min_second_top_feat_deriv_1_image <- load_image(min_second_top_feat_deriv_1)
min_second_top_feat_deriv_2_image <- load_image(min_second_top_feat_deriv_2)
min_second_top_feat_deriv_3_image <- load_image(min_second_top_feat_deriv_3)

# Create list of images
list_of_images <- list(
    max_top_feat_orig_1_image,
    max_top_feat_orig_2_image,
    max_top_feat_orig_3_image,
    max_top_feat_deriv_1_image,
    max_top_feat_deriv_2_image,
    max_top_feat_deriv_3_image,

    min_top_feat_orig_1_image,
    min_top_feat_orig_2_image,
    min_top_feat_orig_3_image,
    min_top_feat_deriv_1_image,
    min_top_feat_deriv_2_image,
    min_top_feat_deriv_3_image,

    max_second_top_feat_orig_1_image,
    max_second_top_feat_orig_2_image,
    max_second_top_feat_orig_3_image,
    max_second_top_feat_deriv_1_image,
    max_second_top_feat_deriv_2_image,
    max_second_top_feat_deriv_3_image,

    min_second_top_feat_orig_1_image,
    min_second_top_feat_orig_2_image,
    min_second_top_feat_orig_3_image,
    min_second_top_feat_deriv_1_image,
    min_second_top_feat_deriv_2_image,
    min_second_top_feat_deriv_3_image
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
## Top coefficient feature
max_top_feat_images_orig <- (
    list_of_images[[1]]
    + list_of_images[[2]]
    + list_of_images[[3]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
max_top_feat_images_deriv <- (
    list_of_images[[4]]
    + list_of_images[[5]]
    + list_of_images[[6]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_top_feat_images_orig <- (
    list_of_images[[7]]
    + list_of_images[[8]]
    + list_of_images[[9]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_top_feat_images_deriv <- (
    list_of_images[[10]]
    + list_of_images[[11]]
    + list_of_images[[12]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
## Second top coefficient feature
max_second_top_feat_images_orig <- (
    list_of_images[[13]]
    + list_of_images[[14]]
    + list_of_images[[15]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
max_second_top_feat_images_deriv <- (
    list_of_images[[16]]
    + list_of_images[[17]]
    + list_of_images[[18]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_second_top_feat_images_orig <- (
    list_of_images[[19]]
    + list_of_images[[20]]
    + list_of_images[[21]]
    + plot_layout(nrow = 1, widths = c(0.5, 0.5, 0.5))
)
min_second_top_feat_images_deriv <- (
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
    + geom_text(aes(x = 0.5, y = 0.5, label = "Second highest coefficient feature:\nDistribution of F-actin stain on edge of nuclei"), size = text_size) 
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
    wrap_elements(max_top_feat_images_orig) + 
    wrap_elements(max_top_feat_images_deriv) + 
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
        wrap_elements(min_top_feat_images_orig) + 
        wrap_elements(min_top_feat_images_deriv) + 
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

# Read image
top_feat_plot_min_img <- image_read("./nuclei_radial_feature_montage_min.png")

# Crop: width 9000, height 5628, starting at y = 972
top_feat_plot_min_img_cropped <- image_crop(top_feat_plot_min_img, "9000x5628+0+972")

# Save result
image_write(top_feat_plot_min_img_cropped, "./nuclei_radial_feature_montage_min.png")

# Create montage
width <- 15
height <- 11

options(repr.plot.width = width, repr.plot.height = height)

# patch the images together
second_top_feat_plot_min <- (
    blank + 
    (
        wrap_elements(full = second_top_patch_text) + 
        wrap_elements(min_second_top_feat_images_orig) + 
        wrap_elements(min_second_top_feat_images_deriv) + 
        plot_layout(ncol = 1, heights = c(0.35, 1, 1))
    ) +
    plot_layout(widths = c(0.1, 1)) # Adjusts the width of the left area for the text
)

second_top_feat_plot_min

# save plot
ggsave(
    file.path(
        paste0(
            "./","actin_radial_feature_montage_min.png"
        )
    ),
    second_top_feat_plot_min, width = width, height = height, dpi = 600
)

# Load the image
second_top_feat_plot_min_img <- image_read("./actin_radial_feature_montage_min.png")

# Crop to 8208px width and 6600px height, starting 792px from the left
second_top_feat_plot_min_img_cropped <- image_crop(second_top_feat_plot_min_img, "8208x6600+792+0")

# Save the cropped image
image_write(second_top_feat_plot_min_img_cropped, "./actin_radial_feature_montage_min.png")

# Create montage
width <- 15
height <- 11

options(repr.plot.width = width, repr.plot.height = height)

# patch the images together
second_top_feat_plot_max <- (
    blank + 
    (
        wrap_elements(full = blank) + 
        wrap_elements(max_second_top_feat_images_orig) + 
        wrap_elements(max_second_top_feat_images_deriv) + 
        plot_layout(ncol = 1, heights = c(0.35, 1, 1))
    ) +
    plot_layout(widths = c(0.1, 1)) # Adjusts the width of the left area for the text
)

second_top_feat_plot_max

# save plot
ggsave(
    file.path(
        paste0(
            "./","actin_radial_feature_montage_max.png"
        )
    ),
    second_top_feat_plot_max, width = width, height = height, dpi = 600
)

# Load the image
second_top_feat_plot_max_img <- image_read("./actin_radial_feature_montage_max.png")

# Crop: width = 8208, height = 5628, from (x=792, y=972)
second_top_feat_plot_max_img_cropped <- image_crop(second_top_feat_plot_max_img, "8208x5628+792+972")

# Save it
image_write(second_top_feat_plot_max_img_cropped, "./actin_radial_feature_montage_max.png")
