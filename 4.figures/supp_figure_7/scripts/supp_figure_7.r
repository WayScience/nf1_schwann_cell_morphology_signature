suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))

figure_dir <- "../figures/supplementary"
output_supp_figure <- file.path(figure_dir, "supp_figure_7_qualitycontrol.png")

# Base directory containing Plate_ folders
plates_dir <- "/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data/1.cellprofiler_ic/image_quality_control/"

# List all Image.csv files in Plate_ folders except Plate_4
image_csv_paths <- list.files(
  path = plates_dir,
  pattern = "Image\\.csv$",
  recursive = TRUE,
  full.names = TRUE
) %>%
  # Keep only paths that contain "Plate_" in the directory and exclude Plate_4
  keep(~ str_detect(.x, "Plate_") && !str_detect(.x, "Plate_4"))

# Read all Image.csv files and bind into one dataframe
image_data_df <- image_csv_paths %>%
  set_names() %>%
  map_dfr(read_csv, .id = "source_file", show_col_types = FALSE)

# Add Metadata_Plate column (capture Plate_3_prime and similar names)
image_data_df <- image_data_df %>%
  mutate(
    Metadata_Plate = str_extract(source_file, "Plate_\\d+(_[A-Za-z0-9]+)?"),
    Metadata_Plate = ifelse(is.na(Metadata_Plate), "Unknown", Metadata_Plate)
  )

# Check result
dim(image_data_df)
head(image_data_df)

# Rename plates in Metadata_Plate
image_data_df <- image_data_df %>%
    mutate(
        Metadata_Plate = recode(
            Metadata_Plate,
            "Plate_3" = "Plate A",
            "Plate_3_prime" = "Plate B",
            "Plate_5" = "Plate C",
            "Plate_6" = "Plate D"
        )
    )

# Check result
unique(image_data_df$Metadata_Plate)

# Add Metadata_Well and Metadata_Site columns
image_data_df <- image_data_df %>%
  mutate(
    Metadata_Well = str_split_fixed(FileName_OrigRFP, "_", n = 6)[, 1],
    Metadata_Site = str_split_fixed(FileName_OrigRFP, "_", n = 6)[, 4]
  )

dim(image_data_df)
head(image_data_df)

qc_results_df <- image_data_df %>%
    select(
        starts_with("Metadata_") 
        | contains("PowerLogLogSlope")
        | contains("PercentMaximal")
    )
dim(qc_results_df)
head(qc_results_df)

channels <- c("OrigDAPI", "OrigGFP", "OrigRFP", "OrigCY5")

qc_df <- bind_rows(lapply(channels, function(channel) {
  qc_results_df %>%
    select(
      starts_with("Metadata_"),
      all_of(paste0("ImageQuality_PowerLogLogSlope_", channel)),
      all_of(paste0("ImageQuality_PercentMaximal_", channel))
    ) %>%
    rename(
      PowerLogLogSlope = paste0("ImageQuality_PowerLogLogSlope_", channel),
      PercentMaximal = paste0("ImageQuality_PercentMaximal_", channel)
    ) %>%
    mutate(Channel = channel)
}))

# Remove Orig prefix from Channel column
qc_df <- qc_df %>%
  mutate(
    Channel = str_replace(Channel, "^Orig", ""),  # Remove 'Orig' prefix from values
    Channel = factor(Channel, levels = c("DAPI", "GFP", "CY5", "RFP"))  # Set factor levels
  )

dim(qc_df)
head(qc_df)

width <- 18
height <- 6
options(repr.plot.width = width, repr.plot.height = height) 

# Threshold values
upper_threshold <- -1.36839084751248
lower_threshold <- -2.4936740062477365

# Create the plot with the specified color mapping
qc_blur_plot <- ggplot(qc_df, aes(x = PowerLogLogSlope, fill = Channel)) +
    geom_density(alpha = 0.5, adjust = 1.5) + 
    scale_fill_manual(
        values = c("DAPI" = "#4472C4",  # Blue
                   "GFP" = "#70AD47",  # Green
                   "CY5" = "#D35B9D",  # Brighter magenta
                   "RFP" = "#b42718")  # Red
    ) +
    geom_vline(aes(xintercept = upper_threshold, color = paste("Upper =", round(upper_threshold, 2))), 
               linetype = "dashed", linewidth = 1.2) + 
    geom_vline(aes(xintercept = lower_threshold, color = paste("Lower =", round(lower_threshold, 2))), 
               linetype = "dashed", linewidth = 1.2) +  
    scale_color_manual(name = "QC thresholds", values = c("Upper = -1.37" = "red", 
                                                          "Lower = -2.49" = "red")) +  # Set the color and label for the legend
    labs(
        x = "PowerLogLogSlope",
        y = "Density"
    ) +
    theme_bw() +
    facet_wrap(~ Metadata_Plate) +  # Facet by plate
    theme(
        strip.text = element_text(size = 16),   # Facet titles
        axis.title = element_text(size = 16),  # Axis titles
        axis.text = element_text(size = 14),   # Axis text
        legend.title = element_text(size = 16), # Legend title
        legend.text = element_text(size = 14)   # Legend text
    )

qc_blur_plot


width <- 15
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

# Take natural logarithm plus one of the percent maximal for better visualization
qc_df <- qc_df %>%
  mutate(LogImageQuality = log1p(PercentMaximal))

# Create a label for the threshold line
percentmaximal_threshold <- 1.4

# Create the plot
qc_saturation_plot <- ggplot(qc_df, aes(x = LogImageQuality)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black") +
  facet_wrap(~ Metadata_Plate) +
  geom_vline(aes(xintercept = log1p(percentmaximal_threshold), 
                 color = paste("Threshold =", percentmaximal_threshold, "%")),
             linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "QC thresholds", values = "red") + # Customize the legend title and color
  labs(
    x = "Natural logarithm of (percent maximal percentage + 1)",
    y = "Count"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16),   # Facet titles
    axis.title = element_text(size = 16),  # Axis titles
    axis.text = element_text(size = 14),   # Axis text
    legend.title = element_text(size = 16), # Legend title
    legend.text = element_text(size = 14)   # Legend text
  )

qc_saturation_plot

align_plot <- (
    qc_blur_plot /
    qc_saturation_plot 
) + plot_layout(heights= c(2,2))

align_plot

supp_fig_gg <- (
  align_plot
) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 25))

# Save or display the plot
ggsave(output_supp_figure, plot = supp_fig_gg, dpi = 500, height = 8, width = 10)

supp_fig_gg

# Define the directory containing Corrected_Images folders
corrected_dir <- file.path("/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data/1.cellprofiler_ic/Corrected_Images")

# List all RunImage.csv files for Plates 3, 3_prime, 5, and 6
runimage_csv_paths <- list.files(
    path = corrected_dir,
    pattern = "RunImage\\.csv$",
    recursive = TRUE,
    full.names = TRUE
) %>%
    keep(~ str_detect(.x, "Plate_3($|/)")
                | str_detect(.x, "Plate_3_prime($|/)")
                | str_detect(.x, "Plate_5($|/)")
                | str_detect(.x, "Plate_6($|/)"))

# Read all RunImage.csv files and bind into one dataframe
runimage_data_df <- runimage_csv_paths %>%
    set_names() %>%
    map_dfr(read_csv, .id = "source_file", show_col_types = FALSE)

# Check result
dim(runimage_data_df)
head(runimage_data_df)

# Count number of FOVs that failed QC
num_failed <- runimage_data_df %>% filter(Image_Quality_Control_QC_Flag == 1) %>% nrow()

# Total number of FOVs
num_total <- nrow(runimage_data_df)

# Percentage failed
percent_failed <- (num_failed / num_total) * 100

cat("Number of FOVs failed QC:", num_failed, "\n")
cat("Total number of FOVs:", num_total, "\n")
cat(sprintf("Percentage failed: %.2f%%\n", percent_failed))
