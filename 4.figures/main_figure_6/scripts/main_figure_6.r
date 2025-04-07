suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(RColorBrewer))

figure_dir <- "../figures"
output_main_figure_6 <- file.path(
    figure_dir, "main_figure_6_cell_line_generalizability.png"
)
results_dir <- file.path(
    "../../3.assess_generalizability/results"
)

# Load in platemap data
platemap_df <- read.csv(
    "/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data/0.download_data/metadata/platemap_NF1_plate6.csv"
)

# Path to UMAP results
UMAP_results_dir <- file.path(
    "/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data/4.analyze_data/notebooks/UMAP/results/qc_profiles_results"
)

# Subset the data frame and rename columns
platemap_df_filtered <- platemap_df[, c("well_position", "Institution")]
colnames(platemap_df_filtered) <- c("Metadata_Well", "Metadata_Institution")

# Load data
UMAP_results_file <- file.path(UMAP_results_dir, "UMAP_Plate_6_sc_only_model_features_qc.tsv")

UMAP_results_df <- readr::read_tsv(UMAP_results_file)

# Merge institution info onto UMAP df
UMAP_results_df <- platemap_df_filtered %>% inner_join(UMAP_results_df, by = "Metadata_Well")

# Add new column for cell line derivative
UMAP_results_df <- UMAP_results_df %>%
    mutate(cell_line_derivative = ifelse(Metadata_Institution == "iNFixion", "original", 
        ifelse(Metadata_Institution == "MGH", "derivative", NA)))

# Add new column for cell line code
UMAP_results_df <- UMAP_results_df %>%
    mutate(cell_line_code = case_when(
        Metadata_Institution == "iNFixion" & Metadata_genotype == "Null" ~ "C04",
        Metadata_Institution == "iNFixion" & Metadata_genotype == "WT" ~ "A3",
        Metadata_Institution == "MGH" & Metadata_genotype == "WT" ~ "GFP 3",
        Metadata_Institution == "MGH" & Metadata_genotype == "Null" ~ "C23",
        TRUE ~ NA_character_
    ))

dim(UMAP_results_df)
head(UMAP_results_df)

umap_fig_gg <- (
    ggplot(UMAP_results_df, aes(x = UMAP0, y = UMAP1)) +
    geom_point(
        aes(color = Metadata_genotype),
        size = 1.0,
        alpha = 0.4
    ) +
    theme_bw() +
    guides(
        color = guide_legend(
            override.aes = list(size = 5)
        )
    ) +
    labs(x = "UMAP0", y = "UMAP1", color = "NF1\ngenotype") +
    facet_wrap(
        cell_line_code ~ cell_line_derivative, 
        labeller = labeller(
            cell_line_code = function(x) paste("Cell line code:", x),
            cell_line_derivative = function(x) paste("ipn02.3 2λ:", x)
        )
    ) +
    coord_fixed(ratio = 1.0) + 
    # Change the text size
    theme(
        strip.text = element_text(size = 17),
        # X and Y axis text size
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        # X and Y axis title size
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        # Legend text size
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22)
    )
)

umap_fig_gg

# Load data (includes optimization in this file)
PR_results_file <- file.path(results_dir, "plate6_precision_recall_final_model_qc.parquet")

PR_results_df <- arrow::read_parquet(PR_results_file)

# Create new column for model using the datasplit prefix
PR_results_df <- PR_results_df %>%
  mutate(shuffled_type = ifelse(grepl("^shuffled", data_type), "TRUE", "FALSE"))

# Add new column for cell line derivative
PR_results_df <- PR_results_df %>%
    mutate(cell_line_derivative = ifelse(Metadata_Institution == "iNFixion", "original", ifelse(Metadata_Institution == "MGH", "derivative", NA)))

dim(PR_results_df)
head(PR_results_df)

width <- 12
height <- 12
options(repr.plot.width = width, repr.plot.height = height)

pr_curve_plot <- (
    ggplot(PR_results_df, aes(x = Recall, y = Precision, color = cell_line_derivative, linetype = shuffled_type))
    + geom_line(aes(linetype = shuffled_type), linewidth = 1)
    + theme_bw()
    # + coord_fixed()
    + labs(color = "ipn02.3 2λ", linetype = "Features\nshuffled", x = "Recall", y = "Precision")
    # change the colors
    + scale_color_manual(values = c(
        "original" = brewer.pal(8, "Dark2")[4],
        "derivative" = brewer.pal(8, "Dark2")[3]
    ))
    + scale_y_continuous(limits = c(0, 1))
    # change the line thickness of the lines in the legend
    + guides(linetype = guide_legend(override.aes = list(size = 1)))
    # change the text size
    + theme(
        # x and y axis text size
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        # x and y axis title size
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        # legend text size
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
    )
)

pr_curve_plot

# Load data
accuracy_results_file <- file.path(results_dir, "plate6_accuracy_final_model_qc.parquet")

accuracy_results_df <- arrow::read_parquet(accuracy_results_file)

# Create new column for model using the datasplit prefix
accuracy_results_df <- accuracy_results_df %>%
  mutate(shuffled_type = ifelse(grepl("^shuffled", data_type), "TRUE", "FALSE"))

# Add new column for cell line derivative
accuracy_results_df <- accuracy_results_df %>%
    mutate(cell_line_derivative = ifelse(Metadata_Institution == "iNFixion", "original", ifelse(Metadata_Institution == "MGH", "derivative", NA)))

dim(accuracy_results_df)
head(accuracy_results_df)

width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)
# bar plot of the accuracy scores
accuracy_score_plot <- (
    ggplot(accuracy_results_df, aes(x = shuffled_type, y = accuracy, fill = cell_line_derivative))
    + geom_bar(stat = "identity", position = "dodge")

    # Add text labels for accuracy scores on top of bars
    + geom_text(
        aes(label = sprintf("%.2f", accuracy)), 
        position = position_dodge(width = 0.9), 
        vjust = -0.5, 
        size = 6
    )

    + ylim(0, 1)
    + theme_bw()
    + ylab("Accuracy")
    + xlab("Features shuffled")
    # change the legend title
    + labs(fill = "ipn02.3 2λ")
    # change the colours
    + scale_fill_manual(values = c(
        "original" = brewer.pal(8, "Dark2")[4],
        "derivative" = brewer.pal(8, "Dark2")[3]
    ))
    # change the text size
    + theme(
        # x and y axis text size
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        # x and y axis title size
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        # legend text size
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
    )
)

accuracy_score_plot

# Load data
kstest_results_file <- file.path("../../3.assess_generalizability/results/ks_test_derivatives_results_qc.parquet")

kstest_results_df <- arrow::read_parquet(kstest_results_file)

# Filter for genotype_comparison == "All"
kstest_results_df <- kstest_results_df %>%
    filter(genotype_comparison == "All")

# Create a new column extracting the first part of 'feature' after the compartment
kstest_results_df$feature_base <- sub("^[^_]+_", "", kstest_results_df$feature)

# Update the channel column where anything other than DAPI, CY5, GFP, or RFP is called "other"
kstest_results_df$channel <- ifelse(kstest_results_df$channel %in% c("DAPI", "CY5", "GFP", "RFP"), kstest_results_df$channel, "Other")

# Update the channel names
kstest_results_df$channel <- recode(kstest_results_df$channel, "DAPI" = "Nuclei", "GFP" = "ER", "CY5" = "Mito", "RFP" = "F-actin")

dim(kstest_results_df)
head(kstest_results_df)

width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

# Update feature group name
kstest_results_df$feature_group <- ifelse(kstest_results_df$feature_group == "RadialDistribution", "RadialDist", kstest_results_df$feature_group)

# Reorder feature_base based on the sum of the base feature across all compartments
kstest_results_df$feature_base <- factor(kstest_results_df$feature_base, 
                                         levels = kstest_results_df %>% 
                                           group_by(feature_base) %>% 
                                           summarise(total_ks = sum(ks_stat)) %>% 
                                           arrange(desc(total_ks)) %>% 
                                           pull(feature_base))

# Create the plot
ks_test_scatter <- (
    ggplot(kstest_results_df, aes(x = feature_base, y = ks_stat))
    + geom_point(aes(color = feature_group, size = feature_importances, shape = channel), alpha = 0.4) 
    + theme_bw()
    + facet_grid(compartment ~ .)
    + theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 23),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        strip.text = element_text(size = 22),
        panel.spacing = unit(1, "lines"),
        legend.position = "right"
    )
    + ylim(0,1)
    + scale_color_discrete(name = "Feature\ngroup")
    + scale_size_continuous(name = "Feature\nimportance", range = c(1, 8)) 
    + scale_shape_manual(name = "Organelle", values = c(16, 17, 15, 18, 7))
    + labs(
        x = "Rank of summed KS statistic per core feature group\n(summed across compartments)",
        y = "KS test statistic"
    )
    + guides(
    shape = guide_legend(override.aes = list(size = 5)), 
    color = guide_legend(override.aes = list(size = 5))
    )
)

ks_test_scatter

# Drop rows where genotype is 'HET'
platemap_df <- platemap_df %>%
    filter(genotype != "HET")

# Add new column for cell line code
platemap_df <- platemap_df %>%
    mutate(cell_line_code = case_when(
        Institution == "iNFixion" & genotype == "Null" ~ "C04",
        Institution == "iNFixion" & genotype == "WT" ~ "A3",
        Institution == "MGH" & genotype == "WT" ~ "GFP 3",
        Institution == "MGH" & genotype == "Null" ~ "C23",
        TRUE ~ NA_character_
    ))

head(platemap_df)

width <- 10
height <- 8
options(repr.plot.width = width, repr.plot.height = height)

# Generate platemap figure
platemap <-
    platetools::raw_map(
        data = platemap_df$genotype,
        well = platemap_df$well_position,
        plate = 96,
        size = 12
    ) +
    coord_fixed(ratio = 1.0) +
    ggplot2::scale_fill_discrete(name = "NF1\ngenotype") +
    ggplot2::geom_point(aes(shape = platemap_df$cell_line_code), size= 3) +
    ggplot2::scale_shape_discrete(name = "Cell line\ncode") +
    theme(
        legend.title = element_text(size = 22),  # Larger legend title
        legend.text = element_text(size = 20),  # Larger legend text
        axis.text = element_text(size = 22),  # Larger axis tick labels
        axis.title = element_text(size = 22)  # Larger axis titles
    ) +
    guides(
        shape = guide_legend(override.aes = list(size = 5)), 
    )

platemap

height <- 8
width <- 6
options(repr.plot.width = width, repr.plot.height = height)

right_plot <- (
    pr_curve_plot /
    accuracy_score_plot
) + plot_layout(heights = c(2,2))

right_plot

height <- 8
width <- 15
options(repr.plot.width = width, repr.plot.height = height)

align_bottom_plot <- (
    right_plot |
    ks_test_scatter
) + plot_layout(widths = c(1,1))

align_bottom_plot

height <- 10
width <- 16
options(repr.plot.width = width, repr.plot.height = height)

platemap_umap <- (
    platemap /
    umap_fig_gg
) + plot_layout(heights = c(1,1.5))

platemap_umap

height <- 12
width <- 25
options(repr.plot.width = width, repr.plot.height = height)

align_plot <- (
    platemap_umap |
    align_bottom_plot
) + plot_layout(widths = c(0.5,1))

align_plot

fig_6_gg <- (
  align_plot
) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 25))

# Save the plot
ggsave(output_main_figure_6, plot = fig_6_gg, dpi = 500, height = height, width = width)
