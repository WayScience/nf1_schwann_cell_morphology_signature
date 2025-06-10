suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(RColorBrewer))

figure_dir <- "../figures/supplementary"
output_supp_figure_7 <- file.path(
    figure_dir, "supp_figure_7_log_reg_pr_curves.png"
)
results_dir <- file.path(
    "./pr_results"
)

# Load data
PR_results_file <- file.path(results_dir, "pr_curve_all_plates.parquet")

PR_results_df <- arrow::read_parquet(PR_results_file)

dim(PR_results_df)
head(PR_results_df)

width <- 12
height <- 12
options(repr.plot.width = width, repr.plot.height = height)
pr_train_test_plot <- (
    ggplot(PR_results_df, aes(x = recall, y = precision, color = split, linetype = model))
    + geom_line(aes(linetype = model), linewidth = 1)
    + theme_bw()
    + coord_fixed()
    + labs(color = "ML model\ndata split", linetype = "Features\nshuffled", x = "Recall", y = "Precision")
    # change the colors
    + scale_color_manual(values = c(
        "test" = brewer.pal(8, "Dark2")[6],
        "train" = brewer.pal(8, "Dark2")[3],
        "holdout" = brewer.pal(8, "Dark2")[1]
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

pr_train_test_plot

# Load data
PR_results_file <- file.path(results_dir, "pr_curve_per_plate.parquet")

PR_results_df <- arrow::read_parquet(PR_results_file)

# Update plate names
PR_results_df <- PR_results_df %>%
    mutate(
        plate = recode(
            plate,
            "Plate_3" = "Plate_A",
            "Plate_3_prime" = "Plate_B",
            "Plate_5" = "Plate_C",
            "Plate_6" = "Plate_D"
        )
    )

dim(PR_results_df)
head(PR_results_df)

width <- 12
height <- 12
options(repr.plot.width = width, repr.plot.height = height)
pr_per_plate_plot <- (
    ggplot(PR_results_df, aes(x = recall, y = precision, color = split, linetype = model))
    + geom_line(aes(linetype = model), linewidth = 1)
    + theme_bw()
    + coord_fixed()
    + facet_wrap(~ plate, ncol = 2)
    + labs(color = "ML model\ndata split", linetype = "Features\nshuffled", x = "Recall", y = "Precision")
    # change the colors
    + scale_color_manual(values = c(
        "test" = brewer.pal(8, "Dark2")[6],
        "train" = brewer.pal(8, "Dark2")[3],
        "holdout" = brewer.pal(8, "Dark2")[1]
    ))
    + scale_y_continuous(limits = c(0, 1))
    # change the line thickness of the lines in the legend
    + guides(linetype = guide_legend(override.aes = list(size = 1)))  
    # change the text size
    + theme(
        # x and y axis text size
        axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 22),
        # x and y axis title size
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        # legend text size
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        # strip text size for facets
        strip.text = element_text(size = 20)
    )
)

pr_per_plate_plot

# Load data
PR_results_file <- file.path(results_dir, "pr_curve_plate6_per_institution.parquet")

PR_results_df <- arrow::read_parquet(PR_results_file)

dim(PR_results_df)
head(PR_results_df)

width <- 12
height <- 12
options(repr.plot.width = width, repr.plot.height = height)
pr_plate6_derivative_plot <- (
    ggplot(PR_results_df, aes(x = recall, y = precision, color = institution, linetype = model))
    + geom_line(aes(linetype = model), linewidth = 1)
    + theme_bw()
    + coord_fixed()
    + labs(color = "Cell line\nderivative", linetype = "Features\nshuffled", x = "Recall", y = "Precision")
    # change the colors
    + scale_color_manual(values = c(
        "MGH" = brewer.pal(8, "Dark2")[7],
        "iNFixion" = brewer.pal(8, "Dark2")[8]
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

pr_plate6_derivative_plot

width <- 22
height <- 6
options(repr.plot.width = width, repr.plot.height = height)

align_plot <- 
    pr_train_test_plot + 
    pr_per_plate_plot + 
    pr_plate6_derivative_plot +
    plot_layout(widths = c(1,1,1))

align_plot

supp_fig_7_gg <- (
  align_plot + plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(plot.tag = element_text(size = 30))
)

# Save or display the plot
ggsave(output_supp_figure_7, plot = supp_fig_7_gg, dpi = 500, height = height, width = width)
