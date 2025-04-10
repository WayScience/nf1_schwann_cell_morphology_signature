suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(RColorBrewer))

figure_dir <- "../figures/supplementary"
output_supp_figure <- file.path(figure_dir, "supp_figure_5_kstestgenotype.png")

# Load data
kstest_results_file <- file.path("../../3.assess_generalizability/results/ks_test_derivatives_results_qc.parquet")

kstest_results_df <- arrow::read_parquet(kstest_results_file)

# Filter only ks test results from each genotype
kstest_results_df <- kstest_results_df %>%
    filter(genotype_comparison != "All")

# Create a new column extracting the first part of 'feature' after the compartment
kstest_results_df$feature_base <- sub("^[^_]+_", "", kstest_results_df$feature)

# Update the channel column where anything other than DAPI, CY5, GFP, or RFP is called "other"
kstest_results_df$channel <- ifelse(kstest_results_df$channel %in% c("DAPI", "CY5", "GFP", "RFP"), kstest_results_df$channel, "Other")

# Update the channel names
kstest_results_df$channel <- recode(kstest_results_df$channel, "DAPI" = "Nuclei", "GFP" = "ER", "CY5" = "Mito", "RFP" = "F-actin")

dim(kstest_results_df)
head(kstest_results_df)

width <- 12
height <- 10
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
supp_fig_gg <- (
    ggplot(kstest_results_df, aes(x = feature_base, y = ks_stat))
    + geom_point(aes(color = feature_group, size = feature_importances, shape = channel), alpha = 0.4) 
    + theme_bw()
    + facet_grid(compartment ~ genotype_comparison)
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

# Save or display the plot
ggsave(output_supp_figure, plot = supp_fig_gg, dpi = 500, height = height, width = width)

supp_fig_gg
