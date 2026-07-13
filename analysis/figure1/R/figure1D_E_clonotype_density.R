# Figure 1D/E: unique clonotype counts per sample and clonotype density per mm2.

source("analysis/figure1/R/00_setup.R")

spot_diameter <- 55
spot_area_mm2 <- pi * (spot_diameter / 2)^2 / 1e6

calc_density <- function(obj, assay_name) {
  n_features <- nrow(obj@assays[[assay_name]])
  n_spots <- ncol(obj@assays[[assay_name]])
  n_features / (n_spots * spot_area_mm2)
}

results <- lapply(seq_len(nrow(sample_config)), function(i) {
  object_name <- sample_config$object_name[i]
  obj <- get(object_name)
  data.frame(
    Sample = sample_config$display_sample[i],
    TRB = calc_density(obj, "aaSeqCDR3_Lv1_UMI3_TRB_20G"),
    IGH = calc_density(obj, "aaSeqCDR3_Lv1_UMI3_IGH_20G")
  )
}) |> bind_rows()

write.csv(results, file.path(output_dir, "figure1D_unique_clonotypes_per_sample.csv"), row.names = FALSE)

df_melt <- reshape2::melt(results, id.vars = "Sample", variable.name = "Receptor", value.name = "Density")
df_melt$Sample <- factor(df_melt$Sample, levels = sample_levels)

p <- ggplot(df_melt, aes(x = Receptor, y = Density)) +
  stat_boxplot(geom = "errorbar", width = 0.2, color = "black") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, fatten = 2, color = "black") +
  stat_summary(fun = mean, geom = "bar", width = 0.6, fill = NA, color = "black") +
  geom_point(aes(color = Sample), size = 3, position = position_jitter(width = 0.15)) +
  scale_color_manual(values = sample_colors, drop = FALSE) +
  labs(y = expression("Unique clonotypes per mm"^2), x = "", title = "Clonotype density per mm2") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

save_panel(p, "figure1E_clonotype_density_per_mm2.pdf", width = 4, height = 4)
