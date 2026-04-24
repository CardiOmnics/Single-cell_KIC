# =============================================================================
# cluster_helpers.R
# Helper functions for subpopulation (k-means) clustering analysis.
# Sourced by 4_plot_and_stats.R and figure_calcium_clusters.R.
# Requires: tidyverse, ggplot2, theme_project(), group_colors (from setup files)
# =============================================================================

# Run k-means subpopulation analysis, produce individual diagnostic plots,
# and return the ggplot objects so figure assembly scripts can reuse them.
run_subpop_analysis <- function(data, param_cols, param_labels,
                                session_label, fig_dir, analysis_date,
                                k = 3L) {

  # --- Prepare matrix: complete cases, scaled ---
  clust_data <- data %>%
    select(cell_uid, well_uid, Condition, Cell_Line, all_of(param_cols)) %>%
    drop_na()

  if (nrow(clust_data) < 50) {
    cat("    Skipped", session_label, "clustering (too few complete cases)\n")
    return(NULL)
  }

  clust_mat <- clust_data %>% select(all_of(param_cols)) %>% scale()

  # --- Fixed k (biologically motivated: healthy / mildly stressed / impaired) ---
  optimal_k <- as.integer(k)
  cat("  Using fixed k =", optimal_k, "for", session_label, "\n")

  # --- k-means ---
  set.seed(42)
  km_fit <- kmeans(clust_mat, centers = optimal_k, nstart = 50, iter.max = 200)
  clust_data$cluster <- factor(km_fit$cluster)

  # Print cluster size summary so user can evaluate balance
  clust_sizes <- sort(table(km_fit$cluster), decreasing = TRUE)
  cat("  Cluster sizes:", paste(names(clust_sizes), clust_sizes, sep = "=", collapse = ", "), "\n")
  cat("  Largest/smallest ratio:", round(max(clust_sizes) / min(clust_sizes), 1), "\n")

  # --- Characterize clusters from k-means centers (guarantees optimal_k rows) ---

  # Z-scored centers: used for heatmap
  profile_z <- as.data.frame(km_fit$centers)
  colnames(profile_z) <- param_cols
  profile_z$cluster <- factor(1:optimal_k)

  # Back-transform centers to original units: used for Excel export
  scale_center <- attr(clust_mat, "scaled:center")
  scale_scale  <- attr(clust_mat, "scaled:scale")
  centers_orig <- sweep(
    sweep(km_fit$centers, 2, scale_scale, "*"),
    2, scale_center, "+"
  )
  cluster_profiles <- as.data.frame(centers_orig)
  colnames(cluster_profiles) <- param_cols
  cluster_profiles$cluster <- factor(1:optimal_k)
  cluster_profiles$n_cells <- as.integer(table(km_fit$cluster))

  # --- Label clusters by dominant characteristic ---
  # Use the parameter with the most extreme mean z-score per cluster
  profile_mat <- profile_z %>% select(all_of(param_cols)) %>% as.matrix()
  rownames(profile_mat) <- paste0("C", profile_z$cluster)

  # Assign descriptive labels based on profiles
  cluster_sizes <- table(clust_data$cluster)
  clust_labels <- paste0("Cluster ", seq_len(optimal_k),
                         " (n=", cluster_sizes[seq_len(optimal_k)], ")")
  clust_data$cluster_label <- clust_labels[as.integer(clust_data$cluster)]

  # --- Plot 1: Cluster profile heatmap ---
  profile_long <- profile_z %>%
    pivot_longer(cols = all_of(param_cols),
                 names_to = "parameter", values_to = "z_score") %>%
    mutate(
      parameter = factor(parameter, levels = param_cols,
                         labels = param_labels),
      cluster = factor(cluster)
    )

  cluster_palette <- c("#2166AC", "#D95F02", "#1B9E77", "#E7298A", "#66A61E", "#7570B3")

  p_profile <- ggplot(profile_long,
                      aes(x = parameter, y = cluster, fill = z_score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", z_score)), size = 2.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         limits = c(-2, 2), oob = scales::squish,
                         name = "Mean z-score") +
    theme_project() +
    labs(title = paste0("Cluster profiles: ", session_label),
         subtitle = paste0("k = ", optimal_k, " clusters, n = ",
                           nrow(clust_data), " cells"),
         x = NULL, y = "Cluster") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.key.height = unit(3, "mm"),
          legend.key.width  = unit(8, "mm"),
          legend.title      = element_text(size = 7),
          legend.text       = element_text(size = 6))
  ggsave(file.path(fig_dir,
    paste0("fig_cluster_profiles_", tolower(gsub(" ", "_", session_label)),
           "_", analysis_date, ".pdf")),
    p_profile, width = 10, height = 4 + optimal_k * 0.5, device = cairo_pdf)
  cat("  Saved: cluster_profiles\n")

  # --- Plot 2: PCA colored by cluster ---
  pca_fit <- prcomp(clust_mat, center = FALSE, scale. = FALSE)
  var_exp <- summary(pca_fit)$importance[2, 1:2] * 100

  pca_scores <- data.frame(
    PC1 = pca_fit$x[, 1], PC2 = pca_fit$x[, 2],
    cluster = clust_data$cluster,
    Condition = clust_data$Condition,
    Cell_Line = clust_data$Cell_Line
  )

  # Trim PCA axis limits to 1st–99th percentile to avoid outlier-driven whitespace
  pc1_lo <- quantile(pca_scores$PC1, 0.01); pc1_hi <- quantile(pca_scores$PC1, 0.99)
  pc2_lo <- quantile(pca_scores$PC2, 0.01); pc2_hi <- quantile(pca_scores$PC2, 0.99)
  pc_pad <- 0.08
  pc1_range <- pc1_hi - pc1_lo; pc2_range <- pc2_hi - pc2_lo

  p_pca_clust <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(alpha = 0.15, size = 0.4) +
    stat_ellipse(level = 0.68, linewidth = 0.8) +
    coord_cartesian(
      xlim = c(pc1_lo - pc_pad * pc1_range, pc1_hi + pc_pad * pc1_range),
      ylim = c(pc2_lo - pc_pad * pc2_range, pc2_hi + pc_pad * pc2_range)
    ) +
    facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
    scale_color_manual(values = cluster_palette[1:optimal_k],
                       name = "Cluster") +
    theme_project() +
    labs(title = paste0("Functional subpopulations: ", session_label),
         subtitle = "k-means clusters projected onto PCA space",
         x = sprintf("PC1 (%.1f%% var)", var_exp[1]),
         y = sprintf("PC2 (%.1f%% var)", var_exp[2])) +
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir,
    paste0("fig_pca_clusters_", tolower(gsub(" ", "_", session_label)),
           "_", analysis_date, ".pdf")),
    p_pca_clust, width = 12, height = 6, device = cairo_pdf)
  cat("  Saved: pca_clusters\n")

  # --- Plot 3: Cluster composition per Condition (THE KEY FIGURE) ---
  composition <- clust_data %>%
    count(Condition, Cell_Line, cluster) %>%
    group_by(Condition, Cell_Line) %>%
    mutate(
      total = sum(n),
      pct = 100 * n / total
    ) %>%
    ungroup()

  p_comp <- ggplot(composition,
                   aes(x = Condition, y = pct, fill = cluster)) +
    geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(pct >= 5, sprintf("%.0f%%", pct), "")),
              position = position_stack(vjust = 0.5), size = 2, color = "white",
              fontface = "bold") +
    facet_wrap(~ Cell_Line, nrow = 1, labeller = cell_line_fig_labeller) +
    scale_fill_manual(values = cluster_palette[1:optimal_k],
                      name = "Cluster") +
    theme_project() +
    labs(title = paste0("Subpopulation shifts across conditions: ", session_label),
         subtitle = "Does lipid overload change the proportion of functional subpopulations?",
         x = NULL, y = "% of cells") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(fig_dir,
    paste0("fig_cluster_composition_", tolower(gsub(" ", "_", session_label)),
           "_", analysis_date, ".pdf")),
    p_comp, width = 10, height = 6, device = cairo_pdf)
  cat("  Saved: cluster_composition\n")

  # --- Plot 4: PCA colored by Condition, split by cluster ---
  p_pca_cond_clust <- ggplot(pca_scores,
                             aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(alpha = 0.15, size = 0.4) +
    stat_ellipse(level = 0.68, linewidth = 0.8) +
    facet_grid(cluster ~ Cell_Line, labeller = labeller(Cell_Line = cell_line_fig_labeller)) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    theme_project() +
    labs(title = paste0("Condition overlap within clusters: ", session_label),
         subtitle = "Clusters may be condition-enriched but not condition-exclusive",
         x = sprintf("PC1 (%.1f%% var)", var_exp[1]),
         y = sprintf("PC2 (%.1f%% var)", var_exp[2])) +
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir,
    paste0("fig_pca_condition_by_cluster_", tolower(gsub(" ", "_", session_label)),
           "_", analysis_date, ".pdf")),
    p_pca_cond_clust, width = 12, height = 4 + optimal_k * 2.5,
    device = cairo_pdf)
  cat("  Saved: pca_condition_by_cluster\n")

  # --- Well-level cluster proportion test (correct unit of analysis) ---
  # Compute % of cells in each cluster per well, then test across conditions.
  # This is the statistically defensible version: well = unit of analysis.
  cat("\n  Well-level cluster proportion test (Kruskal-Wallis):\n")

  well_clust_pct <- clust_data %>%
    count(well_uid, Condition, Cell_Line, cluster, .drop = FALSE) %>%
    group_by(well_uid, Condition, Cell_Line) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  well_kw <- well_clust_pct %>%
    group_by(Cell_Line, cluster) %>%
    group_map(~ {
      if (length(unique(.x$Condition)) < 2) return(NULL)
      kt <- kruskal.test(pct ~ Condition, data = .x)
      tibble(
        Cell_Line = .y$Cell_Line,
        cluster   = .y$cluster,
        KW_chi_sq = round(kt$statistic, 2),
        p_value   = round(kt$p.value, 4),
        sig       = case_when(
          kt$p.value < 0.001 ~ "***",
          kt$p.value < 0.01  ~ "**",
          kt$p.value < 0.05  ~ "*",
          TRUE               ~ "ns"
        ),
        n_wells = n_distinct(.x$well_uid)
      )
    }, .keep = TRUE) %>%
    bind_rows()

  print(as.data.frame(well_kw), row.names = FALSE)

  # Box + jitter plot: well-level % per cluster x condition
  p_well_pct <- ggplot(
    well_clust_pct,
    aes(x = Condition, y = pct, fill = Condition, color = Condition)
  ) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, linewidth = 0.5) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    facet_grid(cluster ~ Cell_Line,
               labeller = labeller(cluster   = function(x) paste0("Cluster ", x),
                                   Cell_Line = cell_line_fig_labeller)) +
    scale_fill_manual(values  = group_colors, drop = FALSE) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    labs(
      title    = paste0("Well-level cluster proportions: ", session_label),
      subtitle = "Each point = one well; Kruskal-Wallis p-values in console",
      x = NULL, y = "% of cells in cluster"
    ) +
    theme_project() +
    theme(legend.position = "bottom")

  ggsave(
    file.path(fig_dir, paste0("fig_cluster_wellpct_",
      tolower(gsub("[^A-Za-z0-9]", "_", session_label)),
      "_", analysis_date, ".pdf")),
    p_well_pct,
    width = 10, height = 3 + optimal_k * 2.5,
    device = cairo_pdf
  )
  cat("  Saved: cluster_wellpct\n")

  # Return results and plot objects for figure assembly
  list(
    clust_data       = clust_data,
    cluster_profiles = cluster_profiles,
    composition      = composition,
    well_kw          = well_kw,
    well_clust_pct   = well_clust_pct,
    optimal_k        = optimal_k,
    clust_mat        = clust_mat,
    km_centers       = km_fit$centers,
    p_profile        = p_profile,
    p_pca_clust      = p_pca_clust,
    p_comp           = p_comp
  )
}

# Plot representative traces per cluster (medoid cells closest to centroid).
# Returns the ggplot object invisibly so callers can capture it for figure assembly.
plot_cluster_traces <- function(clust_result, cell_params, results_list,
                                session_label, y_label,
                                n_per_cluster = 8,
                                fig_dir, analysis_date) {
  if (is.null(clust_result)) return(invisible(NULL))

  clust_data  <- clust_result$clust_data
  clust_mat   <- clust_result$clust_mat
  km_centers  <- clust_result$km_centers
  optimal_k   <- clust_result$optimal_k
  cluster_palette <- c("#2166AC", "#D95F02", "#1B9E77",
                       "#E7298A", "#66A61E", "#7570B3")

  set.seed(42)
  trace_df_all <- bind_rows(lapply(seq_len(optimal_k), function(k) {
    # Indices of cells in this cluster (rows of clust_mat)
    idx_k <- which(clust_data$cluster == k)
    if (length(idx_k) == 0) return(NULL)

    # Distance to cluster centroid in scaled space
    center_k <- km_centers[k, ]
    dists <- apply(clust_mat[idx_k, , drop = FALSE], 1, function(row) {
      sqrt(sum((row - center_k)^2))
    })

    # Pick closest n_per_cluster cells (medoids)
    n_pick  <- min(n_per_cluster, length(idx_k))
    picked  <- idx_k[order(dists)[seq_len(n_pick)]]
    picked_uids <- clust_data$cell_uid[picked]

    # Look up traces
    bind_rows(lapply(seq_along(picked_uids), function(j) {
      uid     <- picked_uids[j]
      row_idx <- which(cell_params$cell_uid == uid)[1]
      if (is.na(row_idx)) return(NULL)
      res <- results_list[[row_idx]]
      if (is.null(res$trace) || nrow(res$trace) == 0) return(NULL)
      tibble(
        cluster  = factor(k),
        trace_id = j,
        time_ms  = res$trace$time_ms,
        signal   = res$trace$signal_bc
      )
    }))
  }))

  if (is.null(trace_df_all) || nrow(trace_df_all) == 0) {
    cat("  Warning: no traces found for cluster gallery\n")
    return(invisible(NULL))
  }

  # Mean trace per cluster
  mean_trace <- trace_df_all %>%
    group_by(cluster, time_ms) %>%
    summarise(signal = mean(signal, na.rm = TRUE), .groups = "drop")

  p <- ggplot(trace_df_all,
              aes(x = time_ms, y = signal,
                  group = interaction(cluster, trace_id),
                  color = cluster)) +
    geom_line(alpha = 0.25, linewidth = 0.3) +
    geom_line(data = mean_trace,
              aes(x = time_ms, y = signal, group = cluster),
              linewidth = 1.0, alpha = 1) +
    facet_wrap(~ cluster, nrow = 1,
               labeller = labeller(cluster = function(x)
                 paste0("Cluster ", x,
                        " (n=", table(clust_data$cluster)[x], ")"))) +
    scale_color_manual(values = cluster_palette[seq_len(optimal_k)],
                       guide = "none") +
    labs(title    = paste("Representative traces per cluster:", session_label),
         subtitle = paste0("Thick line = mean of ", n_per_cluster,
                           " medoid cells; thin lines = individual cells"),
         x = "Time (ms)", y = y_label) +
    theme_project() +
    theme(strip.text = element_text(size = 9))

  pdf_path <- file.path(fig_dir,
    paste0("fig_cluster_traces_",
           tolower(gsub("[^A-Za-z0-9]", "_", session_label)),
           "_", analysis_date, ".pdf"))
  ggsave(pdf_path, p,
         width = 4 * optimal_k, height = 4,
         device = cairo_pdf)
  cat("  Saved: cluster_traces ->", pdf_path, "\n")
  invisible(p)
}
