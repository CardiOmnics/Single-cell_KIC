library(ggplot2)

theme_project <- function() {
  theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 14),
      strip.background = element_blank()
    )
}

# Publication-quality theme for journal figures (Nature-family style).
# Font sizes: axis titles 9 pt, axis text 7 pt, strip labels 9 pt bold.
# Panel labels (A, B, C) are added via cowplot::plot_grid(label_size = 9,
# label_fontface = "bold") or patchwork plot_annotation(tag_levels = "A") with
# theme(plot.tag = element_text(size = 9, face = "bold")).
journal_theme <- theme_bw(base_size = 9) +
  theme(
    axis.text        = element_text(size = 7),
    axis.title       = element_text(size = 9),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 7),
    strip.text       = element_text(size = 9, face = "bold"),
    strip.background = element_blank(),
    axis.line        = element_line(linewidth = 0.5),
    axis.ticks       = element_line(linewidth = 0.5),
    panel.border     = element_rect(linewidth = 0.5, fill = NA)
  )
