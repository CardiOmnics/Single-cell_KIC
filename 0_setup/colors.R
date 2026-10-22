# Unified Okabe-Ito palette (Okabe & Ito 2008; colorblind-safe gold standard).
# Used by the entire manuscript so every figure speaks the same color language.
group_colors <- c(
  "Control"     = "#0072B2",
  "Palmitate"   = "#009E73",
  "Oleate"      = "#D55E00",
  "Combination" = "#CC79A7"
)

group_levels <- c("Control", "Palmitate", "Oleate", "Combination")

# The plate maps and cached RDS files code the control condition as "MM"
# (maturation medium). Figures, tables and text all call it "Control", so the
# label is swapped once, right after the cell tables are loaded.
relabel_condition <- function(df, col = "Condition") {
  if (col %in% names(df)) {
    vals <- as.character(df[[col]])
    vals[vals == "MM"] <- "Control"
    df[[col]] <- vals
  }
  df
}

# -----------------------------------------------------------------------------
# Subpopulation (k-means cluster) palettes — one family per imaging session
# -----------------------------------------------------------------------------
# Calcium clusters and contraction clusters are independent clusterings of
# different cells; sharing one palette implied a correspondence between
# "cluster 2" in each that does not exist. Each session therefore gets its own
# hue family.
#
# Each family is an ANALOGOUS triad: neighbouring hues that also step down in
# lightness. A single-hue ramp was tried first and rejected — in the dense PCA
# scatter three tints of one hue do not separate, and the stacked bars looked
# flat. Rotating the hue ~35° per step keeps the family calm while giving the
# eye a second cue; the lightness step means the steps still survive colour-blind
# simulation and greyscale printing, where hue alone would not.
#
# Calcium runs cool (aqua → blue → indigo), contraction runs warm
# (gold → orange → red). The cool/warm split makes the two figures
# unmistakably different, which is the point: these are independent clusterings,
# and "cluster 2" in one has nothing to do with "cluster 2" in the other.
#
# Steps are OKLCH L = 0.765 / 0.60 / 0.435. Validated with the dataviz
# six-checks validator (light mode, --pairs all — the strict scatter case where
# any two clusters can sit side by side): both families PASS lightness band,
# chroma floor, CVD separation and the normal-vision floor, worst pair ΔE 18.0
# (calcium) and 17.3 (contraction) against a floor of 15.
#
# Known, accepted overlap: the deep indigo sits near the Control condition blue and
# the burnt orange near the Oleate vermillion. Conditions are never colour-coded
# in the same panel as clusters (the composition bars put Condition on the x
# axis), so the two encodings never meet.
#
# Beyond k = 4 an analogous triad stops separating; widen the hue arc and
# re-validate rather than appending steps.
cluster_colors_calcium     <- c("#58C6CC", "#1B8ABD", "#364998")  # aqua→indigo
cluster_colors_contraction <- c("#D7AB5C", "#C65D26", "#8B2B2F")  # gold→red

# Re-spaced 4-step versions, if k is ever raised to 4.
cluster_colors_calcium_k4     <- c("#58C6CC", "#2C9EC4", "#1C72AC", "#364998")
cluster_colors_contraction_k4 <- c("#D7AB5C", "#CE8340", "#C65D26", "#8B2B2F")

# Diverging scale for the cluster-profile heatmap (mean z-score per parameter).
# Deliberately low-chroma: it encodes magnitude, not identity, so it should
# recede behind the cluster colours rather than compete with them. The previous
# saturated blue↔orange also reused the Control and Oleate condition hues.
zscore_diverging <- c(low = "#6E90AD", mid = "#F7F7F5", high = "#C0885A")

# Pick the label ink (near-black or white) that contrasts better with a given
# fill. Needed because these palettes vary in lightness by design: a fixed white
# label is illegible on the lightest step of a ramp. Uses WCAG relative
# luminance, so it stays correct if the palettes are re-stepped.
contrast_ink <- function(fill, dark = "#1A1A1A", light = "#FFFFFF") {
  vapply(fill, function(f) {
    rgb_lin <- grDevices::col2rgb(f)[, 1] / 255
    rgb_lin <- ifelse(rgb_lin <= 0.03928, rgb_lin / 12.92,
                      ((rgb_lin + 0.055) / 1.055)^2.4)
    L <- sum(c(0.2126, 0.7152, 0.0722) * rgb_lin)
    if ((L + 0.05) / 0.05 > 1.05 / (L + 0.05)) dark else light
  }, character(1), USE.NAMES = FALSE)
}