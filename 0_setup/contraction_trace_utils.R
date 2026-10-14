# ============================================================
# Utilities for parsing and analysing cardiomyocyte traces
# Contraction (Calcium.Stretch) + Calcium (Calcium.API)
# ============================================================

# ---- 1. Parse serialized trace string ----
parse_trace_string <- function(x) {
  if (is.na(x) || x == "") return(NULL)

  x <- gsub("^\\{|\\}$", "", x)
  parts <- strsplit(x, "\\|")[[1]]

  header <- strsplit(parts[1], ";")[[1]]
  signal <- as.numeric(strsplit(parts[2], ";")[[1]])

  dt_ms <- as.numeric(header[4])
  time_ms <- seq(0, by = dt_ms, length.out = length(signal))

  data.frame(
    time_ms = time_ms,
    signal  = signal
  )
}

# ---- 2. Preprocess trace (detrend + baseline + smoothing) ----
preprocess_trace <- function(df, smooth_k = 5, detrend = TRUE) {
  if (smooth_k > 1) {
    df$signal <- stats::filter(
      df$signal,
      rep(1 / smooth_k, smooth_k),
      sides = 2
    )
    df$signal <- as.numeric(df$signal)
  }

  # Remove NAs introduced by smoothing before detrending

  df <- df[!is.na(df$signal), ]

  # Linear detrending: removes photobleaching / drift

  if (detrend && nrow(df) > 2) {
    fit <- lm(signal ~ time_ms, data = df)
    df$signal <- df$signal - fitted(fit) + mean(df$signal, na.rm = TRUE)
  }

  baseline <- stats::quantile(df$signal, 0.1, na.rm = TRUE)
  df$signal_bc <- df$signal - baseline

  df <- df[!is.na(df$signal_bc), ]
  df
}

# ---- 3. Beat / peak detection ----
detect_beats <- function(df,
                         min_prominence = 0.2,
                         min_distance_ms = 800) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("Package 'pracma' is required for peak detection")
  }

  df <- df[!is.na(df$signal_bc), ]

  # Need at least a few points
  if (nrow(df) < 3) return(NULL)

  fs <- 1000 / median(diff(df$time_ms))
  min_distance <- round((min_distance_ms / 1000) * fs)

  peaks <- pracma::findpeaks(
    df$signal_bc,
    minpeakheight  = min_prominence * max(df$signal_bc, na.rm = TRUE),
    minpeakdistance = min_distance
  )

  if (is.null(peaks)) return(NULL)

  peaks_df <- data.frame(
    peak_amp  = peaks[, 1],
    peak_idx  = peaks[, 2],
    peak_time = df$time_ms[peaks[, 2]]
  )

  # CRUCIAL: sort by time (or index) so diff(peak_time) is always positive
  peaks_df <- peaks_df[order(peaks_df$peak_idx), ]

  # Optional: drop duplicate indices (rare, but safe)
  peaks_df <- peaks_df[!duplicated(peaks_df$peak_idx), ]

  peaks_df
}

# ---- 4. Extract contraction parameters ----
extract_paced_contraction_parameters <- function(df, peaks) {
  if (is.null(peaks) || nrow(peaks) < 2) return(NULL)

  ibi_ms <- diff(peaks$peak_time)

  # Robust noise estimate from the baseline-corrected signal
  noise_mad <- stats::mad(df$signal_bc, constant = 1, na.rm = TRUE)
  mean_amp  <- mean(peaks$peak_amp, na.rm = TRUE)

  # SNR-like metric
  snr_amp <- mean_amp / (noise_mad + 1e-12)

  # Beat kinetics (rise/decay): take robust summary across detected beats
  kin <- extract_beat_kinetics(df, peaks, rise_frac = 0.2, decay_frac = 0.2)

  rise_time_ms  <- if (!is.null(kin)) stats::median(kin$rise_time_ms, na.rm = TRUE) else NA_real_
  decay_time_ms <- if (!is.null(kin)) stats::median(kin$decay_time_ms, na.rm = TRUE) else NA_real_

  data.frame(
    mean_ibi_ms    = mean(ibi_ms),
    median_ibi_ms  = stats::median(ibi_ms),
    sd_ibi_ms      = stats::sd(ibi_ms),
    ibi_cv         = stats::sd(ibi_ms) / mean(ibi_ms),

    mean_amp       = mean_amp,
    sd_amp         = stats::sd(peaks$peak_amp, na.rm = TRUE),
    n_beats        = nrow(peaks),

    noise_mad      = noise_mad,
    snr_amp        = snr_amp,

    rise_time_ms   = rise_time_ms,
    decay_time_ms  = decay_time_ms
  )
}

# ---- 5. One wrapper for a single well ----
analyze_well_trace <- function(x, normalized = FALSE) {
  df <- parse_trace_string(x)
  if (is.null(df)) return(NULL)

  if (!normalized) {
    df <- preprocess_trace(df)
  } else {
    df$signal_bc <- df$signal
  }

  peaks <- detect_beats(df)
  params <- extract_paced_contraction_parameters(df, peaks)

  list(
    trace  = df,
    peaks  = peaks,
    params = params
  )
}

# ---- 6. Quick diagnostic plot ----
plot_well_trace <- function(res) {
  plot(
    res$trace$time_ms,
    res$trace$signal_bc,
    type = "l",
    xlab = "Time (ms)",
    ylab = "Contraction"
  )

  if (!is.null(res$peaks)) {
    points(
      res$peaks$peak_time,
      res$peaks$peak_amp,
      col = "red",
      pch = 16
    )
  }
}

# ============================================================
# Calcium-specific functions
# ============================================================

# ---- 7. Calcium transient duration at a given recovery fraction ----
# For each beat: time from upstroke crossing to decay crossing at `frac` of peak amplitude
calc_ctd <- function(df, peaks, frac = 0.5) {
  if (is.null(peaks) || nrow(peaks) == 0) return(rep(NA_real_, nrow(peaks)))

  df <- df[!is.na(df$signal_bc) & !is.na(df$time_ms), ]
  threshold_frac <- 1 - frac  # CTD50 means 50% recovery, so threshold = 50% of peak

  sapply(seq_len(nrow(peaks)), function(i) {
    peak_idx <- peaks$peak_idx[i]
    peak_amp <- peaks$peak_amp[i]

    if (is.na(peak_idx) || peak_idx < 2 || peak_idx > nrow(df)) return(NA_real_)
    if (is.na(peak_amp) || peak_amp <= 0) return(NA_real_)

    thr <- threshold_frac * peak_amp

    # Upstroke crossing: last point below threshold before peak
    pre <- df$signal_bc[1:peak_idx]
    up_candidates <- which(pre <= thr)
    if (length(up_candidates) == 0) return(NA_real_)
    up_idx <- tail(up_candidates, 1)

    # Decay crossing: first point below threshold after peak
    post <- df$signal_bc[peak_idx:nrow(df)]
    dn_candidates <- which(post <= thr)
    if (length(dn_candidates) == 0) return(NA_real_)
    dn_idx <- peak_idx + dn_candidates[1] - 1

    df$time_ms[dn_idx] - df$time_ms[up_idx]
  })
}

# ---- 8. Extract calcium transient parameters ----
# Expects df$signal_bc = CyteSeer Trimmed trace (already normalised).
# mean_amp is the amplitude metric — no separate dF/F0 is computed.
extract_calcium_parameters <- function(df, peaks) {
  if (is.null(peaks) || nrow(peaks) < 2) return(NULL)

  ibi_ms <- diff(peaks$peak_time)

  noise_mad <- stats::mad(df$signal_bc, constant = 1, na.rm = TRUE)
  mean_amp  <- mean(peaks$peak_amp, na.rm = TRUE)
  snr_amp   <- mean_amp / (noise_mad + 1e-12)

  # CTD at 50% and 90% recovery
  ctd50_vals <- calc_ctd(df, peaks, frac = 0.5)
  ctd90_vals <- calc_ctd(df, peaks, frac = 0.9)

  # Rise (10% to 90%) and decay (peak to 80% recovery) using existing kinetics function
  kin_rise <- extract_beat_kinetics(df, peaks, rise_frac = 0.1, decay_frac = 0.1)
  kin_decay <- extract_beat_kinetics(df, peaks, rise_frac = 0.1, decay_frac = 0.2)

  rise_time_ms  <- if (!is.null(kin_rise))
    stats::median(kin_rise$rise_time_ms, na.rm = TRUE) else NA_real_
  decay_time_ms <- if (!is.null(kin_decay))
    stats::median(kin_decay$decay_time_ms, na.rm = TRUE) else NA_real_

  data.frame(
    mean_ibi_ms   = mean(ibi_ms),
    median_ibi_ms = stats::median(ibi_ms),
    sd_ibi_ms     = stats::sd(ibi_ms),
    ibi_cv        = stats::sd(ibi_ms) / mean(ibi_ms),

    mean_amp      = mean_amp,
    sd_amp        = stats::sd(peaks$peak_amp, na.rm = TRUE),
    n_beats       = nrow(peaks),

    noise_mad     = noise_mad,
    snr_amp       = snr_amp,

    CTD50_ms      = stats::median(ctd50_vals, na.rm = TRUE),
    CTD90_ms      = stats::median(ctd90_vals, na.rm = TRUE),
    rise_time_ms  = rise_time_ms,
    decay_time_ms = decay_time_ms
  )
}

# ---- 9. Single-cell calcium wrapper ----
# Expects the CyteSeer "Trimmed.Ci.on.Cm.API.F.t." trace, which is already a
# normalised dF/F-style signal. We do NO preprocessing here — no smoothing,
# no detrending, no baseline subtraction.
analyze_cell_calcium <- function(x) {
  df <- parse_trace_string(x)
  if (is.null(df)) return(NULL)

  # No preprocessing: signal_bc is the trace as delivered by CyteSeer.
  df$signal_bc <- df$signal
  df <- df[!is.na(df$signal_bc) & !is.na(df$time_ms), ]

  peaks <- detect_beats(df, min_prominence = 0.15, min_distance_ms = 700)
  params <- extract_calcium_parameters(df, peaks)

  list(
    trace  = df,
    peaks  = peaks,
    params = params
  )
}

# ============================================================
# Shared helper: beat kinetics (used by both contraction and calcium)
# ============================================================

extract_beat_kinetics <- function(df, peaks,
                                  rise_frac = 0.2,
                                  decay_frac = 0.2) {
  if (is.null(peaks) || nrow(peaks) == 0) return(NULL)

  # Ensure no NAs
  df <- df[!is.na(df$signal_bc) & !is.na(df$time_ms), ]

  kinetics <- lapply(seq_len(nrow(peaks)), function(i) {
    peak_idx <- peaks$peak_idx[i]
    peak_amp <- peaks$peak_amp[i]

    # Guard against weird indices
    if (is.na(peak_idx) || peak_idx < 2 || peak_idx > nrow(df)) return(NULL)
    if (is.na(peak_amp) || peak_amp <= 0) return(NULL)

    rise_thr  <- rise_frac  * peak_amp
    decay_thr <- decay_frac * peak_amp

    # Search backward for last point below threshold
    pre <- df$signal_bc[1:peak_idx]
    rise_candidates <- which(pre <= rise_thr)
    if (length(rise_candidates) == 0) return(NULL)
    rise_idx <- tail(rise_candidates, 1)

    # Search forward for first point below threshold
    post <- df$signal_bc[peak_idx:length(df$signal_bc)]
    decay_candidates <- which(post <= decay_thr)
    if (length(decay_candidates) == 0) return(NULL)
    decay_offset <- decay_candidates[1]
    decay_idx <- peak_idx + decay_offset - 1

    data.frame(
      rise_time_ms  = df$time_ms[peak_idx] - df$time_ms[rise_idx],
      decay_time_ms = df$time_ms[decay_idx] - df$time_ms[peak_idx]
    )
  })

  out <- do.call(rbind, kinetics)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out
}

# ============================================================
# QC: pacing fidelity score + data-driven threshold
# ============================================================
# AC1 = normalised autocorrelation of a detrended trace at lag equal to
# one pacing period (default 1 s). Waveform-agnostic single-number score
# that quantifies how well a cell follows the pacemaker. Used together
# with a hard upper cap on n_beats (max beats delivered per recording)
# to classify cells as good_paced vs irregular.
#
# The threshold on AC1 is data-driven: Otsu's method on the per-modality
# distribution picks the cutoff at the valley between the noise and
# paced modes, so no hand-picked number enters the QC.

# ---- Pacing autocorrelation score (one trace) ----
# tr: data.frame with columns time_ms and signal_bc
# pace_period_ms: expected beat-to-beat interval (1000 ms at 1 Hz)
score_pacing_ac1 <- function(tr, pace_period_ms = 1000) {
  if (is.null(tr) || !is.data.frame(tr) || nrow(tr) < 8) return(NA_real_)
  t_ms <- tr$time_ms
  y    <- tr$signal_bc
  ok   <- is.finite(t_ms) & is.finite(y)
  t_ms <- t_ms[ok]; y <- y[ok]
  if (length(y) < 8) return(NA_real_)

  # Linear detrend + DC removal — matches the prototype preprocessing
  idx <- seq_along(y)
  y   <- residuals(stats::lm(y ~ idx))
  y   <- y - mean(y)
  if (stats::sd(y) == 0) return(NA_real_)

  fs    <- 1000 / stats::median(diff(t_ms))   # samples per second
  lag_k <- round(pace_period_ms / 1000 * fs)
  if (lag_k < 1 || lag_k >= length(y)) return(NA_real_)
  ac <- stats::acf(y, lag.max = lag_k, plot = FALSE,
                   na.action = stats::na.pass)$acf
  as.numeric(ac[lag_k + 1])
}

# ---- Otsu threshold (data-driven bimodal split) ----
# Picks the cutoff that minimises within-class variance. Suitable when
# the score distribution is bimodal (noise mode + paced mode).
otsu_threshold <- function(x, nbins = 256) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NA_real_)
  rng <- range(x)
  if (diff(rng) == 0) return(NA_real_)
  br    <- seq(rng[1], rng[2], length.out = nbins + 1)
  h     <- graphics::hist(x, breaks = br, plot = FALSE)
  p     <- h$counts / sum(h$counts)
  mids  <- h$mids
  omega <- cumsum(p)
  mu    <- cumsum(p * mids)
  mu_t  <- mu[length(mu)]
  sigma_b2 <- (mu_t * omega - mu)^2 / (omega * (1 - omega) + 1e-12)
  mids[which.max(sigma_b2)]
}
