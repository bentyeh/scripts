library(tidyverse)

mzR_to_MassSpectrum <- function(spectrum) {
  # Create MALDIquant::MassSpectrum from mzR spectrum
  #
  # Args
  # - spectrum: matrix
  #     Mass spectrum, such as that returned by mzR::peaks(mzR_object)
  #     - Column 1: mass value
  #     - Column 2: intensity
  # 
  # Returns: MassSpectrum
  require(MALDIquant)
  MALDIquant::createMassSpectrum(
    mass = spectrum[,1],
    intensity = spectrum[,2]
  )
}

mzR_to_MSImagingExperiment <- function(spectrum) {
  # Create Cardinal::MSImagingExperiment from mzR spectrum
  #
  # Args
  # - spectrum: matrix
  #     Mass spectrum, such as that returned by mzR::peaks(mzR_object)
  #     - Column 1: mass value
  #     - Column 2: intensity
  # 
  # Returns: MSImagingExperiment
  require(Cardinal)
  Cardinal::MSImagingExperiment(
    imageData = ImageArrayList(matrix(spectrum[,2], nrow = nrow(spectrum), ncol = 1)),
    featureData = MassDataFrame(mz = spectrum[,1]),
    pixelData = PositionDataFrame(coord = expand.grid(x = 1, y = 1))
  )
}

mzR_to_MSImageSet <- function(spectrum) {
  # Create Cardinal::MSImageSet from mzR spectrum
  #
  # Args
  # - spectrum: matrix
  #     Mass spectrum, such as that returned by mzR::peaks(mzR_object)
  #     - Column 1: mass value
  #     - Column 2: intensity
  # 
  # Returns: MSImageSet
  require(Cardinal)
  Cardinal::MSImageSet(
    spectra = matrix(spectrum[,2]),
    mz = spectrum[,1]
  )
}

getMassPeak <- function(spectrum, mass, tol_ppm = 1e3) {
  # Get mass and intensity of peak with highest intensity near a mass of interest.
  # 
  # Args
  # - spectrum: matrix
  #     Mass spectrum, such as that returned by mzR::peaks(mzR_object)
  #     - Column 1: mass value
  #     - Column 2: intensity
  # - mass: numeric
  #     Mass of peak of interest
  # - tol_ppm: numeric. default=1e3
  #     Mass tolerance "window". Look for peaks within mass * tol_ppm/1e6 of mass.
  # 
  # Returns: vector, numeric of length 2
  #   Row in `spectrum` corresponding to mass peak

  tol <- mass * tol_ppm / 1e6
  idx <- which((spectrum[, 1] > mass - tol) & (spectrum[, 1] < mass + tol))
  mass_peak_idx <- idx[which.max(spectrum[idx, 2])]
  spectrum[mass_peak_idx, ]
}

noiseRunMed <- function(spectrum, mass, tol_ppm = 1e4, summary = median) {
  # Estimate the noise around a mass value of interest as a function (e.g.,
  # median) of a running median of intensities in a region around that mass
  # value.
  # 
  # Args
  # - spectrum: matrix
  #     Mass spectrum, such as that returned by mzR::peaks(mzR_object)
  #     - Column 1: mass value
  #     - Column 2: intensity
  # - mass: numeric
  #     Mass of peak of interest
  # - tol_ppm: numeric. default=1e4
  #     Mass tolerance "window"
  # - summary: function. default = median
  #     Must take a numeric vector as an argument and return a single numeric value.
  #     Examples: min, max, mean, median
  # 
  # Returns: numeric
  #   Estimated noise intensity
  # 
  # Alternatives
  # - stats::smooth(): Tukey's running median with a window of 3
  # - stats::lowess(): locally-weighted polynomial regression smoother

  tol <- mass * tol_ppm/1e6
  idx <- which((spectrum[,1] > mass - tol) &
               (spectrum[,1] < mass + tol))
  avg_resolution <- 2 * tol / length(idx)
  width <- ceiling(tol/avg_resolution)
  if (width %% 2 == 0) {
    width <- width + 1
  }
  summary(runmed(spectrum[idx, 2], width))
}

calibrate_spectrum <- function(
  spectrum,
  ref_masses,
  SNR = 10,
  noise_method = "runmed",
  min_peaks = NULL,
  tol = 1e3,
  tol_units = "ppm",
  model = "auto",
  bg = NULL,
  bg_ignore = NULL,
  force = NULL,
  verbose = TRUE) {
  # Calibrate spectrum based on reference masses.
  # 
  # Args
  # - spectrum: matrix
  #     Mass spectrum, such as that returned by mzR::peaks(mzR_object)
  #     - Column 1: mass value
  #     - Column 2: intensity
  # - ref_masses: vector, numeric
  #     Peaks to align to
  # - SNR: numeric. default = 10
  #     Signal to noise ratio threshold for peak detection
  # - noise_method: character or function. default = "runmed"
  #     Noise estimation method.
  #     If character: uses the following functions with default arguments
  #     - "runmed": noiseRunMed()
  #     - "simple": Cardinal:::.estimateNoiseSimple()
  #     - "adaptive": Cardinal:::.estimateNoiseAdapative()
  #     - "mad": Cardinal:::.estimateNoiseMAD()
  #     If function: must be of the form f(spectrum, mass) -> numeric
  # - min_peaks: integer. default = NULL
  #     Minimum required number of peaks in ref_masses to detect in spectrum. If not satisfied,
  #     returns the spectrum as is (or returns the identity function if `return_function`` is TRUE.)
  #     Defaults to length(ref_masses).
  # - tol: numeric. default = 1e3. length = 1 or length(ref_masses)
  #     Mass tolerance window for peak identification.
  # - tol_units: character. default = "ppm". length = length(tol)
  #     "mz" or "ppm"
  # - model: character or function. default = "auto"
  #     Built-in options: "loess", "lm", "spline", "quad", "auto", "identity"
  #       - "auto" uses "lm" if number of matched peaks
  #     To use a custom model, pass in a function that takes a data.frame (x = peaks, y = ref_masses)
  #       as an argument and returns a vectorized function mapping original to calibrated masses.
  # - bg: matrix. default = NULL
  #     Mass spectrum, such as that returned by mzR::peaks(mzR_object)
  #     - Column 1: mass value
  #     - Column 2: intensity
  #     Background mass spectrum for additional validation.
  # - bg_ignore: vector, logical. default = NULL
  #     Reference masses to ignore in the background (e.g., if using the matrix peak as a reference mass)
  # - force: vector, logical. default = NULL
  #     Force use of certain reference mass peaks
  # - verbose: logical. default = TRUE
  # 
  # Returns: list
  # - warp: function
  #     Calibration function mapping numeric vector of uncalibrated masses to a numeric vector of calibrated masses
  # - spectrum: matrix
  #     Calibrated spectrum
  # - noise_method: character or function
  #     (same as argument)
  # - min_peaks: integer
  #     Resolved minimum number of peaks (i.e., if min_peaks argument was NULL)
  # - SNR: numeric
  #     (same as argument)
  # - model: character or function
  #     Resolved model (i.e., if model argument was "auto", or min_peaks threshold was not satisfied)
  # - df: data.frame
  #     ref_mass (numeric): reference masses
  #     detected_mass (numeric): detected mass peaks of reference masses
  #     tol: vector, numeric
  #       Tolerance used for each reference mass in ppm
  #     snr (numeric): SNR of detected mass peaks of reference masses, relative to the local region in the
  #       calibration spectrum
  #     detected_mass_bg (numeric; only present if bg != NULL): detected mass peaks of reference masses in
  #       background spectrum
  #     snr_bg (numeric; only present if bg != NULL): SNR of detected mass peaks of reference masses,
  #       relative to the detected mass peaks in the background spectrum
  #     passed (logical): whether the reference mass peak was ultimately used for calibration
  #     residuals_uncalibrated (numeric): residuals between detected_mass and reference_mass
  #     residuals_calibrated (numeric): residuals after calibration
  # - mse_uncal: numeric
  #     Mean squared error of uncalibrated reference masses
  # - mse_cal: numeric
  #     Mean squared error of calibrated reference masses
  # - model_summary: summary.lm or summary.loess. Only present if model is "lm", "loess", or "quad"
  #     Model summary as returned by summary() on a lm() or loess() model object
  # 
  # Notes
  # - Arguments inspired by AB/SCIEX 5800 Series Explorer's Processing Method parameters.
  # - Calibration algorithm based on Cardinal::mzAlign(). See
  #   https://rdrr.io/bioc/Cardinal/src/R/process2-mzAlign.R

  if (is.null(min_peaks)) {
    min_peaks <- length(ref_masses)
  }
  min_peaks <- as.integer(min_peaks)
  if (is.null(bg_ignore)) {
    bg_ignore <- rep(FALSE, length(ref_masses))
  }
  if (is.null(force)) {
    force <- rep(FALSE, length(ref_masses))
  }
  stopifnot(tol_units %in% c("mz", "ppm"))
  stopifnot(min_peaks <= length(ref_masses))
  stopifnot(length(tol) == length(tol_units))
  stopifnot(length(tol) %in% c(1, length(ref_masses)))
  stopifnot(rlang::is_callable(noise_method) || noise_method %in% c("simple", "adaptive", "mad", "runmed"))
  stopifnot(rlang::is_callable(model) || model %in% c("loess", "lm", "quad", "spline", "auto", "identity"))

  noise_method_fun <- noise_method
  if (identical(noise_method_fun, "runmed")) {
    noise_method_fun <- noiseRunMed
  }

  n <- length(ref_masses)

  # convert all tolerances to units "ppm"
  if (length(tol) != n) {
    tol <- rep(tol, n)
    tol_units <- rep(tol_units, n)
  }
  for (i in 1:n) {
    if (tol_units[i] != "ppm") {
      tol[i] <- tol[i] / ref_masses[i] * 1e6
      tol_units[i] <- "ppm"
    }
  }

  peaks <- matrix(NA_real_, nrow = n, ncol = 2)
  for (i in 1:n) {
    peaks[i, ] <- getMassPeak(spectrum, ref_masses[i], tol_ppm = tol[i])
  }

  noise <- vector("numeric", n)
  if (rlang::is_callable(noise_method_fun)) {
    for (i in 1:n) {
      noise[i] <- noise_method_fun(spectrum, peaks[i, 1])
    }
  } else {
    noise_method_fun <- list(
      simple = Cardinal:::.estimateNoiseSimple,
      adaptive = Cardinal:::.estimateNoiseAdaptive,
      mad = Cardinal:::.estimateNoiseMAD
    )[[noise_method_fun]]
    noise_all <- noise_method_fun(spectrum[, 2])
    for (i in 1:n) {
      noise[i] <- noise_all[spectrum[,1] == peaks[i, 1]]
    }
  }
  snr <- peaks[, 2] / noise

  if (is.null(bg)) {
    detected_peaks <- snr >= SNR
    df <- data.frame(
      ref_mass = ref_masses,
      detected_mass = peaks[, 1],
      tol = tol,
      snr = snr,
      passed = detected_peaks
    )
  } else {
    peaks_bg <- matrix(NA_real_, nrow = n, ncol = 2)
    for (i in 1:n) {
      peaks_bg[i, ] <- getMassPeak(bg, ref_masses[i], tol_ppm = tol[i])
    }

    snr <- peaks[, 2] / noise
    snr_bg <- (peaks[, 2] / sum(spectrum[, 2])) / (peaks_bg[, 2] / sum(bg[, 2]))
    detected_peaks <- ((snr >= SNR) & ((snr_bg >= SNR) | bg_ignore)) | force
    n_detected_peaks <- sum(detected_peaks)
    df <- data.frame(
      ref_mass = ref_masses,
      detected_mass = peaks[, 1],
      tol = tol,
      snr = snr,
      detected_mass_bg = peaks_bg[, 1],
      snr_bg = snr_bg,
      passed = detected_peaks
    )
  }
  if (tibble::has_rownames(df)) {
    df <- df %>% tibble::rownames_to_column(var = "ref_name")
  }
  n_detected_peaks <- sum(detected_peaks)

  if (verbose) {
    print(stringr::str_glue("Detected {n_detected_peaks} / {n} reference peaks."))
  }

  if (identical(model, "auto") || n_detected_peaks < min_peaks) {
    if (n_detected_peaks < min_peaks) {
      message(stringr::str_glue(
        "Number of detected peaks ({n_detected_peaks}) < min_peaks ({min_peaks}). Not calibrating."
      ))
      model <- "identity"
    } else {
      model <- ifelse(n_detected_peaks > 3, "spline", "lm")
    }
    if (verbose) {
      print(stringr::str_glue("Using {model} model."))
    }
  }

  ref_masses <- ref_masses[detected_peaks]
  peaks <- peaks[detected_peaks, , drop = FALSE]

  df_xy <- data.frame(x = peaks[, 1], y = ref_masses)
  model_obj <- NULL
  if (is.character(model) && model %in% c("loess", "lm", "quad")) {
    if (model == "loess") {
      model_obj <- loess(y ~ x, data = df_xy, span = 1, control = loess.control(surface = "direct"))
    } else if (model == "lm") {
      model_obj <- lm(y ~ x, data = df_xy)
    } else {
      model_obj <- lm(y ~ poly(x, 2), data = df_xy)
    }
    warp <- function(x) pmax(predict(model_obj, data.frame(x = x)), 0)
  } else if (identical(model, "spline")) {
    warp <- splinefun(peaks[, 1], ref_masses, "natural")
  } else if (identical(model, "identity")) {
    warp <- identity
  } else {
    warp <- model(df_xy)
  }
  spectrum[, 1] <- warp(spectrum[, 1])

  df[["residuals_uncalibrated"]] <- NA_real_
  df[["residuals_uncalibrated"]][detected_peaks] <- peaks[, 1] - ref_masses
  df[["residuals_calibrated"]] <- NA_real_
  df[["residuals_calibrated"]][detected_peaks] <- warp(peaks[,1]) - ref_masses
  if (n_detected_peaks > 0) {
    mse_uncal <- sum(df[["residuals_uncalibrated"]] ^ 2, na.rm = TRUE) / n_detected_peaks
    mse_cal <- sum(df[["residuals_calibrated"]] ^ 2, na.rm = TRUE) / n_detected_peaks
  } else {
    mse_uncal <- NA
    mse_cal <- NA
  }

  if (verbose) {
    print(df)
    print(stringr::str_glue("MSE (uncalibrated): {mse_uncal}"))
    print(stringr::str_glue("MSE (calibrated): {mse_cal}"))
  }

  retval <- list(
    warp = warp,
    spectrum = spectrum,
    min_peaks = min_peaks,
    noise_method = noise_method,
    model = model,
    df = df,
    SNR = SNR,
    mse_uncal = mse_uncal,
    mse_cal = mse_cal
  )

  if (!is.null(model_obj)) {
    retval[["model_summary"]] <- summary(model_obj)
    if (verbose) {
      print(retval[["model_summary"]])
    }
  }

  return(retval)
}

sliding_distance <- function(v1, v2, norm_fun = "1", range = NULL, verbose = FALSE) {
  # Compute distance between 2 vectors across a range of sliding offsets.
  # 
  # Offsets are given as the number of positions v2 is slided right relative to v1.
  # Distances are given for the aligned subsets of v1 and v2 according to a distance metric.
  # 
  # Args
  # - v1: vector of numeric
  #     
  # - v2: vector of numeric
  #     Must be of the same length as `v1`
  # - norm: character or callable. default = "1"
  #     Distance metric.
  #     If character, specifies the norm type as used by base::norm()
  #     If callable, a function used to compute a vector norm.
  # - range: vector of integer, length = 2. default = NULL
  #     Range of offsets to try. If NULL, defaults to all possible offets.
  # - verbose: logical, default = FALSE
  # 
  # Returns: numeric. length = range[2] - range[1] + 1. names = range[1]:range[2]
  #   Distance between aligned subsets of vectors v1 and v2 at different offsets.

  n <- length(v1)
  stopifnot(n == length(v2))
  if (is.null(range)) {
    range <- c(-n + 1, n - 1)
  }
  stopifnot(length(range) == 2)
  stopifnot(all(abs(range) < n))
  stopifnot(all(range%%1 == 0))
  stopifnot(range[2] >= range[1])

  if (is.character(norm_fun)) {
    type <- norm_fun
    norm_fun <- function(x) base::norm(as.matrix(x), type = type)/length(x)
  }
  stopifnot(rlang::is_callable(norm_fun))

  d <- numeric(range[2] - range[1] + 1)
  names(d) <- range[1]:range[2]
  for (i in range[1]:range[2]) {
    sub_v1 <- v1[max(1, i + 1):min(n + i, n)]
    sub_v2 <- v2[max(1, 1 - i):min(n, n - i)]
    stopifnot(length(sub_v1) == length(sub_v2))
    d[as.character(i)] <- norm_fun(sub_v1 - sub_v2)
  }
  if (verbose) {
    offsets <- as.integer(names(d)[which(d == min(d))])
    offset <- offsets[which.min(abs(offsets))]
    v1_aligned <- v1[max(1, offset + 1):min(n + offset, n)]
    v2_aligned <- v2[max(1, 1 - offset):min(n, n - offset)]
    print("Difference after alignment")
    print(summary(v1_aligned - v2_aligned))
  }
  return(d)
}

align_vectors <- function(...) {
  # Minimum offset necessary to optimally align 2 vectors. A positive offset indicates that the
  # second vector should be slid right relative to the first vector.
  # 
  # offset = align_vectors(v1, v2)
  # n = length(v1)
  # v1_aligned = v1[max(1, offset + 1):min(n + offset, n)]
  # v2_aligned = v2[max(1, 1 - offset):min(n, n - offset)]

  d <- sliding_distance(...)
  offsets <- as.integer(names(d)[which(d == min(d))])
  return(offsets[which.min(abs(offsets))])
}

align_spectra <- function(spectra_list, trim = TRUE, ...) {
  # Subset each spectrum such that the subset masses are as closely aligned as possible.
  # 
  # Args
  # - spectra_list: list of 2-column matrices
  #     First spectrum in the list is assumed to be the "reference" spectrum
  # - trim: logical. default = TRUE
  #     Keep only the region where all spectra are aligned. Otherwise, fill in unaligned regions
  #     with NA.
  # - ...: arguments to pass to align_vectors()
  # 
  # Returns: list of 2-column matrices
  #   Each spectrum subset and aligned to first spectrum.
  #   Note: masses (first column of each spectra) are not necessarily identical.

  # assert that all spectra have the same number of rows
  stopifnot(identical(purrr::map_int(spectra_list, nrow) %>% unique() %>% length(), 1L))

  n <- nrow(spectra_list[[1]])
  offsets <- c(
    0,
    purrr::map_int(
      spectra_list[2:length(spectra_list)],
      function(x, ...) align_vectors(x[,1], spectra_list[[1]][,1], ...),
      ...
    )
  )
  if (trim) {
    max_offset <- max(offsets, 0)
    min_offset <- min(offsets, 0)
    for (i in seq_along(spectra_list)) {
      # spectra_list[[i]] = spectra_list[[i]][max(1, offsets[i] + 1):min(n + offsets[i], n),]
      spectra_list[[i]] = spectra_list[[i]][(max_offset - offsets[i] + 1):(n - offsets[i] + min_offset),]
    }
  }
  return(spectra_list)
}

aggregate_spectra <- function(spectra_list, agg_fun = median, verbose = TRUE, ...) {
  # Args
  # - spectra_list: list of 2-column matrices
  #     First spectrum in the list is assumed to be the "reference" spectrum
  # - agg_fun: callable
  #     The callable should take a single vector argument and return a numeric scalar
  # - verbose: logical
  # - ...: arguments to pass to align_spectra()

  stopifnot(rlang::is_callable(agg_fun))

  aligned_spectra <- align_spectra(spectra_list, ...)
  stopifnot(identical(purrr::map_int(aligned_spectra, nrow) %>% unique() %>% length(), 1L))
  masses <- purrr::map(aligned_spectra, ~ .x[,1]) %>%
    purrr::pmap_dbl(function(...) agg_fun(c(...)))
  intensities <- purrr::map(aligned_spectra, ~ .x[,2]) %>%
    purrr::pmap_dbl(function(...) agg_fun(c(...)))
  return(unname(cbind(masses, intensities)))
}

focusSpectrum <- function(spectrum, min = NULL, max = NULL, mass = NULL, tol_ppm = 1e4) {
  # Extract region of spectrum based on min/max m/z values or a given mass and tolerance (in ppm).
  # 
  # Args
  # - spectrum: matrix
  #     Mass spectrum, such as that returned by mzR::peaks(mzR_object)
  #     - Column 1: mass value
  #     - Column 2: intensity
  # - min: numeric. default = NULL
  #     Minimum m/z value to plot
  # - max: numeric. default = NULL
  #     Maximum m/z value to plot
  # - mass: numeric. default = NULL
  #     Mass of peak of interest.
  # - tol_ppm: numeric. default=1e4
  #     Mass tolerance "window."
  # 
  # Notes: The arguments min and max take precedence over mass and tol_ppm.
  # 
  # Returns: matrix

  if (!is.null(min)) {
    spectrum <- spectrum[spectrum[,1] >= min,]
  }

  if (!is.null(max)) {
    spectrum <- spectrum[spectrum[,1] <= max,]
  }

  if (is.null(min) && is.null(max) && !is.null(mass)) {
    tol <- mass * tol_ppm/1e6
    idx <- which((spectrum[,1] > mass - tol) &
                   (spectrum[,1] < mass + tol))
    spectrum <- spectrum[idx,]
  }

  return(spectrum)
}

plotSpectrum <- function(spectrum, min = NULL, max = NULL, mass = NULL, tol_ppm = 1e4,
                         show_noise = FALSE, noise_col = "red", plotter = "baseR", overlay = FALSE,
                         ...) {
  # Create basic plot (bar chart) of mass spectrum. By default, labels the
  # x-axis "m/z" and the y-axis "count".
  # 
  # Args
  # - spectrum: matrix
  #     Mass spectrum, such as that returned by mzR::peaks(mzR_object)
  #     - Column 1: mass value
  #     - Column 2: intensity
  # - min, max, mass, tol_ppm: numeric.
  #     See focusSpectrum().
  # - show_noise: bool or numeric. default = TRUE
  #     TRUE: compute and plot noise line
  #     numeric: plot given noise line
  # - plotter: character. default = "baseR"
  #     "baseR": create base R plot
  #     "ggplot": create ggplot2 plot
  # - overlay: logical or ggplot. default = FALSE
  #     Overlay current plot on existing plot.
  # - ...
  #     Arguments to pass to base R plotter or aes() within geom_col()
  # 
  # Returns: NULL (`plotter` == 'baseR') or ggplot (`plotter` == 'ggplot')
  # 
  # Examples
  # 1. Overlay 2 spectra with different colors using baseR plotter
  #   plotSpectrum(spectrum1, mass = 100, col = "red")
  #   plotSpectrum(spectrum2, mass = 100, col = "blue", overlay = TRUE)
  # 2. Overlay 2 spectra with different colors using ggplot plotter
  #   g <- plotSpectrum(spectrum1, mass = 100, plotter = "ggplot", fill = "red", alpha = 0.5)
  #   g <- plotSpectrum(spectrum2, mass = 100, plotter = "ggplot", fill = "blue", alpha = 0.5, overlay = g)
  #   g + guides(alpha = "none", fill = guide_legend(title = "spectra"))

  spectrum <- focusSpectrum(spectrum, min, max, mass, tol_ppm)
  xlim <- c(min(spectrum[,1]), max(spectrum[,1]))

  # default parameters
  default_args <- list(
    xlab = "m/z",
    ylab = "count",
    main = NULL,
    type = "h", # 'histogram' like vertical lines,
    xlim = xlim,
    ylim = c(0, max(spectrum[,2])),
    yaxs = "i", # tightly fit axis to data range
    xaxs = "i", # tightly fit axis to data range
    bty = "l"   # draw only the left and bottom axes lines
  )

  # set parameters based on input arguments
  args <- list(...)
  default_args_to_use <- names(default_args)[!names(default_args) %in% names(args)]
  args[default_args_to_use] <- default_args[default_args_to_use]

  args[["x"]] <- spectrum[,1]
  args[["y"]] <- spectrum[,2]

  if (identical(show_noise, TRUE)) {
    if (is.null(min) && is.null(max) && !is.null(mass)) {
      show_noise <- noiseRunMed(spectrum, mass, tol_ppm)
    } else {
      show_noise <- noiseRunMed(spectrum, mean(xlim), tol_ppm)
    }
  }

  if (plotter == "baseR") {
    if (overlay) {
      tryCatch(do.call(lines, args), error = function(e) do.call(plot, args))
    } else {
      do.call(plot, args)
    }
    if (is.numeric(show_noise)) {
      abline(h = show_noise, col = noise_col)
    }
  } else {
    colnames(spectrum) <- c("m/z", "count")
    if ("ggplot" %in% class(overlay)) {
      g <-
        overlay + 
        geom_col(aes(`m/z`, count, ...), data = as_tibble(spectrum)) + 
        labs(x = "m/z")
    } else {
      g <-
        spectrum %>%
        as_tibble() %>%
        ggplot(aes(`m/z`, count, ...)) +
        geom_col() +
        labs(x = "m/z")
    }
    if (is.numeric(show_noise)) {
      g <- 
        g +
        geom_hline(yintercept = show_noise, color = "red")
    }
    print(g)
    return(g)
  }
}