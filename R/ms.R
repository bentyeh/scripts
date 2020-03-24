library(tidyverse)
library(Cardinal)
library(MALDIquant)

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
  MALDIquant::createMassSpectrum(
    mass = spectrum[,1],
    intensity = spectrum[,2]
  )
}

MassSpectrum_to_mzR <- function(MS) {
  # Create spectrum matrix (column 1 = m/z, column 2 = intensity) from MALDIquant::MassSpectrum
  cbind(MALDIquant::mass(MS), MALDIquant::intensity(MS))
}

add_colnames <- function(spectrum, names = c("m/z", "count")) {
  # Add column names to spectrum matrix
  colnames(spectrum) <- names
  return(spectrum)
}

mzR_to_MSImagingExperiment <- function(spectrum) {
  # Create Cardinal::MSImagingExperiment from mzR spectrum
  Cardinal::MSImagingExperiment(
    imageData = ImageArrayList(matrix(spectrum[,2], nrow = nrow(spectrum), ncol = 1)),
    featureData = MassDataFrame(mz = spectrum[,1]),
    pixelData = PositionDataFrame(coord = expand.grid(x = 1, y = 1))
  )
}

mzR_to_MSImageSet <- function(spectrum) {
  # Create Cardinal::MSImageSet from mzR spectrum
  Cardinal::MSImageSet(
    spectra = matrix(spectrum[,2]),
    mz = spectrum[,1]
  )
}

normalize_tic <- function(spectrum) {
  spectrum[, 2] <- spectrum[, 2] / sum(spectrum[, 2])
  return(spectrum)
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
  stopifnot(
    rlang::is_callable(model) ||
    model %in% c("loess", "lm", "spline", "auto", "identity") ||
    str_detect(model, "^poly[\\d.]+")
  )

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
      model <- ifelse(n_detected_peaks > 3, "loess", "lm")
    }
    if (verbose) {
      print(stringr::str_glue("Using {model} model."))
    }
  }

  ref_masses <- ref_masses[detected_peaks]
  peaks <- peaks[detected_peaks, , drop = FALSE]

  df_xy <- data.frame(x = peaks[, 1], y = ref_masses)
  model_obj <- NULL
  if (identical(model, "loess")) {
    model_obj <- loess(y ~ x, data = df_xy, span = 1, control = loess.control(surface = "direct"))
    warp <- function(x) pmax(predict(model_obj, data.frame(x = x)), 0)
  } else if (identical(model, "lm")) {
    model_obj <- lm(y ~ x, data = df_xy)
    warp <- function(x) pmax(predict(model_obj, data.frame(x = x)), 0)
  } else if (is.character(model) && str_detect(model, "^poly[\\d.]+")) {
    degree <- as.numeric(str_match(model, "poly([\\d.]+)")[2])
    df_xy$z <- df_xy$x ^ degree
    model_obj <- lm(y ~ x + z, data = df_xy)
    warp <- function(x) pmax(predict(model_obj, data.frame(x = x, z = x ^ degree)), 0)
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

bin_spectrum <- function(spectrum, decimals = 2, agg_fun = mean) {
  binned_spectrum <- spectrum
  binned_spectrum[,1] <- round(binned_spectrum[,1], decimals)
  colnames(binned_spectrum) <- c("mass", "intensity")
  binned_spectrum <- binned_spectrum %>%
    as_tibble() %>% 
    group_by(mass) %>% 
    summarise(intensity = agg_fun(intensity)) %>% 
    as.matrix()
  colnames(binned_spectrum) <- colnames(spectrum)
  return(binned_spectrum)
}

snr_mask <- function(mse, mask, mass, tol = 0.5) {
  # Signal to noise ratio of mean intensities over all pixels in the mask at the specified mass
  # 
  # Args
  # - mse: MSImagingExperiment
  # - mask: vector, logical or integer
  #     Logical or index (integer) mask of pixels.
  #     length(mask) must equal ncol(mse)
  # - mass: numeric
  #     Mass of interest. The closest mass in mse is used.
  # - tol: numeric. default = 0.5
  #     A warning is raised if the closest mass in mse is more than tol away from mass.
  # 
  # Returns: numeric
  diff <- Cardinal::mz(mse) - mass
  mass_idx <- which(abs(diff) == min(abs(diff)))[1]
  if (min(abs(diff)) > tol) {
    warning(str_glue("Closest mass in spectra ({Cardinal::mz(mse)[mass_idx]}) is far from desired mass ({mass})."))
  }
  signal <- Cardinal::spectra(mse)[mass_idx, mask] %>% mean()
  noise <-
    cbind(
      Cardinal::mz(mse),
      Cardinal::spectra(mse)[, mask] %>% rowMeans()
    ) %>%
    noiseRunMed(mass)
  return(signal / noise)
}

sparse_matc_to_sparse_matr <- function(mat) {
  stopifnot(is(mat, "sparse_matc"))
  keys <- mat@data[[1]]
  values <- mat@data[[2]]
  keys_all <- mat@keys
  
  keys_new <- rep(list(NULL), nrow(mat))
  values_new <- rep(list(NULL), nrow(mat))
  keys_all_new <- 1:ncol(mat)
  
  for (i in 1:ncol(mat)) {
    if (i %% 100 == 0) {print(i)}
    if (is.null(keys_all)) {
      idx <- match(keys[[i]], keys_all)
    } else {
      idx <- keys[[i]]
    }
    for (j in 1:length(idx)) {
      values_new[[idx[j]]] <- c(values_new[[idx[j]]], values[[i]][[j]])
      keys_new[[idx[j]]] <- c(keys_new[[idx[j]]], i)
    }
  }
  matter::sparse_mat(
    data = list(keys_new, values_new),
    nrow = nrow(mat),
    ncol = ncol(mat),
    keys = keys_all_new,
    rowMaj = TRUE
  )
}

sparse_matc_to_sparseMatrix <- function(mat) {
  stopifnot(is(mat, "sparse_matc"))
  keys <- mat@data[[1]]
  values <- mat@data[[2]]
  keys_all <- mat@keys
  
  total_values <- sum(sapply(values, length))
  row <- integer(total_values)
  col <- integer(total_values)
  x <- numeric(total_values)
  
  k <- 1
  cache <- list()
  for (i in 1:ncol(mat)) {
    if (i %% 10 == 0) {print(i)}
    if (is.null(keys_all)) {
      idx <- keys[[i]]
    } else {
      idx <- sapply(keys[[i]], function(x) which.min(abs(keys_all - x)))
    }
    for (j in 1:length(idx)) {
      row[[k]] <- idx[[j]]
      col[[k]] <- i
      x[[k]] <- values[[i]][[j]]
      k <- k + 1
    }
  }
  Matrix::sparseMatrix(i = row, j = col, x = x, dims = dim(mat), dimnames = dimnames(mat), giveCsparse = TRUE)
}

normalize_tic_mse <- function(mse) {
  # Normalize spectra of Cardinal::MSImagingExperiment.
  # 
  # Each pixel is normalized such that the TIC for that pixel is nrow(mse),
  # following what Cardinal::normalize(method = "tic") does.
  s <- Cardinal::spectra(mse)
  if (matter::is.sparse(s)) {
    if (is(s, "sparse_matc")) { # column major order sparse matrix
      s@data[[2]] <- lapply(X = s@data[[2]], FUN = function(x) {if (sum(x) > 0) {x / sum(x) * nrow(mse)} else {0}})
      Cardinal::spectra(mse) <- s
    } else {
      stop(str_glue("normalize_tic_mse() not implemented for {class(Cardinal::spectra(mse))} matrices."))
    }
  } else if (is.matrix(s)) {
    Cardinal::spectra(mse) <- scale(s, center = FALSE, scale = colSums(s) / nrow(mse))
  } else {
    stop(str_glue("normalize_tic_mse() not implemented for {class(Cardinal::spectra(mse))} matrices."))
  }
  return(mse)
}

image_mse <- function(mse, mz = NULL, plusminus = 0, fun = colMeans, contrast.enhance = NULL) {
  # Create heatmap of mass intensity using ggplot2::geom_raster()
  # 
  # Args
  # - mse: Cardinal::MSImagingExperiment
  # - mz: numeric. default = NULL
  # - plusminus: numeric. default = 0
  #     mass range in units of m/z
  # - fun: callable. default = colMeans
  #     function to aggregate intensities in the mass range
  # - contrast.enhance: str. default = NULL
  #     "histogram" or "suppression" -- see Cardinal::image() documentation
  # 
  # Returns: ggplot
  # - fixed aspect ratio (coord_fixed())
  # - no grid lines (theme_classic())
  # - viridis color palette (scale_fill_continuous(type = "viridis"))
  if (!is.null(mz)) {
    mass <- mz
    mse <- Cardinal::filter(mse, mz >= mass - plusminus, mz <= mass + plusminus)
  }
  df <- as_tibble(Cardinal::pixelData(mse)) %>% 
    mutate(intensity = fun(Cardinal::spectra(mse)))

  if (contrast.enhance == "histogram") {
    df <- mutate(df, intensity = Cardinal:::contrast.enhance.histogram(intensity))
  } else if (contrast.enhance == "suppression") {
    df <- mutate(df, intensity = Cardinal:::contrast.enhance.histogram(intensity))
  }

  g <- df %>% 
    ggplot(aes(x, y, fill = intensity)) +
    geom_raster() + 
    scale_fill_continuous(type = "viridis") +
    coord_fixed() +
    theme_classic()
  return(g)
}