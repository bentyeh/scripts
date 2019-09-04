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
  return_function = FALSE,
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
  #     Built-in options: "loess", "lm", "spline", "quad", "auto"
  #       - "auto" uses "lm" if number of matched peaks
  #     To use a custom model, pass in a function that takes a data.frame (x = peaks, y = ref_masses)
  #       as an argument and returns a vectorized function mapping original to calibrated masses.
  # - return_function: logical. default = FALSE
  #     Return calibration function, as opposed to returning calibrated spectrum
  # - verbose: logical. default = TRUE
  # 
  # Returns: matrix or function
  #   If `return_function` == TRUE: returns a function that takes a spectrum (2-column matrix)
  #     and returns a calibrated spectrum.
  #   Otherwise, returns a calibrated spectrum (2-column matrix).
  # 
  # Notes
  # - Arguments inspired by AB/SCIEX 5800 Series Explorer's Processing Method parameters.
  # - Calibration algorithm based on Cardinal::mzAlign(). See
  #   https://rdrr.io/bioc/Cardinal/src/R/process2-mzAlign.R

  if (is.null(min_peaks)) {
    min_peaks <- length(ref_masses)
  }
  if (identical(noise_method, "runmed")) {
    noise_method <- noiseRunMed
  }
  stopifnot(tol_units %in% c("mz", "ppm"))
  stopifnot(min_peaks <= length(ref_masses))
  stopifnot(length(tol) == length(tol_units))
  stopifnot(length(tol) %in% c(1, length(ref_masses)))
  stopifnot(rlang::is_callable(noise_method) || noise_method %in% c("simple", "adaptive", "mad"))
  stopifnot(rlang::is_callable(model) || model %in% c("loess", "lm", "quad", "spline", "auto"))

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
  if (rlang::is_callable(noise_method)) {
    for (i in 1:n) {
      noise[i] <- noise_method(spectrum, peaks[i, 1])
    }
  } else {
    noise_method <- list(
      simple = Cardinal:::.estimateNoiseSimple,
      adaptive = Cardinal:::.estimateNoiseAdaptive,
      mad = Cardinal:::.estimateNoiseMAD
    )[[noise_method]]
    noise_all <- noise_method(spectrum[, 2])
    for (i in 1:n) {
      noise[i] <- noise_all[spectrum[,1] == peaks[i, 1]]
    }
  }

  snr <- peaks[, 2] / noise
  detected_peaks <- snr >= SNR
  n_detected_peaks <- sum(detected_peaks)
  if (verbose) {
    print(stringr::str_glue("Detected peaks ({n_detected_peaks}):"))
    print(data.frame(ref = ref_masses, detected_mass = peaks[,1], snr = snr, passed = detected_peaks))
  }
  if (n_detected_peaks < min_peaks) {
    message(stringr::str_glue(
      "Number of detected peaks ({n_detected_peaks}) < min_peaks ",
      "({min_peaks}). Not calibrating."
    ))
    if (verbose) {
      cat("Using identity model.\n")
    }
    model <- function(x) {return(identity)}
  }

  if (identical(model, "auto")) {
    model <- ifelse(n_detected_peaks > 3, "spline", "lm")
    if (verbose) {
      print(stringr::str_glue("Using {model} model."))
    }
  }

  ref_masses <- ref_masses[detected_peaks]
  peaks <- peaks[detected_peaks, , drop = FALSE]

  df <- data.frame(x = peaks[, 1], y = ref_masses)
  if (is.character(model) && model %in% c("loess", "lm", "quad")) {
    model <- list(
      loess = loess(y ~ x, data = df, span = 1, control = loess.control(surface = "direct")),
      lm = lm(y ~ x, data = df),
      quad = lm(y ~ poly(x, 2), data = df)
    )[[model]]
    warp <- function(x) pmax(predict(model, data.frame(x = x)), 0)
  } else if (identical(model, "spline")) {
    warp <- splinefun(peaks[, 1], ref_masses, "natural")
  } else {
    warp <- model(df)
  }

  if (verbose) {
    residuals <- peaks[, 1] - ref_masses
    if (length(residuals) == 0) {
      residuals <- NA
    }
    rss <- sum(residuals ^ 2)
    cat("Residuals (original):\n")
    cat(residuals)
    cat("\n\n")
    print(stringr::str_glue("RSS (original): {rss}"))

    residuals <- warp(peaks[,1]) - ref_masses
    if (length(residuals) == 0) {
      residuals <- NA
    }
    if (!is.character(model)) {
      tmp <- try(print(summary(model)), silent = TRUE)
    } 
    if (is.character(model) || identical(class(tmp), "try-error")) {
      cat("Residuals (fitted):\n")
      cat(residuals)
      cat("\n\n")
    }
    rss <- sum(residuals ^ 2)
    print(stringr::str_glue("RSS (fitted): {rss}"))
  }

  if (return_function) {
    return(warp)
  } else {
    spectrum[, 1] <- warp(spectrum[, 1])
    return(spectrum)
  }
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
                         show_noise = FALSE, noise_col = "red", plotter = 'baseR', ...) {
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
  # - plotter: character. default = 'baseR'
  #     'baseR': create base R plot
  #     'ggplot': create ggplot2 plot
  # - ...
  #     Arguments to pass to base R plotter
  # 
  # Returns: NULL (`plotter` == 'baseR') or ggplot (`plotter` == 'ggplot')

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

  if (plotter == 'baseR') {
    do.call(plot, args)
    if (is.numeric(show_noise)) {
      abline(h = show_noise, col = noise_col)
    }
  } else {
    colnames(spectrum) <- c("m/z", "count")
    g <-
      spectrum %>%
      as_tibble() %>%
      ggplot(aes(`m/z`, count)) +
      geom_col() +
      labs(x = "m/z")
    if (is.numeric(show_noise)) {
      g <- 
        g +
        geom_hline(yintercept = show_noise, color = "red")
    }
    print(g)
    return(g)
  }
}