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
  
  tol = mass * tol_ppm/1e6
  idx = which((spectrum[,1] > mass - tol) &
              (spectrum[,1] < mass + tol))
  mass_peak_idx = idx[which.max(spectrum[idx, 2])]
  spectrum[mass_peak_idx,]
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
  # - summary: function. default=median
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
    width = width + 1
  }
  summary(runmed(spectrum[idx, 2], width))
}

sliding_distance <- function(v1, v2, norm_fun = "1", range = NULL, verbose = FALSE) {
  # Args
  # - v1: vector of numeric
  # - v2: vector of numeric
  # - norm: character or callable. default = "1"
  #     If character, specifies the norm type as used by base::norm()
  #     If callable, a function used to compute a vector norm.
  # - range: vector of integer, length = 2. default = NULL
  # - verbose: logical, default = FALSE

  n = length(v1)
  stopifnot(n == length(v2))
  if (is.null(range)) {
    range = c(-n + 1, n - 1)
  }
  stopifnot(length(range) == 2)
  stopifnot(all(abs(range) < n))
  stopifnot(all(range %% 1 == 0))
  stopifnot(range[2] >= range[1])
  
  if (is.character(norm_fun)) {
    type = norm_fun
    norm_fun = function(x) base::norm(as.matrix(x), type = type) / length(x)
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
    v1_aligned = v1[max(1, offset + 1):min(n + offset, n)]
    v2_aligned = v2[max(1, 1 - offset):min(n, n - offset)]
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
  # Args
  # - spectra_list: list of 2-column matrices
  #     First spectrum in the list is assumed to be the "reference" spectrum
  # - trim: logical. default = TRUE
  #     Keep only the region where all spectra are aligned. Otherwise, fill in unaligned regions
  #     with NA.
  # - ...: arguments to pass to align_vectors()
  
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
                         plotter = 'baseR', ...) {
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
  
  if (plotter == 'baseR') {
    do.call(plot, args)
  } else {
    colnames(spectrum) <- c("m/z", "count")
    spectrum %>%
      as_tibble() %>%
      ggplot(aes(`m/z`, count)) +
      geom_col() +
      labs(x = "m/z")
  }
}