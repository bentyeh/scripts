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
  # - min: numeric. default = NULL
  #     Minimum m/z value to plot
  # - max: numeric. default = NULL
  #     Maximum m/z value to plot
  # - mass: numeric. default = NULL
  #     Mass of peak of interest.
  # - tol_ppm: numeric. default=1e4
  #     Mass tolerance "window."
  # - plotter: character. default = 'baseR'
  #     'baseR': create base R plot
  #     'ggplot': create ggplot2 plot
  # - ...
  #     Arguments to pass to base R plotter
  # 
  # Returns: NULL (`plotter` == 'baseR') or ggplot (`plotter` == 'ggplot')
  # 
  # Notes: The arguments min and max take precedence over mass and tol_ppm.
  
  stopifnot(plotter %in% c('baseR', 'ggplot'))
  xlim <- c(min(spectrum[,1]), max(spectrum[,1]))
  
  if (!is.null(min)) {
    spectrum <- spectrum[spectrum[,1] >= min,]
    xlim[1] <- min
  }
  
  if (!is.null(max)) {
    spectrum <- spectrum[spectrum[,1] <= max,]
    xlim[2] <- max
  }
  
  if (is.null(min) && is.null(max) && !is.null(mass)) {
    tol <- mass * tol_ppm/1e6
    idx <- which((spectrum[,1] > mass - tol) &
                 (spectrum[,1] < mass + tol))
    spectrum <- spectrum[idx,]
    xlim <- c(min(spectrum[,1]), max(spectrum[,1]))
  }
  
  if (plotter == 'baseR') {
    plot(spectrum, type = 'h', xlab = 'm/z', ylab = 'count',
         ylim = c(0, max(spectrum[,2])), yaxs = 'i',
         xlim = xlim, xaxs = 'i', bty = 'l', ...)
  } else {
    colnames(spectrum) <- c('m/z', 'count')
    spectrum %>%
      as_tibble() %>%
      ggplot(aes(`m/z`, count)) +
      geom_col() +
      labs(x = 'm/z')
  }
}