quiet <- function(x) {
  # Evaluate an object or expression while suppressing output
  # Credit: Hadley Wickham
  # Source: https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Wait on child processes
tryCatch(
  {
    wait <- Rcpp::cppFunction(
      "void wait() {int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};}",
      includes = "#include <sys/wait.h>"
    )
  },
  error = function(e) {
    e
    message("Rcpp may not be installed.")
  },
  NULL
)

arrayToDf <- function(ar) {
  # Convert an array to a data frame.
  # 
  # Args
  # - ar: array
  # 
  # Returns: data.frame. dim = [prod(dim(ar)), length(dim(ar)) + 1]
  # - Rows: one row per element of the array
  #   - Rows are "entered" into the data frame in column major order
  # - Columns: one column per dimension of the array, plus one column for the values in the array
  #   - Column names are retained from names of array dimensions (`names(dimnames(ar)`).
  #     Where non-existant, column names are given as 'd#', where # is the number of the dimension.
  #   - For each column, if the corresponding array dimension was named, then the column is
  #     of class 'character'. Otherwise, the column is of 'numeric' class.
  #   - The last column (the values from the array) is named 'value'.
  # 
  # References
  # - as.data.frame.table()
  # - https://stackoverflow.com/a/42810479

  dims = dim(ar)
  ndims = length(dims)

  if (is.null(dimnames(ar))) {
    nullDims = 1:ndims
    nullDimNames = 1:ndims
  } else {
    nullDims = which(sapply(dimnames(ar), is.null))
    if (is.null(names(dimnames(ar)))) {
      nullDimNames = 1:ndims
    } else {
      nullDimNames = which(sapply(names(dimnames(ar)), function(x) identical(x, "")))
    }
  }
  namedDims = setdiff(1:ndims, nullDims)

  df = as.data.frame.table(ar)
  df[namedDims] = lapply(df[namedDims], as.character)
  df[nullDims] = lapply(df[nullDims], as.numeric)

  colnames(df)[nullDimNames] = paste0('d', nullDimNames)
  colnames(df)[ncol(df)] = 'value'

  return(df)
}

dfToArray <- function(df, dimOrders = NULL) {
  # Convert a data frame to an array.
  # 
  # Args
  # - df: data.frame. dim=[nr, nc]
  #     The first (nc-1) columns represent dimensions. The last column gives values in the array.
  # - dimOrders: list of vector. default = NULL
  #     Ordering (i.e., factor levels) of each dimension in the array.
  #     Must be a fully named list (arbitrary length) or fully unnamed list (length = nc).
  #       If a column ordering is given for a column, any values present in that column but
  #       not in the ordering assumes a value of NA.
  #     Any columns for which an ordering is not given is sorted lexicographically.
  #     The class of each vector should match the class of the corresponding column.
  # 
  # Returns: array
  # 
  # References
  # - https://stackoverflow.com/a/9617424
  # - https://stackoverflow.com/a/46129338

  nDim = ncol(df) - 1
  stopifnot(is.null(dimOrders) || !is.null(names(dimOrders)) || length(dimOrders) == nDim)

  df = unique(df)

  # convert columns in df to factors
  if (is.null(dimOrders) || !is.null(names(dimOrders))) {
    namedCols = intersect(names(df[1:nDim]), names(dimOrders))
    unnamedCols = setdiff(names(df[1:nDim]), namedCols)
    for (c in namedCols) {
      df[[c]] = factor(df[[c]], levels = dimOrders[[c]])
    }
    for (c in unnamedCols) {
      df[[c]] = factor(df[[c]])
    }
  } else {
    for (c in 1:nDim) {
      df[[c]] = factor(df[[c]], levels = dimOrders[[c]])
    }
  }

  # initialize array
  ar = array(dim = sapply(df[1:nDim], function(c) length(levels(c))),
             dimnames = sapply(df[1:nDim], levels))

  # input values into array
  ar[do.call(cbind, df[1:nDim])] = df[[nDim + 1]]

  return(ar)
}

orderArray <- function(ar, dims = NULL, metric = 'cor') {
  # Order each dimension of array by hierarchical clustering of pairwise correlations.
  # 
  # Args
  # - ar: array
  # - dims: vector, integer. default=NULL.
  #     Dimensions to order. If NULL, all dims are ordered.
  # - metric: character. default='cor'
  #     Metric by which to compare hyperplanes along dimension. Hyperplanes are flattened
  #     into vectors for comparison.
  #     - 'cor': Pairwise correlation
  #         Any NA values in the pairwise correlation matrices (e.g. too many missing entries or
  #         zero variance) are replaced with the maximum distance value of 1.
  #     - 'dist': Euclidean distance
  # 
  # Returns: list, vector, integer. length = length(dim(ar))
  #   Each vector gives the permutation of the corresponding array dimension.

  nDims = length(dim(ar))
  if (is.null(dims)) {
    dims = 1:nDims
  }

  # validate that dims to order are present in array
  stopifnot(all(dims %in% 1:nDims))

  # validate metric
  stopifnot(metric %in% c('cor', 'dist'))

  # create named dimOrder list
  dimOrder = vector('list', nDims)
  if (!is.null(names(dimnames(ar)))) {
    names(dimOrder) = sapply(names(dimnames(ar)), function(x) ifelse(identical(x, ''), NULL, x))
  }

  # order specified dims
  for (d in dims) {
    # construct matrix to pass into dist()
    # - each row represents a hyperplane along dimension d
    nValues = dim(ar)[d]
    mat = matrix(nrow = nValues, ncol = prod(dim(ar)[-d]))

    # flatten each hyperplane in dimension d into a row vector
    for (value in 1:nValues) {
      mat[value, ] = as.vector(indexArray(ar, d, value))
    }

    if (metric == 'cor') {
      distMat = (1 - cor(t(mat), use = "pairwise.complete.obs"))/2
      distMat[is.na(distMat)] = 1
      distMat = as.dist(distMat)
    } else { # metric == 'dist'
      distMat = dist(mat)
    }
    dimOrder[[d]] = hclust(distMat)$order
    # alternatively, order.dendrogram(as.dendrogram(hclust(distMat)))
  }

  # maintain original orders for other dims
  unOrderedDims = setdiff(1:nDims, dims)
  dimOrder[unOrderedDims] = lapply(unOrderedDims, function(x) 1:(dim(ar)[x]))

  return(dimOrder)
}

indexArray <- function(x, dim, value, drop = FALSE, verbose = FALSE) {
  # Select along one dimension in multidimensional array
  # 
  # Examples: Let A = array(data=1:24, dim=c(2, 3, 4), dimnames=list(c('x1', 'x2'), c('y1', 'y2'), NULL)).
  # - indexArray(A, 1, 1) == A[1, , ]
  # - indexArray(A, 1, 2) == A[2, , ]
  # - indexArray(A, 2, 3) == A[, 3, ]
  # 
  # Args
  # - x: array
  # - dim: integer
  #     Dimension along which to select one hyperplane
  # - value: integer
  #     Hyperplane to select from given dimension
  # - drop: logical. default=FALSE
  #     For matrices and arrays. If TRUE the result is coerced to the lowest possible dimension.
  # - verbose: logical. default=FALSE
  #     Print out equivalent command using `[` operator.
  # 
  # Returns: array or vector
  # 
  # Source: https://stackoverflow.com/a/14502298

  # Create list representing arguments supplied to [
  # bquote() creates an object corresponding to a missing argument
  indices <- rep(list(bquote()), length(dim(x)))
  indices[[dim]] <- value

  # Generate the call to [
  call <- as.call(c(
    list(as.name("["), quote(x)),
    indices,
    list(drop = drop)))

  if (verbose) {
    print(call)
  }

  eval(call)
}