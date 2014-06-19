#' Read the M-bias.txt file produced from Bismark.
#' 
#' Read the M-bias.txt file produced from Bismark and return as a 
#' \code{\link[base]{data.frame}}.
#' 
#' @param file The name of the M-bias file which the data are to be read from. 
#' 
#' @return A \code{\link{data.frame}}. The columns are:
#' \itemize{
#'  \item \code{read}: \code{SE} if the data are from single-end sequencing; 
#' \code{R1} or \code{R2} if the data are read_1 or read_2 from paired-end 
#' sequencing.
#'  \item \code{read_pos}: the position along the read.
#'  \item \code{methylation_type}: the methylation type. Either "CG", "CHG", 
#' "CHH".
#'  \item \code{M}: the number of methylated base calls.
#'  \item \code{U}: the number of unmethylated base calls.
#'  \item \code{total}: the total number of base calls.
#'  \item \code{perc_M}: the percentage of methylated base calls.
#' }
#' 
#' @note \strong{WARNING}: The "position along the read" is only the same as 
#' the "sequencing cycle" if the reads did not undergo any pre-alignment 
#' trimming. 
#' 
#' @export
#' 
#' @examples
#' cat("TODO")
read.mbias <- function(file) {
  
  # Read the file.
  x <- read.table(file, sep = '\n', as.is = TRUE)
  
  # Find 'header' lines, which are dispersed throughout the file.
  header_lines <- sapply(X = x, function(xx) {
    grepl(pattern = 'context', x = xx)
    })
  y <- vector(mode = "list", length = sum(header_lines))
  names(y) <- x[header_lines]
  
  # z indexes the lines containing M-bias values.
  z <- sapply(X = x, function(xx) {
    grep(pattern = "context|^=|position", x= xx, invert = TRUE)
    })
  dz <- diff(z)
  next_context_idx <- c(which(dz > 1), nrow(z))
  colnames_idx <- sapply(X = x, function(xx) {
    grep(pattern = "position", x = xx)
    })[1]
  colnames <- strsplit(x = x[colnames_idx, ], split = "\t")[[1]]
  
  # Loop over the lines in the file.
  for (i in seq_along(next_context_idx)) {
    # Note the "+4" to skip the intermediate 'header' rows
    start_idx <- ifelse(i == 1, 4, z[next_context_idx[i - 1]] + 4) 
    end_idx <- z[next_context_idx[i]]  
    yy <- x[seq(from = start_idx, to = end_idx, by = 1), ]
    y[[i]] <- do.call("rbind", lapply(yy, function(yyy, colnames) {
      tmp <- strsplit(yyy, split = "\t")
      val <- data.frame(as.integer(tmp[[1]][1]), 
                        as.integer(tmp[[1]][2]), 
                        as.integer(tmp[[1]][3]), 
                        as.numeric(tmp[[1]][4]), 
                        as.integer(tmp[[1]][5]))
      colnames(val) <- colnames
      return(val)
    }, colnames = colnames))
  }
  y <- do.call("rbind", y)
  y$Context <- strtrim(rownames(y), width = 3)
  
  # Add R1/R2 annotation if paired-end data
  if (grepl(pattern = "R1", x = rownames(y)[1])){
    y$Read <- ifelse(grepl(pattern = "R1", x = rownames(y)), "R1", "R2")
  } else{
    y$Read <- "SE"
  }
  rownames(y) <- NULL
  
  # Simplify column names
  colnames(y) <- c('read_pos', 'M', 'U', 'perc_M', 'total', 
                   'methylation_type', 'read')
  
  # Re-order columns
  y <- y[, c('read', 'read_pos', 'methylation_type', 'M', 'U', 'total', 'perc_M')]
  
  # Replace CpG with CG
  y$methylation_type <- ifelse(y$methylation_type == 'CpG', 'CG', 
                               y$methylation_type)
  
  # Convert to an MBias object.
  return(y)
}