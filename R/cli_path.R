#' Get path to the standard bcftools binary
#'
#' Returns the path to the standard bcftools executable included with RBCFLib.
#' This is the same bcftools binary used by the bcftools project, allowing
#' you to avoid spurious warnings and use the standard bcftools interface.
#'
#' @return Character string containing the path to the bcftools executable
#'
#' @examples
#' \dontrun{
#' # Get path to bcftools binary
#' bcftools_bin <- BCFToolsBinaryPath()
#'
#' # Use system2 to run standard bcftools commands
#' system2(bcftools_bin, c("view", "--help"))
#' system2(bcftools_bin, c("stats", "input.vcf"))
#' }
#'
#' @export
BCFToolsBinaryPath <- function() {
  bcftools_path <- system.file("bin", "bcftools", package = "RBCFLib")

  if (bcftools_path == "") {
    stop("bcftools binary not found. Please reinstall RBCFLib.")
  }

  return(bcftools_path)
}

#' Get path to the assoc_plot.R script
#' @return Character string containing the path to the assoc_plot.R script
#'
#' @examples
#' \dontrun{
#' # Get path to assoc_plot.R script
#' assoc_plot_script <- assocPlotRscript()
#'
#'  Run the script on the termninal
#'  requires ggplot2 and data.table
#' }
#'
#' @export
assocPlotRscript <- function() {
  assoc_plot <- system.file("bin", "assoc_plot.R", package = "RBCFLib")
  return(assoc_plot)
}
