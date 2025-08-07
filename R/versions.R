#' Get htslib Version
#'
#' Returns the version of the htslib library used by RBCFLib
#'
#' @return A character string with the htslib version
#' @export
HTSLibVersion <- function() {
  .Call(RC_HTSLibVersion)
}

#' Get bcftools Version
#'
#' Returns the version of the bcftools library used by RBCFLib
#'
#' @return A character string with the bcftools version
#' @export
BCFToolsVersion <- function() {
  .Call(RC_BCFToolsVersion)
}
