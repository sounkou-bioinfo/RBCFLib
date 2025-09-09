#' Get htslib Version
#'
#' Returns the version of the htslib library used by RBCFLib
#'
#' @return A character string with the htslib version
#' @export
HTSLibVersion <- function() {
  .Call(RC_HTSLibVersion, PACKAGE = "RBCFLib")
}

#' Get bcftools Version
#'
#' Returns the version of the bcftools library used by RBCFLib
#'
#' @return A character string with the bcftools version
#' @export
BCFToolsVersion <- function() {
  .Call(RC_BCFToolsVersion, PACKAGE = "RBCFLib")
}


#' Get bcftools score Version
#'
#' Returns the version of the bcftools score plugin used by RBCFLib
#''
#' @return A character string with the bcftools score version
#' @export
BCFToolsScoreVersion <- function() {
  .Call(RC_BCFToolsScoreVersion, PACKAGE = "RBCFLib")
}
