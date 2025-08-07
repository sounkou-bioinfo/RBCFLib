#' Get path to BCFTools CLI executables
#'
#' Returns the path to the command-line interface executables provided by RBCFLib.
#' These CLI tools allow users to run bcftools commands from the command line
#' using the RBCFLib package's BCFToolsRun function.
#'
#' @param cli Character string; which CLI to get the path for.
#'   Either "RBcftools" or "BCFToolsCli". Default is "RBcftools".
#'
#' @return Character string containing the path to the requested CLI executable
#'
#' @examples
#' \dontrun{
#' # Get path to RBcftools
#' bcftools_cli <- BCFToolsCLIPath()
#'
#' # Get path to BCFToolsCli (PascalCase version)
#' bcftools_pascal_cli <- BCFToolsCLIPath("BCFToolsCli")
#'
#' # Use system2 to run commands
#' system2(bcftools_cli, c("view", "-h", vcfFile))
#' }
#'
#' @export
BCFToolsCLIPath <- function(cli = c("RBcftools", "BCFToolsCli")) {
  cli <- match.arg(cli)
  cli_path <- system.file("bin", cli, package = "RBCFLib")

  if (cli_path == "") {
    stop("CLI executable not found. Please reinstall RBCFLib.")
  }

  return(cli_path)
}
