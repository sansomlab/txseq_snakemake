#!/usr/bin/env Rscript
# Install the txseq R package and all R dependencies of the txseq pipelines.
#
# All dependencies are declared in txseq/DESCRIPTION:
#   Imports  - packages used by the txseq R package itself
#   Suggests - packages used by the pipeline scripts and reports
#   Remotes  - GitHub sources for packages not on CRAN or Bioconductor
#
# pak resolves CRAN, Bioconductor and GitHub packages in a single step and
# installs only what is missing. Existing packages are not upgraded.
#
# Usage (run on a login node; compute nodes have no internet access):
#   Rscript src/txseq/R/txseq/install.R [--lib=/path/to/library] [--no-suggests]

args <- commandArgs(trailingOnly = TRUE)
lib_arg <- sub("^--lib=", "", grep("^--lib=", args, value = TRUE))
with_suggests <- !("--no-suggests" %in% args)

# Locate the package relative to this script, so it can be run from any directory
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}
# Works whether this script sits inside the package directory or beside it
candidates <- c(script_dir, file.path(script_dir, "txseq"))
pkg_dir <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))][1]
if (is.na(pkg_dir)) {
  stop("Cannot find the txseq package (no DESCRIPTION file) in ",
       paste(candidates, collapse = " or "))
}

# Target library: --lib if given, otherwise the first writable library path
lib <- if (length(lib_arg)) lib_arg[1] else .libPaths()[1]
if (!dir.exists(lib)) dir.create(lib, recursive = TRUE, showWarnings = FALSE)
if (file.access(lib, 2) != 0) {
  lib <- Sys.getenv("R_LIBS_USER")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(lib, .libPaths()))

# Bootstrap pak (self-contained; it brings no further dependencies)
if (!requireNamespace("pak", quietly = TRUE)) {
  message("Installing pak into ", lib)
  install.packages("pak", lib = lib, repos = "https://cloud.r-project.org")
}

# TRUE  = hard dependencies plus Suggests of txseq (not Suggests of its dependencies)
# NA    = hard dependencies only
pak::local_install(
  pkg_dir,
  lib = lib,
  dependencies = if (with_suggests) TRUE else NA,
  upgrade = FALSE,
  ask = FALSE
)

if (!requireNamespace("txseq", quietly = TRUE)) {
  stop("txseq failed to install. See the output above for the cause.")
}
message("txseq ", packageVersion("txseq", lib.loc = lib), " installed to ", lib)
