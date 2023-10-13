CRANdep <- c("ranger","rlang","pbkrtest","BiocManager")
versions <- c("0.14.1","1.0.6","0.5.1","1.30.19")
for (i in 1:length(CRANdep)) {
  if (!(CRANdep[i] %in% installed.packages()[,"Package"]) | !(versions[i] %in% installed.packages()[installed.packages()[,"Package"]==CRANdep[i],"Version"])) {
      devtools::install_version(CRANdep[i], version = versions[i], repos = "https://cloud.r-project.org")
    }
}

BioCdep <- c("GenomeInfoDb","GenomicRanges","SummarizedExperiment","SingleCellExperiment","destiny")
versions <- c("1.26.7","1.42.0","1.20.0", "1.12.0", "3.4.0")
for (i in 1:length(BioCdep)) {
  if (!(BioCdep[i] %in% installed.packages()[,"Package"]) | !(versions[i] %in% installed.packages()[installed.packages()[,"Package"]==BioCdep[i],"Version"])) {
      BiocManager::install(BioCdep[i], version = "3.12")
    }
}

devtools::install_github("aet21/SCENT")

devtools::install_github(
  "bicciatolab/FIERCE",
  ref = "HEAD",
  subdir = NULL,
  auth_token = "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",
  host = "api.github.com",
  dependencies = NA,
  upgrade = c("default", "ask", "always", "never"),
  force = FALSE,
  quiet = FALSE,
  build = TRUE,
  build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"),
  build_manual = FALSE,
  build_vignettes = FALSE,
  repos = getOption("repos"),
  type = getOption("pkgType")
)