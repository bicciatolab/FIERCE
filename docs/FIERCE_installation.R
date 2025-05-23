CRANdep <- c("ranger","rlang","pbkrtest","BiocManager")
versions <- c("0.14.1","1.1.0","0.5.1","1.30.19")
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

if(!"qlcMatrix"%in% installed.packages()[,"Package"]){devtools::install_github("cysouw/qlcMatrix")}
devtools::install_github("aet21/SCENT")

devtools::install_github("bicciatolab/FIERCE")

