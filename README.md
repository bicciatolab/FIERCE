# FIERCE

#### Contact:

mattia.forcato@unimore.it; luca.calderoni@unimore.it

#### Citation:

L. Calderoni, F. Grandi, S. Bicciato, O. Romano and M. Forcato

# Table of Contents

- [System requirements](https://github.com/bicciatolab/FIERCE#system-requirements)
- [Installation Python packages](https://github.com/bicciatolab/FIERCE#install-python-packages)
- [Installation in R](https://github.com/bicciatolab/FIERCE#install-packages-in-r)
- [Tutorial on example data](https://github.com/bicciatolab/FIERCE/main/docs/FIERCE_tutorial.html)

## System requirements

* R version: >= 4.0.5
* Dependencies: *reticulate*, *SCENT*, *Seurat*, *SeuratObject*, *ggplot2*, *reshape2*, *scales*, *usethis*, *tools*, *usethis*.

* Python version: >= 3.10
* Dependencies: *numpy*, *scipy*, *cython*, *numba*, *matplotlib*, *scikit-learn*, *h5py*, *click*, *velocyto*, *scvelo*,*pandas*.

## Installation

To install `FIERCE` package, especially if working on a Linux machine, we strongly recommend the use of a virtual Anaconda environment to avoid conflicts between package's dependencies versions.

If Anaconda is already installed on your system, a virtual environment for `FIERCE` can be set either manually, installing all packages one by one, or automatically following the instruction provided below.

If you are not working on a Linux machine, all required packages must be manually installed on your R and Python.

#### Create a `FIERCE` environment and install all anaconda available packages automatically

On a Linux system, a working Anaconda environment for `FIERCE` can be set exploiting two different strategies:

1- Generating the `FIERCE` environment from a .yml file

According to this installation procedure, a provided pre-configured working environment can be downloaded from your terminal using the command

```bash
wget https://github.com/bicciatolab/FIERCE/blob/main/docs/FIERCE.yml
```

and then cloned with 

```bash
conda env create -n FIERCE --file FIERCE.yml
```
or 

2- Creating the `FIERCE` environment using a single pre-defined Anaconda command

To install the environment open your terminal and run:

```bash
conda create --name FIERCE -c conda-forge r-base=4.0.5 scanpy=1.9.1 python-louvain=0.15 tqdm=4.64.1 pandas=1.5.1 scipy=1.9.3 numba=0.56.3 matplotlib-base=3.6.2 h5py=3.7.0 click=8.1.3 r-ggplot2=3.3.6 r-reshape2=1.4.4 r-scales=1.2.1 anndata=0.8.0 r-rgeos=3.11.0 r-igraph=1.3.4 python-igraph=0.10.2 r-leiden=0.4.3 r-rcurl=1.98_1.8 r-devtools=2.4.4 r-reticulate=1.15
```

Whatever the strategy applied, once created the environment, it may be accessed through the command:

```bash
conda activate FIERCE
```

and the last dependencies must be installed as explained below

#### Install Python packages

accessed the environment use the following lines of code on your terminal to install the remaining mandatory python packages:

```python
python -m pip install cython
#pip3 install Cython==0.29.32 
pip3 install velocyto==0.17.17 scvelo scikit-learn
```

Since not all required R packages are provided in anaconda.org, some packages must be installed directly from the console, as described in [install packages in R](https://github.com/bicciatolab/FIERCE#install-packages-in-r).

#### Install packages in R

Thus, before installing `FIERCE`, users must open their R console and run the following code to install package's dependencies, which are deployed by CRAN, Bioconductor, and Github.

```r
CRANdep <- c("spatstat.core","ranger","pbkrtest","Seurat", "SeuratObject", "BiocManager")
versions <- c("2.4-4","0.14.1","0.5.1","4.2.0","4.1.2","1.30.19")
for (i in 1:length(CRANdep)) {
  if (!(CRANdep[i] %in% installed.packages()[,"Package"]) | !(versions[i] %in% installed.packages()[installed.packages()[,"Package"]==CRANdep[i],"Version"])) {
      devtools::install_version(CRANdep[i], version = versions[i], repos = "https://cloud.r-project.org")
    }
}
```

```r
BioCdep <- c("GenomeInfoDb","GenomicRanges","SummarizedExperiment","SingleCellExperiment","destiny")
versions <- c("1.26.7","1.42.0","1.20.0", "1.12.0", "3.4.0")
for (i in 1:length(BioCdep)) {
  if (!(BioCdep[i] %in% installed.packages()[,"Package"]) | !(versions[i] %in% installed.packages()[installed.packages()[,"Package"]==BioCdep[i],"Version"])) {
      BiocManager::install(BioCdep[i], version = "3.12")
    }
}
```

```r
if(!"SCENT"%in% installed.packages()[,"Package"]){devtools::install_github("aet21/SCENT")}

```

Once installed all dependencies, `FIERCE` can be installed with the following script:

```r
devtools::install_github("bicciatolab/FIERCE")
```

In case of any issue with `FIERCE` installation through `install_github`, it is possible to download the package.tar.gz and install `FIERCE` from a local repository using the following script:

```r
install.packages("/path/to/package_directory", repos = NULL, type="source")
```

#### Automatic installation

More bravely, `FIERCE` can be automatically installed downloading its installation file by

```bash
wget https://github.com/bicciatolab/FIERCE/blob/main/docs/FIERCE_istallation.sh
```

and running it with

```bash
sh /path/to/Installation/file/FIERCE_istallation.sh
```
or

```bash
source /path/to/Installation/file/FIERCE_istallation.sh
```
