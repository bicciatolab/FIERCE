% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIERCE.R
\name{save_h5ad}
\alias{save_h5ad}
\title{save_h5ad}
\usage{
save_h5ad(adata, file_name)
}
\arguments{
\item{adata}{anndata object}

\item{file_name}{name of h5ad output file (including its path)}
}
\value{
No objects are returned. The h5ad file containing the anndata object will be saved in the destination specified in the file name
}
\description{
Saves anndata object to h5ad file
}
\examples{
save_h5ad(adata, "./anndata_object.h5ad")

}
