% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIERCE.R
\name{compute_entropy_UMAP}
\alias{compute_entropy_UMAP}
\title{compute_entropy_UMAP}
\usage{
compute_entropy_UMAP(
  adata,
  project_dir = "./FIERCE_results",
  n_neighbors = 30,
  n_pcs = NULL,
  plot_PCA_heatmap = TRUE,
  perform_clustering = FALSE,
  color_as = NULL,
  lab_order = NULL,
  palette = NULL,
  legend_loc = "right margin",
  alpha = 1,
  add_outline = FALSE,
  redo_from_scratch = FALSE,
  adata_copy = FALSE
)
}
\arguments{
\item{adata}{anndata object with "partial_entropies_observed" layer (produced by the "compute_signaling_entropy" function)}

\item{project_dir}{name of the directory containing the results of the main FIERCE analysis (including the path). If it does not exist, it will be created. The default name is "./FIERCE_results". The UMAP plots will be saved in the "UMAP_from_entropy" sub-directory, in separate folders for each combination of PCs and neighbors specified by the user}

\item{n_neighbors}{number of nearest neighbors to compute for each cell in the nearest neighbor graph that will be used for UMAP embedding computation. If the user wishes to try different values for this parameter, an integer vector containing all the desired values can be specified (the results will be saved into dedicated folders). Default is 30 neighbors}

\item{n_pcs}{number of PCs to consider for the nearest neighbor graph computation. If the user wishes to try different values for this parameter, an integer vector containing all the desired values can be specified (the results will be saved into dedicated folders). Default (NULL) is all the computed PCs (50). The PCs are computed directly from the observed entropies with the default PCA procedure of Scanpy}

\item{plot_PCA_heatmap}{boolean; whether to plot the heatmap of the entropies of the most significant genes for each PC. Default is TRUE}

\item{perform_clustering}{boolean; whether to perform clustering analysis on entropy values ("Leiden" algorithm). Default is FALSE. WARNING: a nearest neighbor graph must be computed first}

\item{color_as}{character vector containing the names of the cell annotations in the "obs" slot to visualize on the PCA and UMAP plots. If present, the results of the clustering analysis will be visualized by default}

\item{lab_order}{optional; list containing, for each cell annotation specified in "color_as", a character vector specifying the order by which the respective labels should be printed in the legends of the PCA and UMAP plots. If it is not necessary to specify a particular order for a specific annotation, just write NULL. If it is not necessary to specify any order for any annotation, do not change the default value of this parameter}

\item{palette}{optional; list containing, for each cell annotation specified in "color_as", a character vector specifying the colors of the respective labels. If it is not necessary to specify any color for a specific annotation, just write NULL; in this case, a default color palette will be used. If it is not necessary to specify any color for any annotation, do not change the default value of this parameter; in this case, annotations will be colored according either to a default palette, or to the colors that are already stored in the anndata object ("uns" slot)}

\item{legend_loc}{localization of the legend on the PCA and UMAP plots. Default is "right margin"}

\item{alpha}{transparency of the dots on the PCA and UMAP plots, ranging from 0 (fully transparent) to 1 (fully opaque). Default is 1}

\item{add_outline}{whether to add an outline around groups of dots on the PCA and UMAP plots. Default is FALSE}

\item{redo_from_scratch}{boolean, set it to TRUE to re-compute the PCs and UMAP coordinates and re-draw all the plots from scratch (otherwise, all the steps that have been already performed will be skipped if the function is launched a second time)}

\item{adata_copy}{boolean; if TRUE, create a new anndata object with the results of the analysis; if FALSE, update the provided anndata object. Default is FALSE}
}
\value{
If adata_copy=TRUE, a new anndata object with the results of the analysis will be returned. Otherwise (default), the provided anndata object will be updated with the results of the analysis. In both cases, the PCs and the UMAP coordinates computed from observed entropies will be saved in the "obsm" slot. The nearest neighbor graph (the last that has been computed) will be saved in the "uns" and "obsp" slots
}
\description{
Computes UMAP embeddings from observed entropies with Scanpy. First, a PCA is performed directly on entropy values, then a nearest neighbor graph is computed from resulting PCs, and finally the UMAP embedding is built on the neighbor graph. Optionally, a clustering analysis is also performed on the neighbor graph. The user can specify any combination of PCs and neighbors for UMAP computation, the results will be saved into separate folders
}
\examples{
#Compute UMAP embeddings from observed entropies with default parameters. Separate UMAP embeddings are computed with 10, 30 and 100 neighbors. Additionally, the neighbor graph is computed multiple times considering 15, 30 and 50 PCs. In this example, the clusters, the cell types and the cell cycle phase are visualized on the PCA and UMAP plots. A specific label order and a specific color for each label is specified for the cell cycle phase only

compute_entropy_UMAP(adata, n_pcs=c(15,30,50), n_neighbors=c(10,30,100), color_as=c("clusters", "cell_types", "phase"), lab_order=list(NULL, NULL, c("G1","S","G2M")), palette=list(NULL, NULL, c("green", "red", "blue")), adata_copy=FALSE)

}
