% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIERCE.R
\name{plot_velocity}
\alias{plot_velocity}
\title{plot_velocity}
\usage{
plot_velocity(
  adata,
  project_dir = "./FIERCE_results",
  sqrt_transform = TRUE,
  embedding_basis = "umap",
  force_graph_recalc = FALSE,
  compute_latent_time = FALSE,
  color_as = NULL,
  lab_order = NULL,
  palette = NULL,
  legend_loc = "right margin",
  alpha = 0.3,
  add_outline = TRUE,
  min_mass = 4,
  n_cores = NULL,
  adata_copy = FALSE
)
}
\arguments{
\item{adata}{anndata object containing the results of the "compute_velocity" function}

\item{project_dir}{name of the directory containing the results of the main FIERCE analysis (including the path). If it does not exist, it will be created. The default name is "./FIERCE_results". The streamplot will be saved in the "velocity_field_streamplots" sub-directory}

\item{sqrt_transform}{boolean; whether to apply the variance-stabilizing transformation during transition probabilities computation. It helps to obtain a smoother streamplot. Default is TRUE}

\item{embedding_basis}{name of the embedding to use for the streamplot. Default is "umap"}

\item{force_graph_recalc}{whether to force the recalculation of the cell-cell transition probabilities matrix, if already present. Default is FALSE}

\item{compute_latent_time}{whether to additionally compute and plot the latent time of each cell in the dynamical model. Default is FALSE}

\item{color_as}{character vector containing the names of the cell annotations in the "obs" slot to visualize on the streamplot}

\item{lab_order}{optional; list containing, for each cell annotation specified in "color_as", a character vector specifying the order by which the respective labels should be printed in the legends of the streamplots. If it is not necessary to specify a particular order for a specific annotation, just write NULL. If it is not necessary to specify any order for any annotation, do not change the default value of this parameter}

\item{palette}{optional; list containing, for each cell annotation specified in "color_as", a character vector specifying the colors of the respective labels. If it is not necessary to specify any color for a specific annotation, just write NULL; in this case, a default color palette will be used. If it is not necessary to specify any color for any annotation, do not change the default value of this parameter; in this case, annotations will be colored according either to a default palette, or to the colors that are already stored in the anndata object ("uns" slot)}

\item{legend_loc}{localization of the legend on the streamplot. Default is "right margin"}

\item{alpha}{transparency of the dots on the streamplots, ranging from 0 (fully transparent) to 1 (fully opaque). Default is 0.3}

\item{add_outline}{whether to add an outline around groups of dots on the streamplots. Default is TRUE}

\item{min_mass}{minimum magnitude required for velocity vectors to be drawn on the streamplot, ranging from 0 (all vectors) to 5 (most intense vectors only). Default is 4}

\item{n_cores}{number of cores to use to compute the cell-cell transition probability matrix (most computationally expensive step). If NULL (default), only 1 will be used}

\item{adata_copy}{boolean; if TRUE, create a new anndata object with the results of the analysis; if FALSE, update the provided anndata object. Default is FALSE}
}
\value{
If adata_copy=TRUE, a new anndata object with the results of the analysis will be returned. Otherwise (default), the provided anndata object will be updated with the results of the analysis. In both cases, the velocity graph (correlations used for transition probabilities computation) will be saved in the "uns" slot, and the predicted displacements of cells (used for streamplot drawing) will be saved in the "obsm" slot
}
\description{
Computes the RNA velocity vector field and draws the streamplot on an embedding of choice
}
\examples{
#Draw the streamplot on the UMAP embedding. In this example, the clusters, the cell types and the cell cycle phase are visualized on the streamplot. A specific label order and a specific color for each label is specified for the cell cycle phase only

plot_velocity(adata, color_as=c("clusters", "cell_types", "phase"), lab_order=list(NULL, NULL, c("G1","S","G2M")), palette=list(NULL, NULL, c("green", "red", "blue")))

}
