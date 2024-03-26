#' Pancreas endocrinogenesis test dataset
#'
#' A subset of 499 cells from the mouse pancreas endocrinogenesis dataset sequenced by Bastidas-Ponce et al (2019)
#' This subset was drawn from a bigger dataset that was sequenced by Bastidas-Ponce et al (2019) to study mouse pancreas endocrinogenesis, a developmental process that leads to the differentiation of endocrine cells from multiple pluripotent cell subpopulations. In particular, this subset (included in the scVelo python package) includes the differentiation lineage of endocrine cells from ductal cells, passing through endocrine progenitors (EP) and pre-endorcine cells. The anndata object includes spliced and unspliced counts in separate layers, and pre-computed PC and UMAP coordinates
#'
#' @format ## `pancreas`
#' An AnnData object with 499 rows and 27998 columns:
#' \describe{
#'   \item{country}{Germany}
#'   \item{iso2, iso3}{2 & 3 letter ISO country codes}
#'   \item{year}{2019}
#'   ...
#' }
#' @source <https://github.com/theislab/scvelo>
"pancreas"


#' Dentate gyrus neurogenesis test dataset
#'
#' A subset of 499 cells from the mouse dentate gyrus neurogenesis dataset sequenced by Hochgerner et al (2018)
#' This subset was drawn from a bigger dataset that was sequenced by Hochgerner et al (2018) to study mouse dentate gyrus neurogenesis, a developmental process that leads to the differentiation of neural cells from radial glia and intermediate neuronal progenitor cells (nIPCs). In particular, this subset (included in the scVelo python package) includes the differentiation lineage of granule cells from nIPCs, passing through neuroblasts. The anndata object includes spliced and unspliced counts in separate layers, and pre-computed tSNE coordinates
#'
#' @format ## `dentate gyrus`
#' An AnnData object with 499 rows and 27998 columns:
#' \describe{
#'   \item{country}{Sweden}
#'   \item{iso2, iso3}{2 & 3 letter ISO country codes}
#'   \item{year}{2018}
#'   ...
#' }
#' @source <http://pklab.med.harvard.edu/velocyto/DentateGyrus/DentateGyrus.loom>
"dentate_gyrus"
