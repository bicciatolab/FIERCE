#set your working directory  
#cd path/to/your/working/directory 

#R-4.0.5

library(FIERCE)
library(reticulate)

#use_python("/home/luca/.conda/envs/FIERCE/bin/python3.10")

use_condaenv("FIERCE_test")

source_python(system.file("inter_environment_utils.py", package = "FIERCE", lib.loc= .libPaths()))

adata <- load_test_dataset("dentate_gyrus")

load(system.file("mouse.cc.genes.2021.RData", package = "FIERCE", lib.loc= .libPaths()))

perform_preprocessing(adata, MT_prefix='mt-', cell_cycle_scoring=TRUE, s_genes=mouse.cc.genes.2021$s.genes, g2m_genes=mouse.cc.genes.2021$g2m.genes, perform_regression=FALSE, n_pcs=6, perform_clustering=FALSE, color_as=c("n_genes_by_counts","total_counts","phase","ClusterName"), lab_order=list(NULL,NULL,NULL,c("RadialGlia2","RadialGlia","nIPC","Nbl1","Nbl2","ImmGranule1","ImmGranule2","Granule")),palette=list(NULL,NULL,c("#6baed6","#fcbba1","red"),c('#e31a1c','orange','#ffeda0','#fcc5c0','#c6dbef','#6baed6','#2171b5','#08306b')))

adata <- compute_velocity(adata)

adata <- plot_velocity(adata, color_as=c("ClusterName","phase"))

adata <- compute_signaling_entropy(adata, species="Mouse", compute_potency_states=TRUE, phenotype_annotation="ClusterName")

plot_entropy_results(adata, phenotype_annotation="ClusterName")

adata <- compute_entropy_UMAP(adata, n_neighbors=c(10,30,100), n_pcs=c(15,20,30,50), perform_clustering=FALSE, color_as ="ClusterName")

adata <- compute_graph_and_stream(adata, n_neighbors_graph=100, n_pcs_graph=50, n_neighbors_emb=100, n_pcs_emb=50, color_as=c("ClusterName","phase"))

save_h5ad(adata,"adata.h5ad")


