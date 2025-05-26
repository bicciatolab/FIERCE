from tqdm import tqdm
import ipywidgets
import scvelo as scv
import scanpy as sc
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
import pandas as pd
import numpy as np
import anndata as ad
from scipy.sparse import issparse


def subset_anndata_genes(adata, geneset, adata_copy=True):
	if adata_copy==True:
		adata = adata.copy()
		adata = adata[:, geneset]
		adata.obs["value"] = 0
	else:
		adata = adata[:, geneset]
		adata.obs["value"] = 0
	return adata


def subset_anndata_cells(adata, cellset, adata_copy=True):
	if adata_copy==True:
		adata = adata.copy()
		adata = adata[cellset, :]
		adata.obs["value"] = 0
	else:
		adata = adata[cellset, :]
		adata.obs["value"] = 0
	return adata


def add_layer(adata, new_layer, new_layer_name, transpose=False, compress=False):
	if transpose == True:
		new_layer = new_layer.T
	if compress == True:
		from scipy import sparse
		new_layer = sparse.csc_matrix(new_layer)
	adata.layers[new_layer_name] = new_layer
	return None


def add_annotation_obs(adata, new_annot, new_annot_name):
	adata.obs[new_annot_name] = new_annot
	return None


def add_obsm(adata, new_obsm, new_obsm_name):
	adata.obsm[new_obsm_name] = new_obsm
	return None


def add_X(adata, new_data, transpose=False, compress=False):
	if transpose == True:
		new_data = new_data.T
	if compress == True:
		from scipy import sparse
		new_data = sparse.csc_matrix(new_data)
	adata.X = new_data
	return None


def add_uns(adata, new_uns, new_uns_name, nested=False, sub_slot=None):
	if nested == False:
		adata.uns[new_uns_name] = new_uns
	else:
		if sub_slot is None:
			print("Error: name of uns sub-slot not specified")
		else:
			adata.uns[new_uns_name][sub_slot] = new_uns
	return None


def delete_annotation_obs(adata, annotation_to_delete):
	del(adata.obs[annotation_to_delete])
	return None


def delete_annotation_var(adata, annotation_to_delete):
	del(adata.var[annotation_to_delete])
	return None


def delete_layer(adata, layer_to_delete):
	del(adata.layers[layer_to_delete])
	return None


def delete_uns(adata, uns_to_delete):
	del(adata.uns[uns_to_delete])
	return None


def delete_obsm(adata, obsm_to_delete):
	del(adata.obsm[obsm_to_delete])
	return None


def add_annotation_var(adata, new_annot, new_annot_name):
	adata.var[new_annot_name] = new_annot
	return None


def add_varm(adata, new_varm, new_varm_name):
	adata.varm[new_varm_name] = new_varm
	return None


def add_obsp(adata, new_obsp, new_obsp_name):
	adata.obsp[new_obsp_name] = new_obsp
	return None


def replace_cell_names(adata, new_cell_names):
	adata.obs_names = new_cell_names


def check_layers(adata, query):
	if query in adata.layers.keys():
		return "Yes"
	else:
		return "No"


def check_uns(adata, query):
	if query in adata.uns.keys():
		return "Yes"
	else:
		return "No"


def extract_obsm_keys(adata):
	keys_list = []
	for i in adata.obsm.keys():
		keys_list.append(i)
	return keys_list


def add_shelves(my_shelf, new_shelf_name, new_shelf):
	my_shelf[new_shelf_name] = new_shelf
	return None


def subset_shelf_genes(my_shelf, shelf_name, gene_list):
	my_shelf[shelf_name] = my_shelf[shelf_name][:,gene_list]
	return None


def key_in_shelf(my_shelf, key):
	if key in my_shelf:
		return "yes"
	else:
		return "no"


def log1p_correction(adata):
	adata.uns['log1p'] = {'base': None}
	return None


def plot_multiple_scatters(adata, embedding_basis, genes, layer, legend_loc, filename, color_map):
	pp = PdfPages(filename)
	for gene in genes:
		sc.pl.scatter(adata, basis=embedding_basis, color=gene, layers=layer, color_map=color_map, frameon=False, legend_loc=legend_loc)
		plt.savefig(pp, format='pdf', bbox_inches='tight')
	pp.close()


def pca_heatmap(adata, components, use_raw=None, layer=None, filename="pca_heatmap.pdf"):
	adata.obs['fake'] = np.repeat('fake', adata.shape[0]).tolist()
	attr = 'varm'
	keys = 'PCs'
	components = list(range(components)) 
	pp = PdfPages(filename)
	for component in components:
		scores = getattr(adata, attr)[keys][:, component]
		dd = pd.DataFrame(scores, index=adata.var_names)
		var_names_pos = dd.sort_values(0, ascending=False).index[:20]
		var_names_neg = dd.sort_values(0, ascending=True).index[:20]
		pd2 = pd.DataFrame(adata.obsm['X_pca'][:, component], index=adata.obs.index)
		bottom_cells = pd2.sort_values(0).index[:300].tolist()
		top_cells = pd2.sort_values(0, ascending=False).index[:300].tolist()
#		ax = plt.subplot()
		plt.rcParams.update({'font.size' : 8})
		sc.pl.heatmap(
			adata[top_cells+bottom_cells], list(var_names_pos) + list(var_names_neg),
			show_gene_labels=True,
			swap_axes=True, cmap='viridis',
			use_raw=False, layer=layer, figsize=(3, 3), groupby='fake')
#		for label in (ax.get_xticklabels() + ax.get_yticklabels()):
#    		label.set_fontproperties(font_prop)
#			label.set_fontsize(8) # Size here overrides font_prop
		plt.title('PC_' + str(component+1))
		plt.savefig(pp, format='pdf', bbox_inches='tight')
	pp.close()
	del(adata.obs['fake'])


def draw_embedding(adata, file_name, emb='umap', stream=False, c_as=None, c_type=None, pal=None, c_ord=None, vkey='velocity_of_the_entropy', legend_loc="right margin", alpha=0.3, min_mass=4, add_outline=False, outline_width=(0.1, 0.01)):
#	from matplotlib.backends.backend_pdf import PdfPages
	if ((c_as is not None) and (c_type == 'categorical')):
		if str(adata.obs[c_as].dtype) != 'category':
			adata.obs[c_as] = pd.Series(adata.obs[c_as], dtype ="category")
		if c_ord is not None:
			if pal is not None:
				adata.uns[str(c_as + '_colors')] = np.array(pal, dtype=object)
			elif str(c_as + '_colors') in adata.uns.keys():
#				inds = [j for j, val in enumerate(adata.obs[c_as].cat.categories == c_ord[0]) if val][0]
#				new_colors = adata.uns[str(c_as + '_colors')][inds]
#				new_colors = new_colors[adata.obs[c_as].cat.categories == c_ord[0]]
				new_colors = []
				for i in list(range(len(c_ord))):
					inds = [j for j, val in enumerate(adata.obs[c_as].cat.categories == c_ord[i]) if val][0]
					new_c = adata.uns[str(c_as + '_colors')][inds]
					new_colors.append(new_c)
#					new_colors.append(adata.uns[str(c_as + '_colors')][adata.obs[c_as].cat.categories == c_ord[i]][0])
				adata.uns[str(c_as + '_colors')] = np.array(new_colors, dtype=object)
			else:
				my_palette = ["#2171b5","#6baed6","#c6dbef","#cb181d","yellow","#fcbba1","orange","#74c476","#c7e9c0","#525252","#969696","#d9d9d9", "#762a83","#9970ab","#c2a5cf","#e7d4e8","#d9f0d3","#a6dba0","#5aae61","#1b7837","#f7f4f9","#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f","#7f2704"]
				adata.uns[str(c_as + '_colors')] = np.array(my_palette[0:len(adata.obs[c_as].cat.categories)], dtype=object)
			adata.obs[c_as] = adata.obs[c_as].cat.set_categories(c_ord, ordered=True)
		elif pal is not None:
			adata.uns[str(c_as + '_colors')] = np.array(pal, dtype=object)
		else:
			if str(c_as + '_colors') not in adata.uns.keys():
				my_palette = ["#2171b5","#6baed6","#c6dbef","#cb181d","yellow","#fcbba1","orange","#74c476","#c7e9c0","#525252","#969696","#d9d9d9", "#762a83","#9970ab","#c2a5cf","#e7d4e8","#d9f0d3","#a6dba0","#5aae61","#1b7837","#f7f4f9","#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f","#7f2704"]
				adata.uns[str(c_as + '_colors')] = np.array(my_palette[0:len(adata.obs[c_as].cat.categories)], dtype=object)
	elif ((c_as is not None) and (c_type == 'numeric')):
		colmap = 'viridis'
	else:
		colmap = None
	if stream==False:
		plt.rcdefaults()
		plt.rcParams.update({'legend.markerscale':1})
		if emb=='umap':
			sc.pl.umap(adata, color=c_as, ncols=1, frameon=False, legend_loc=legend_loc, legend_fontsize=10, size=2*round(120000/len(adata.obs_names.tolist())), save=str(file_name+'.pdf'), add_outline=add_outline, outline_width=outline_width, alpha=alpha, color_map=colmap)
		elif ('umap_entropy' in emb):
			adata_temp = ad.AnnData(adata.layers['partial_entropies_observed'])
			adata_temp.obs_names = adata.obs_names
			adata_temp.var_names = adata.var_names
			adata_temp.obsm['X_umap'] = adata.obsm[str('X_'+emb)]
			adata_temp.uns['umap'] = adata.uns[emb]
			adata_temp.obs[c_as] = adata.obs[c_as]
			if (c_type == 'categorical'):
				adata_temp.uns[str(c_as + '_colors')] = adata.uns[str(c_as + '_colors')]
			sc.pl.umap(adata_temp, color=c_as, ncols=1, frameon=False, legend_loc=legend_loc, legend_fontsize=10, size=2*round(120000/len(adata.obs_names.tolist())), save=str(file_name+'.pdf'), add_outline=add_outline, outline_width=outline_width, alpha=alpha, color_map=colmap)
			del adata_temp
		elif emb=='pca':
			sc.pl.pca(adata, color=c_as, ncols=1, frameon=False, legend_loc=legend_loc, legend_fontsize=10, size=2*round(120000/len(adata.obs_names.tolist())), save=str(file_name+'.pdf'), add_outline=add_outline, outline_width=outline_width, alpha=alpha, color_map=colmap)
		else:
			adata_temp = ad.AnnData(adata.layers['partial_entropies_observed'])
			adata_temp.obs_names = adata.obs_names
			adata_temp.var_names = adata.var_names
			adata_temp.obsm['X_pca'] = adata.obsm['X_pca_entropy']
			adata_temp.uns['pca'] = adata.uns['pca_entropy']
			adata_temp.varm['PCs'] = adata.varm['PCs_entropy']
			adata_temp.obs[c_as] = adata.obs[c_as]
			if (c_type == 'categorical'):
				adata_temp.uns[str(c_as + '_colors')] = adata.uns[str(c_as + '_colors')]
			sc.pl.pca(adata_temp, color=c_as, ncols=1, frameon=False, legend_loc=legend_loc, legend_fontsize=10, size=2*round(120000/len(adata.obs_names.tolist())), save=str(file_name+'.pdf'), add_outline=add_outline, outline_width=outline_width, alpha=alpha, color_map=colmap)
			del adata_temp
	else:
		if legend_loc=="right margin":
			legend_loc = 'right'
		if legend_loc!="on data":
			plt.rcParams.update({'legend.loc':legend_loc})
		scv.pl.velocity_embedding_stream(adata, vkey=vkey, basis=emb, color=c_as, ncols=1, legend_loc=legend_loc, save=str(file_name+'.png'), alpha=alpha, min_mass=min_mass, size=2*round(120000/len(adata.obs_names.tolist())), add_outline=add_outline, outline_width=outline_width, dpi=500, color_map=colmap)


def set_annot_levels(adata, c_as, c_ord, pal=None):
	if c_ord is not None:
		if str(adata.obs[c_as].dtype) != 'category':
			adata.obs[c_as] = pd.Series(adata.obs[c_as], dtype ="category")
		if pal is not None:
			adata.uns[str(c_as + '_colors')] = np.array(pal, dtype=object)
		elif str(c_as + '_colors') in adata.uns.keys():
	#		inds = [j for j, val in enumerate(adata.obs[c_as].cat.categories == c_ord[0]) if val][0]
	#		new_colors = adata.uns[str(c_as + '_colors')][inds]
	#		new_colors = new_colors[adata.obs[c_as].cat.categories == c_ord[0]]
			new_colors = []
			for i in list(range(len(c_ord))):
				inds = [j for j, val in enumerate(adata.obs[c_as].cat.categories == c_ord[i]) if val][0]
				new_c = adata.uns[str(c_as + '_colors')][inds]
				new_colors.append(new_c)
	#			new_colors.append(adata.uns[str(c_as + '_colors')][adata.obs[c_as].cat.categories == c_ord[i]][0])
			adata.uns[str(c_as + '_colors')] = np.array(new_colors, dtype=object)
		else:
			my_palette = ["#2171b5","#6baed6","#c6dbef","#cb181d","yellow","#fcbba1","orange","#74c476","#c7e9c0","#525252","#969696","#d9d9d9", "#762a83","#9970ab","#c2a5cf","#e7d4e8","#d9f0d3","#a6dba0","#5aae61","#1b7837","#f7f4f9","#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f","#7f2704"]
			adata.uns[str(c_as + '_colors')] = np.array(my_palette[0:len(adata.obs[c_as].cat.categories)], dtype=object)
		adata.obs[c_as] = adata.obs[c_as].cat.set_categories(c_ord, ordered=True)
	elif pal is not None:
		adata.uns[str(c_as + '_colors')] = np.array(pal, dtype=object)
	else:
		if str(c_as + '_colors') not in adata.uns.keys():
			my_palette = ["#2171b5","#6baed6","#c6dbef","#cb181d","yellow","#fcbba1","orange","#74c476","#c7e9c0","#525252","#969696","#d9d9d9", "#762a83","#9970ab","#c2a5cf","#e7d4e8","#d9f0d3","#a6dba0","#5aae61","#1b7837","#f7f4f9","#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f","#7f2704"]
			adata.uns[str(c_as + '_colors')] = np.array(my_palette[0:len(adata.obs[c_as].cat.categories)], dtype=object)


