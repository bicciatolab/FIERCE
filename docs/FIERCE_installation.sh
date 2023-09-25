# FIERCE environment installation

 conda create --name FIERCE -c conda-forge r-base=4.0.5 python=3.10.6 scanpy=1.9.1 python-louvain=0.15 tqdm=4.64.1 pandas=1.5.1 scipy=1.9.3 numba=0.56.3 matplotlib-base=3.6.2 h5py=3.7.0 click=8.1.3 r-ggplot2=3.3.6 r-reshape2=1.4.4 r-scales=1.2.1 anndata=0.8.0 r-rgeos=0.5_9 r-igraph=1.3.4 python-igraph=0.10.2 r-leiden=0.4.3 r-rcurl=1.98_1.8 r-devtools=2.4.4 r-reticulate=1.15
#eval "$(conda shell.bash hook)"

conda activate FIERCE

python -m pip install cython
#pip3 install Cython==0.29.32 
pip3 install velocyto==0.17.17 scvelo scikit-learn

#wget https://github.com/bicciatolab/FIERCE/blob/main/docs/FIERCE_installation.R

#sudo apt install cmake    #needed to install SCENT and its dependencies

chmod +x ./FIERCE_installation.R

Rscript ./FIERCE_installation.R