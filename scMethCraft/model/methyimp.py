import pandas as pd
import scanpy as sc
import episcanpy as epi
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn_extra.cluster import KMedoids
import sys
import os

def evaluation(true_label,pred_label,adata = None):
    result = dict()
    result["ARI"] =  metrics.adjusted_rand_score(true_label,pred_label)
    result["AMI"] =  metrics.adjusted_mutual_info_score(true_label,pred_label)
    result["NMI"] =  metrics.normalized_mutual_info_score(true_label,pred_label)
    result["homo"] =  metrics.homogeneity_score(true_label,pred_label)
    result["FMI"] =  metrics.fowlkes_mallows_score(true_label,pred_label)
    return result

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
        
def create_umap_pca(adata,color,show = True):
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata,color = color,show = show)

def create_umap_embedding(adata,color,embedding,show = True):
    sc.pp.neighbors(adata,use_rep = embedding)
    sc.tl.umap(adata)
    sc.pl.umap(adata,color = color,show = show)

def cluster(data,method,n = None):
    with HiddenPrints():
        if 'Dleiden' in method:
            epi.tl.leiden(data, key_added='Dleiden')
            label = data.obs['Dleiden']
        if 'Dlouvain' in method:
            epi.tl.louvain(data, key_added='Dlouvain')
            label = data.obs['Dlouvain']
        if 'leiden' in method:
            epi.tl.getNClusters(data, n, method='leiden')
            label = data.obs['leiden']
        if 'louvain' in method:
            epi.tl.getNClusters(data, n, method='louvain')
            label = data.obs['louvain']
        if 'kmedoids' in method:
            kmedoids = KMedoids(n_clusters=n)
            ydata = kmedoids.fit_predict(data.obsm["X_pca"])
            label = kmedoids.labels_
            adata.obs["kmedoids"] = label.values
        if 'kmeans' in method:
            clf = KMeans(n_clusters=n)
            ydata = clf.fit_predict(data.obsm["X_pca"])
            label = clf.labels_
            adata.obs["kmeans"] = label.values




