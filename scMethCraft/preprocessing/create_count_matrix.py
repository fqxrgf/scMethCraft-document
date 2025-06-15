import pandas as pd
import numpy as np
from scipy import spatial
from tqdm import tqdm
import scanpy as sc
import time
import gzip
from functools import reduce
import os
import h5py
import random
import pysam
import anndata as ad
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

MOUSE = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',  
        '11', '12', '13', '14', '15', '16', '17', '18', '19','X', 'Y', 'M']

mm9_size = [197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172, 129993255, 121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430, 166650296, 15902555, 16299]

mm10_size = [195471971, 182113224, 160039680, 156508116,
              151834684, 149736546, 145441459, 129401213, 124595110,
              130694993, 122082543, 120129022, 120421639, 124902244,
              104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698, 16299]

HUMAN = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',  
        '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22','X', 'Y', 'M']

hg38_size = [248956422, 242193529, 198295559, 190214555,
              181538259, 170805979, 159345973, 145138636, 138394717,
              133797422, 135086622, 133275309, 114364328, 107043718,
              101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
              46709983, 50818468, 156040895, 57227415, 16569]

hg19_size = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,  146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983,63025520, 48129895, 51304566, 155270560,59373566, 16571]

def name_features(loaded_features):
    """
    From EpiScanpy
    Extract the names of the loaded features, specifying the chromosome they originated from.
    It also contain the feature coordinates and an unique identifier.
    """
    feat_names = []
    i = 0
    for c in loaded_features.keys():
        for name in loaded_features[c]:
            #add_name = '_'.join(['chr', c, name[-1].rstrip('\n'), str(i)])
            add_name = '_'.join([c, str(name[0]), str(name[1])])
            add_name ='chr' + add_name
            # the feature names will be like chr1_1234_1245
            if add_name[-1] =='\n':
                add_name = add_name[:-1]
            feat_names.append(add_name)
            i += 1
    return(feat_names)
    
def process_cell(args):
    
    cell = args[0]
    path = args[1]
    chromosome = args[2]
    meth_context = args[3]
    annotation = args[4]
    nb_annotation = args[5]
    threshold = args[6]
    output_file = args[7]
    data_header = args[8]
    
    try:
        if meth_context == 'CG':
            tmp_file = read_meth_fileCG(cell, path, chromosome,**data_header)

        elif meth_context == 'CH':
            tmp_file = read_meth_fileCH(cell, path, chromosome)
        else:
            return None

    except Exception as exc:
        print(cell, 'failed with exception:', exc)
        return None

    result = []
    for index_annot in range(nb_annotation):
        meth_level_annot = methylation_level(tmp_file, annotation[index_annot], chromosome, threshold[index_annot])
        if type(output_file) == list:
            write_methlevel(meth_level_annot, output_file[index_annot], cell, writing_option[index_annot], feature_names[index_annot])
        else:
            result.append(np.matrix(meth_level_annot))
    print(cell,'processed')
    return result

def read_meth_fileCG(sample_name, path, chromosome, chrom=0, pos=1, met=4, tot=5, status=3):
    """
    Read file from which you want to extract the methylation level and
    (assuming it is like the Ecker/Methylpy format) extract the number of
    methylated read and the total number of read for the cytosines covered and
    in the right genomic context (CG or CH)
    Parameters
    ----------
    sample_name:
        name of the file to read to extract key information.
    meth_type:
        CG, CH or not specified
    head: if there is header that you don't want to read. An annotation in the
        file you are reading. The default value is the Methylpy/Ecker header
    path: path of the to access the file of the sample you want to read. 
    chromosome: chromosomes if the species you are considering. default value
        is the human genome (including mitochondrial and sexual chromosomes)
    """
    reduced_cyt = {key: [] for key in chromosome} # to store cyt for each chrom (intermed output)

    with gzip.open(path+sample_name) as sample:
        for line in sample:
            line = str(line, encoding="utf-8").split('\t')
            if len(line) < 3:
                print("The separator of input may not be \t or coutains wrong columns")
            break
        
        
    with gzip.open(path+sample_name) as sample:
        j=0
        for line in sample:
            try:
                line = str(line, encoding="utf-8").split('\t')
                try:
                    if (line[status] in ['CGG', 'CGC', 'CGA', 'CGT']):
                        if 'chr' not in line[chrom]:
                            line[chrom] = 'chr'+line[chrom]
                        if(line[chrom][3:] in chromosome):
                            reduced_cyt[line[chrom][3:]].append((int(line[pos]), int(line[met]), int(line[tot])))
#                             print((int(line[pos]), int(line[met]), int(line[tot])))
                except:
                    j+=1
            except:
                j+=1
    return(reduced_cyt)

def methylation_level(reduced_cyt, feature, chromosome, threshold=1):
    """
    Measure the methylation for the feature you give as input using the reduce
    representation of the sample cytosine (output of read_methylation_file)
    Parameters
    ----------
    reduced_cyt: datatype that contained processed sample file. It only contains
        the cytosines that were in the genomic context you wanted to filter for.
        (output of read_methylation_file function).
    feature: the feature in the right datatype for which you want to determine
        the methylation level.
    chromosome: chromosomes if the species you are considering. default value
        is the human genome (including mitochondrial and sexual chromosomes).
    """
    ## actually, to write sparse matrix I need a list, not a dictionary
    #meth_levels_bins = {key:[] for key in chromosome}
    meth_levels_bins = []
    for c in chromosome:
        meth_reads = np.zeros(len(feature[c]))
        tot_reads = np.zeros(len(feature[c]))
        nb_cyt = np.zeros(len(feature[c]))
        cytosines = reduced_cyt[c] 
        i = 0
        for j in range(len(feature[c])): # for every bins in a given chrom
            meth_reads = 0.0
            tot_reads = 0.0
            nb_cyt = 0
            # I am skipping cytosine that are before the beginning 
            # of my current bin. 
            while (i < len(cytosines)) and cytosines[i][0] < feature[c][j][0]:
                i += 1
            # Once I got one cytosine that is after the beginning of my feature 
            # I need to check if this feature is within the enhancer limits
            # so if the position of the cytosine is not > to the end of the feature
            if i<len(cytosines) and cytosines[i][0] <= feature[c][j][1]:
                meth_reads += cytosines[i][-2] # meth cyt read
                tot_reads += cytosines[i][-1] # tot cyt read
                nb_cyt += 1 # nb of cyt

            # check if the next cytosine fall into the current feature is important
            # to this end, I have another pointer/iterator k. 
            # at the next feature I will have to check from i but for the current 
            # feature I need to check the next cytosine and I use the variable k for 
            # this.
            k = i+1
            while k < len(cytosines) and cytosines[k][0] <= feature[c][j][1]:
                meth_reads += cytosines[k][-2] # meth cyt read
                tot_reads += cytosines[k][-1] # tot cyt read
                nb_cyt += 1  # nb of cyt
                k += 1
            ## actually, to write sparse matrix I need a list, not a dictionary
            if nb_cyt >= threshold:
                #meth_levels_bins[c].append(format(meth_reads/tot_reads, '.3f'))
                meth_levels_bins.append(format(meth_reads/tot_reads, '.3f'))
            else:
                #meth_levels_bins[c].append(np.nan)
                meth_levels_bins.append(np.nan)
    return(meth_levels_bins)

def make_windows(size,chromosomes,chromosome_sizes):
    """
    Adapted from EpiScanpy
    Generate windows/bins of the given size for the appropriate genome (default
    choice is human). 
    
    """
    features_chrom = {}
        
    for c in range(len(chromosomes)):
        start = range(1, chromosome_sizes[c] - size, size)
        end = range(size, chromosome_sizes[c], size)
        features_chrom[chromosomes[c]] = [[start[i], end[i], ''.join(["chr", chromosomes[c], "_", str(start[i]), "_", str(end[i])])] for i in range(len(end))]
        
    return(features_chrom)

def build_count_mtx(cells, annotation, path="", output_file=None, writing_option="a",
                    meth_context="CG", chromosome=None, feature_names=None,
                   threshold=1, ct_mtx=None, sparse=False,data_header=None):

    #verbosity
    i = 0
    
    #################################
    if type(annotation) != list:
        annotation = [annotation]
        output_file = [output_file]
        ct_mtx = [ct_mtx]
        feature_names = [feature_names]
    nb_annotation = len(annotation)
    
    if type(writing_option) != list:
        writing_option = [writing_option for x in range(nb_annotation)]
    if type(threshold) != list:
        threshold = [threshold for x in range(nb_annotation)]
    if (output_file != None):
        if (type(output_file) != list):
            output_file = [output_file]
        
    file_tuples = []
    for cell in cells:
        file_tuples.append((cell,path, chromosome,meth_context,annotation,nb_annotation,threshold,output_file,data_header))

    #################################
    
    with ProcessPoolExecutor(max_workers = 4) as executor:
        results = executor.map(process_cell, file_tuples)
    
    i = 0
    ct_mtx = None
    for cell, result in zip(cells, results):
        if result is not None:
            for index_annot, meth_level_annot in enumerate(result):
                if type(output_file) == list:
                    write_methlevel(meth_level_annot, output_file[index_annot], cell, writing_option[index_annot], feature_names[index_annot])
                else:
                    if ct_mtx is None:
                        ct_mtx = [np.matrix(meth_level_annot)]
                    elif index_annot >= len(ct_mtx):
                        ct_mtx.append(np.matrix(meth_level_annot))
                    else:
                        ct_mtx[index_annot] = np.vstack([ct_mtx[index_annot], meth_level_annot])

            i += 1
                
    if ct_mtx != None:
        return(ct_mtx)
    else:
        return()
    
def replace_bedgz(name):
    replace_name = name.split('.', 1)[0]
    return replace_name

def count_matrix(bed_path,bed_files = "total",windows_width=10000,refgenome='hg38',data_header={"chrom":0, "pos":1, "met":4, "tot":5, "status":3}):
    
    #定义参考基因组size
    if refgenome =='mm9':
        refgenome_size = mm9_size
        organism = 'MOUSE'
    elif refgenome =='mm10' :
        refgenome_size = mm10_size
        organism = 'MOUSE'
    elif refgenome =='hg19':
        refgenome_size = hg19_size
        organism = 'HUMAN'
    elif refgenome =='hg38':
        refgenome_size = hg38_size
        organism = 'HUMAN'
                
    else:
        raise Exception("illegal reference genome")

    if bed_files == "total":
        bed_dir = os.listdir(bed_path)
    else:
        bed_dir == bed_files

    if '.ipynb_checkpoints' in bed_dir:
        bed_dir.remove('.ipynb_checkpoints')

    windows = make_windows(windows_width,chromosomes = eval(organism),chromosome_sizes=refgenome_size)
    w_names =name_features(windows) 


    w_mtx = build_count_mtx(cells=bed_dir,
                               annotation=[windows],
                               path=bed_path,
                               chromosome=HUMAN,
                               output_file=None,
                               meth_context='CG',
                                data_header = data_header,
                               threshold=[1])# minimum number of cytosine/reads to have at any given feature to 
                                # not consider the feature to have a missing methylation level
    
    cell_name = list(map(replace_bedgz,bed_dir))
    
    w_mtx=np.array(w_mtx[0]).astype('float64')
    
    legal_region = ((~np.isnan(w_mtx)).sum(axis = 0) != 0)
    
    w_mtx= w_mtx[:,legal_region]
    
    try:
        adata_w = ad.AnnData(w_mtx,  var=pd.DataFrame(index=w_names)[legal_region],obs=pd.DataFrame(index=cell_name))
    except:
        adata_w = ad.AnnData(w_mtx,  var=pd.DataFrame(index=w_names)[legal_region])
    
    adata_w.var[['chr','start','end']] = pd.Series(adata_w.var.index).str.split("_",expand=True).values
    
    return adata_w

def assemble_meta(data,meta):
    assemble_meta = meta.loc[data.obs.index,:]
    if assemble_meta.shape[0] == data.shape[0]:
        data.obs = assemble_meta   
        print("Successfully matched meta information")
    else:
        raise Exception("Meta information not match")
        
def filter_region(adata,threshold_na = 0.8,threshold_var = 5,threshold_sim = 0.99):
    adata = adata[:,(np.isnan(adata.X).sum(axis = 0)< adata.shape[0]*threshold_na)].copy()
    variance = np.nanvar(adata.X,axis = 0)
    adata = adata[:,variance>np.percentile(variance,threshold_var)].copy()

    adata_temp = adata.copy().T
    medians  = np.nanmedian(adata_temp.X,axis = 1)
    adata_temp.X = np.nan_to_num(adata_temp.X-medians.reshape(-1,1))+medians.reshape(-1,1)
    sc.pp.pca(adata_temp)

    feature_hash_list = adata_temp.obsm["X_pca"].tolist()

    i = 1
    data_index = list(range(len(feature_hash_list)))

    while True:    
        if i >= len(feature_hash_list):
            break
        for j in range(max(0,i-100),min(len(feature_hash_list),i+100)):
            if i != j:
                if  1-spatial.distance.cosine(feature_hash_list[i],feature_hash_list[j]) > threshold_sim:
                    feature_hash_list.pop(i)
                    data_index.pop(i)
                    i = i-1
                    break
        i = i+1

    return(adata[:,data_index].copy())

