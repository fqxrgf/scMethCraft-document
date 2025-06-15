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

def obtain_seq(input_fasta,input_ad,output_path):
    h5_name = '%s/all_seqs.h5'%output_path
    seq_len=10000
    batch_size = 1000
    K = 8
    os.makedirs(output_path, exist_ok=True)
    
    ad = input_ad
    ad.write('%s/adata.h5ad'%output_path)
    ad_all = ad.copy()
    ad.var.loc[:,['chr','start','end']].to_csv('%s/peaks.bed'%output_path, sep='\t', header=False, index=False)
    print('successful writing bed file.')
    
    n_peaks = ad.shape[1]
    bed_df = ad.var.loc[:,['chr','start','end']] # bed file
    bed_df.index = np.arange(bed_df.shape[0])
    n_batch = int(np.floor(n_peaks/batch_size))
    batches = np.array_split(np.arange(n_peaks), n_batch)
    
    m_all = ad.X.T
    np.save('%s/m_all'%output_path, m_all)
    print('successful writing sparse m.')
    make_h5_sparse(ad_all,h5_name, input_fasta,seq_len=seq_len,batch_size=1000)
    
    f = h5py.File(h5_name, "r+")
    pos_ = ad.var[["chr","start","end"]].copy()
    pos_[pos_.columns[0]] = pd.Series([x[3:] for x in pos_["chr"]]).values
    pos_[pos_.columns[0]] = pos_[pos_.columns[0]].map(lambda x: {"X":23,"Y":24,"M":25}.get(x, x))
    f["Pos"] = pos_.astype(int)
    
    f.create_dataset(
        "Kmer",
        (n_peaks, 4**K),
        dtype="int8",
    )

    for i in tqdm(range(len(batches))):
        idx = batches[i]
        seqs_dna,_ = make_bed_seqs_from_df(
            bed_df.iloc[idx,:],
            fasta_file=input_fasta,
            seq_len=seq_len,
        )
        dna_array_dense = {x:summary_kmers(build_kmers(seqs_dna[x-min(idx)],K)) for x in idx}
        dna_array_dense = pd.DataFrame(dna_array_dense).T
        dna_array_dense =pd.DataFrame(np.nan_to_num(dna_array_dense),columns = dna_array_dense.columns,index = dna_array_dense.index)
        result = dna_array_dense


        a = reduce(lambda x,y: [i+j for i in x for j in y], [['A','T','C','G']] * K)
        result = result.reindex(columns=a, fill_value=0)
        result[~pd.notna(result)] = 0
        result = result.astype(int)
        f["Kmer"][idx] = result

    f.close()
    
def make_bed_seqs_from_df(input_bed, fasta_file, seq_len, stranded=False): 
    """Return BED regions as sequences and regions as a list of coordinate
    tuples, extended to a specified length."""
    """Extract and extend BED sequences to seq_len."""
    fasta_open = pysam.Fastafile(fasta_file)

    seqs_dna = []
    seqs_coords = []

    for i in range(input_bed.shape[0]):
        chrm = input_bed.iloc[i,0]
        start = int(input_bed.iloc[i,1])
        end = int(input_bed.iloc[i,2])
        strand = "+"

        # determine sequence limits
        mid = (start + end) // 2
        seq_start = mid - seq_len // 2
        seq_end = seq_start + seq_len

        # save
        if stranded:
            seqs_coords.append((chrm, seq_start, seq_end, strand))
        else:
            seqs_coords.append((chrm, seq_start, seq_end))
        # initialize sequence
        seq_dna = ""
        # add N's for left over reach
        if seq_start < 0:
            print(
                "Adding %d Ns to %s:%d-%s" % (-seq_start, chrm, start, end),
                file=sys.stderr,
            )
            seq_dna = "N" * (-seq_start)
            seq_start = 0

        # get dna
        seq_dna += fasta_open.fetch(chrm, seq_start, seq_end).upper()

        # add N's for right over reach
        if len(seq_dna) < seq_len:
            print(
                "Adding %d Ns to %s:%d-%s" % (seq_len - len(seq_dna), chrm, start, end),
                file=sys.stderr,
            )
            seq_dna += "N" * (seq_len - len(seq_dna))
        # append
        seqs_dna.append(seq_dna)
    fasta_open.close()
    return seqs_dna, seqs_coords

def dna_1hot_2vec(seq, seq_len=None):
    """dna_1hot
    Args:
      seq:       nucleotide sequence.
      seq_len:   length to extend/trim sequences to.
      n_uniform: represent N's as 0.25, forcing float16,
                 rather than sampling.
    Returns:
      seq_code: length by nucleotides array representation.
    """
    if seq_len is None:
        seq_len = len(seq)
        seq_start = 0
    else:
        if seq_len <= len(seq):
            # trim the sequence
            seq_trim = (len(seq) - seq_len) // 2
            seq = seq[seq_trim : seq_trim + seq_len]
            seq_start = 0
        else:
            seq_start = (seq_len - len(seq)) // 2
    seq = seq.upper()

    # map nt's to a matrix len(seq)x4 of 0's and 1's.
    seq_code = np.zeros((seq_len, ), dtype="int8")

    for i in range(seq_len):
        if i >= seq_start and i - seq_start < len(seq):
            nt = seq[i - seq_start]
            if nt == "A":
                seq_code[i] = 0
            elif nt == "C":
                seq_code[i] = 1
            elif nt == "G":
                seq_code[i] = 2
            elif nt == "T":
                seq_code[i] = 3
            else:
                seq_code[i] =  random.randint(0, 3)
    return seq_code

def split_train_test_val(ids, seed=10, train_ratio=0.9):
    np.random.seed(seed)
    test_val_ids = np.random.choice(
        ids,
        int(len(ids) * (1 - train_ratio)),
        replace=False,
    )
    train_ids = np.setdiff1d(ids, test_val_ids)
    val_ids = np.random.choice(
        test_val_ids,
        int(len(test_val_ids) / 2),
        replace=False,
    )
    test_ids = np.setdiff1d(test_val_ids, val_ids)
    return train_ids, test_ids, val_ids

def make_h5_sparse(tmp_ad, h5_name, input_fasta, seq_len=1344, batch_size=1000):
    ## batch_size: how many peaks to process at a time
    ## tmp_ad.var must have columns chr, start, end
    
    t0 = time.time()
    
    m = tmp_ad.X
    #m = m.tocoo().transpose().tocsr()
    n_peaks = tmp_ad.shape[1]
    bed_df = tmp_ad.var.loc[:,['chr','start','end']] # bed file
    bed_df.index = np.arange(bed_df.shape[0])
    n_batch = int(np.floor(n_peaks/batch_size))
    batches = np.array_split(np.arange(n_peaks), n_batch) # split all peaks to process in batches
    
    ### create h5 file
    # X is a matrix of n_peaks * 1344
    f = h5py.File(h5_name, "w")
    
    ds_X = f.create_dataset(
        "X",
        (n_peaks, seq_len),
        dtype="int8",
    )

    # save to h5 file
    for i in range(len(batches)):
        
        idx = batches[i]
        # write X to h5 file
        seqs_dna,_ = make_bed_seqs_from_df(
            bed_df.iloc[idx,:],
            fasta_file=input_fasta,
            seq_len=seq_len,
        )
        dna_array_dense = [dna_1hot_2vec(x) for x in seqs_dna]
        dna_array_dense = np.array(dna_array_dense)
        ds_X[idx] = dna_array_dense
            
        t1 = time.time()
        total = t1-t0
        print('process %d peaks takes %.1f s' %((i+1)*batch_size, total))
    
    f.close()

def summary_kmers(kmers):
    """a function to summarize the kmers"""
    kmers_stat = dict(Counter(kmers))
    return kmers_stat

def build_kmers(seq, k_size):
    """a function to calculate kmers from seq"""
    kmers = []  # k-mer存储在列表中
    n_kmers = len(seq) - k_size + 1
    
    for i in range(n_kmers):
        kmer = seq[i:i + k_size]
        kmers.append(kmer)
        
    return kmers