Training Data Generation
========================

.. contents::
   :local:
   :depth: 2

Environment Setup and Library Import
------------------------------------

First, import the necessary modules and set up the Python environment.

.. code-block:: python

    from scMethCraft.preprocessing.create_count_matrix import *
    from scMethCraft.preprocessing.retrive_sequence import *

Import two core modules from scMethCraft:  
- ``create_count_matrix``: For creating count matrices  
- ``retrive_sequence``: For retrieving genomic sequences

Loading Sample Metadata
------------------------------------

Read the sample metadata file.

.. code-block:: python

    meta = pd.read_csv("../sample_data/meta.csv", index_col=0)

Use pandas to read the metadata file in CSV format.  
- ``index_col=0`` indicates using the first column as row index  
- File path is relative: ``../sample_data/meta.csv``

Creating Methylation Count Matrix
------------------------------------

Generate a methylation count matrix from BED files.

.. code-block:: python

    count_matrix = count_matrix("../sample_data/bed/")

Call the ``count_matrix`` function  
- Parameter: the directory path containing BED files  
- Results are stored in the ``count_matrix`` variable

Integrating Metadata with Count Matrix
--------------------------------------------

Associate metadata information with the count matrix.

.. code-block:: python

    assemble_meta(count_matrix, meta)

Call the ``assemble_meta`` function  
- Parameter 1: the count matrix generated in the previous step  
- Parameter 2: the loaded metadata DataFrame  
- This function will integrate metadata into the count matrix

Filtering Regions
------------------------------------

.. code-block:: python

    count_matrix = filter_region(count_matrix)

Use ``filter_region`` to remove low-quality or invalid regions from the count matrix.

Retrieving Genomic Sequences
------------------------------------

Extract sequences of specific regions from the reference genome.

.. code-block:: python

    obtain_seq("/home/sccasimp/data/songzhuofan/hg38.fa", 
               filter_reigon(count_matrix),
               "../sample_data/genome/")

Call the ``obtain_seq`` function with three parameters:  

- Reference genome file path: ``hg38.fa``  
- Filtered genomic regions: obtained by processing ``count_matrix`` with ``filter_reigon``  
- Output directory path  

This step extracts sequences corresponding to methylation sites from the reference genome.
