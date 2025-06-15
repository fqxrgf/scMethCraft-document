Batch Integration
================================================

Environment Setup and Library Import
------------------------------------

Import necessary modules for single-cell methylation analysis.

.. code-block:: python

   import scanpy as sc
   from scMethCraft.model.scmethcraft_model import *
   from scMethCraft.model.scmethcraft_trainning import *
   import scMethCraft.model.methyimp as mp
   from scMethCraft.function.embedding import *
   from scMethCraft.function.batch import *

Data Loading
------------

Load the preprocessed single-cell methylation data.

.. code-block:: python

   input_path = f"../project/sample_data/genome/"
   raw_adata = sc.read(f"{input_path}/adata.h5ad")

Parameter Initialization
------------------------

Set key parameters for the analysis.

.. code-block:: python

   cell = raw_adata.shape[0]
   kmer_k = 8
   seq_length = 10000

Model Loading
-------------

Initialize and load pre-trained scMethCraft models.

.. code-block:: python

   scMethCraft_part1 = Sequence_extraction(cell, K=kmer_k, genomic_seq_length=seq_length).to(device)
   scMethCraft_part2 = Similarity_weighting(cell, dropout_rate=0.5).to(device)
   modelpath = f"../project/sample_data/output/"
   scMethCraft_part1.load_state_dict(torch.load(f"{modelpath}/scMethCraft_part1.pth"))
   scMethCraft_part2.load_state_dict(torch.load(f"{modelpath}/scMethCraft_part2.pth"))



Batch Effect Correction
-----------------------

.. code-block:: python

   adata = output_batch_integration(raw_adata, scMethCraft_part1, scMethCraft_part2)

