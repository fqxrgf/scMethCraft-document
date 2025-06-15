Data Enhancement
============================================

Environment Setup and Library Import
------------------------------------

.. code-block:: python

   import scanpy as sc
   from scMethCraft.model.scmethcraft_model import *
   from scMethCraft.model.scmethcraft_trainning import *
   import scMethCraft.model.methyimp as mp
   from scMethCraft.model.utils_model import *
   from scMethCraft.function.embedding import *
   from scMethCraft.model.compute_pos import return_pos
   from scMethCraft.function.enhancement import *
   import torch

Data Loading
------------

Load the raw methylation data from an AnnData file.

.. code-block:: python

   input_path = f"../project/sample_data/genome/"
   raw_adata = sc.read(f"{input_path}/adata.h5ad")

Parameter Configuration
-----------------------

Set key parameters for the analysis pipeline.

.. code-block:: python

   cell = raw_adata.shape[0]

   seq_length = 10000
   kmer_k = 8
   work_dir = './'
   device = "cuda:2"
   dataset = "Test_dataset"
   input_path = f"../project/sample_data/genome/"

Data Preparation
----------------

Load and prepare sequence data for model processing.

.. code-block:: python

   train_onehot,train_kmer,pos = load_seq(input_path,"all_seqs.h5",False,"both")
   train_pos = return_pos(pos)
   train_state = load_state(input_path,"m_all.npy",False)
   train_data = MethyDataset(train_onehot,train_kmer,train_state,train_pos)
   cell = train_state.shape[1]
   train_dataloader = torch.utils.data.DataLoader(train_data,batch_size=64,shuffle=True,num_workers=10,pin_memory=True)
   test_dataloader = torch.utils.data.DataLoader(train_data,batch_size=64,shuffle=False,num_workers=10,pin_memory=True)
   del train_onehot,train_kmer,train_state,train_pos

Model Initialization and Loading
--------------------------------

Initialize the scMethCraft model and load pre-trained weights.

.. code-block:: python

   scMethCraft_part1 = Sequence_extraction(cell,K=kmer_k,genomic_seq_length=seq_length).to(device)
   scMethCraft_part2 = Similarity_weighting(cell,dropout_rate=0.5).to(device)
   modelpath = f"../project/sample_data/output/"
   scMethCraft_part1.load_state_dict(torch.load(f"{modelpath}/scMethCraft_part1.pth"))
   scMethCraft_part2.load_state_dict(torch.load(f"{modelpath}/scMethCraft_part2.pth"))

Data Enhancement
----------------

Generate enhanced methylation data using the trained model.

.. code-block:: python

   adata = raw_adata.copy()
   adata.X = output_enhanced_data(scMethCraft_part1,scMethCraft_part2,train_dataloader,cell,device)
