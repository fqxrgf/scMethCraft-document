Model Training
========================

.. contents::
   :local:
   :depth: 2

Environment Setup and Library Import
------------------------------------

Import necessary modules and set up the Python environment for methylation analysis.

.. code-block:: python

   from scMethCraft.model.scmethcraft_model import *
   from scMethCraft.model.utils_model import *
   from scMethCraft.model.scmethcraft_trainning import *
   from scMethCraft.model.compute_pos import return_pos

   import scMethCraft.benchmark.methyimp as mp
   import pandas as pd
   import sklearn.metrics as metrics
   import scanpy as sc
   import anndata as ad
   from scipy.special import expit
   import sys
   import numpy
   import random

Seed Initialization
-------------------

Set random seeds for reproducibility across different runs.

.. code-block:: python

   def seed_everything(seed=11):
       random.seed(seed)
       np.random.seed(seed)
       torch.manual_seed(seed)
       torch.cuda.manual_seed(seed)
       torch.cuda.manual_seed_all(seed)
       torch.backends.cudnn.deterministic = True
       torch.backends.cudnn.benchmark = False

   seed_everything()

Parameter Configuration
-----------------------

Set key parameters for the model and analysis.

.. code-block:: python

   seq_length = 10000
   kmer_k = 8
   work_dir = './'
   device = "cuda:2"
   dataset = "Test_dataset"
   input_path = f"../project/sample_data/genome/"

Data Loading and Preparation
----------------------------

Load and prepare sequence data for model training.

.. code-block:: python

   train_onehot, train_kmer, pos = load_seq(input_path, "all_seqs.h5", False, "both")
   train_pos = return_pos(pos)
   train_state = load_state(input_path, "m_all.npy", False)
   train_data = MethyDataset(train_onehot, train_kmer, train_state, train_pos)
   cell = train_state.shape[1]
   train_dataloader = torch.utils.data.DataLoader(train_data, batch_size=64, shuffle=True, num_workers=10, pin_memory=True)
   test_dataloader = torch.utils.data.DataLoader(train_data, batch_size=64, shuffle=False, num_workers=10, pin_memory=True)
   del train_onehot, train_kmer, train_state, train_pos

Model Initialization
---------------------

Initialize the two-part scMethCraft model.

.. code-block:: python

   scMethCraft_part1 = Sequence_extraction(cell, K=kmer_k, genomic_seq_length=seq_length).to(device)
   scMethCraft_part2 = Similarity_weighting(cell, dropout_rate=0.5).to(device)

Training Setup
--------------

Configure loss function, optimizer, and training parameters.

.. code-block:: python

   loss_fn = nn.BCEWithLogitsLoss().to(device)
   learning_rate = 1e-2
   optimizer1 = torch.optim.Adam(scMethCraft_part1.parameters(), lr=learning_rate)
   optimizer2 = torch.optim.Adam(scMethCraft_part2.parameters(), lr=learning_rate)
   epoch = 5

Model Training
--------------

Execute the training loop for 100 epochs.

.. code-block:: python

   train_losses1 = []
   train_losses2 = []

   for i in range(epoch):
       print("-------epoch  {} -------".format(i+1))
       MethyBasset_part1.train()
       train_loss1 = 0
       train_loss2 = 0
       for step, [onehot, targets, kmer, pos] in enumerate(train_dataloader):
           onehot = onehot.to(device).float()
           kmer = kmer.to(device).float()
           pos = pos.to(device).float()
           targets = targets.to(device).float()

           outputs_part1 = MethyBasset_part1(onehot, kmer, pos).float()
           fla = nn.Flatten(0)

           reconstructed_matrix = torch.sigmoid(outputs_part1)
           imputed_matrix = reconstructed_matrix * (torch.isnan(targets)).int() + torch.nan_to_num(targets)

           outputs_part1 = fla(outputs_part1)
           targets = fla(targets)

           is_loss = ~torch.isnan(targets)

           loss1 = loss_fn(outputs_part1[is_loss], targets[is_loss])
           train_loss1 += loss1.item()

           optimizer1.zero_grad()
           loss1.backward()
           optimizer1.step()

           outputs_part2 = MethyBasset_part2(imputed_matrix.detach())
           outputs_part2 = fla(outputs_part2)
           loss2 = loss_fn(outputs_part2[is_loss], targets[is_loss])
           train_loss2 += loss2.item()

           optimizer2.zero_grad()
           loss2.backward()
           optimizer2.step()

       train_losses1.append(train_loss1)
       train_losses2.append(train_loss2)
       print(f"Loss1: {train_loss1}, Loss2: {train_loss2}")

Execute Model Saving
--------------------

Run the output function to save the trained model.

.. code-block:: python

   output_model(scMethCraft_part1, scMethCraft_part2, savepath=f"../project/sample_data/output/")
