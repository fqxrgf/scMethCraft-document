Cell type annotation
=================================================================

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
   from scMethCraft.function.annotation import *
   import scipy

Data Loading
------------

Load the batch integrated data.

.. code-block:: python

   adata = sc.read("/home/sccasimp/data/methyimp/dataset/allchr/mix_total/total_newmodel_0301.h5ad")

Cell Type Annotation
--------------------

Annotate cell types based on training and testing batches.

.. code-block:: python

   train_index = adata.obs["batch"] == "batch1"
   test_index = adata.obs["batch"] != "batch1"

   predict = cell_annotation(adata, train_index, test_index, type_column="MajorType")

.. code-block:: python

   predict
