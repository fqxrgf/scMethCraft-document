# scMethCraft: A versatile toolkit for analysis of single-cell DNA methylation data
![](scMethCraft_figure.jpg)

## Installation  

It's prefered to create a new environment for scMethCraft

```
conda create -n scMethCraft python==3.9
conda activate scMethCraft
```

You can install scMethCraft from GitHub via
	
	git clone git://github.com/BioX-NKU/scMethCraft.git
	cd scMethCraft
	python setup.py install

The dependencies will be automatically installed along with scMethCraft   

This process will take approximately 5 to 10 minutes, depending on the user's computer device and internet connectivition.



## Quick Start

Using scMethCraft to analyze single-cell DNA methylation data involves three main steps:  
* **First**, the methylation site-level data is processed into a **cell-by-region matrix** and subsequently preprocessed into sequence-based inputs compatible with scMethCraft.  
* **Second**, a **scMethCraft model** is trained on the processed data, and the trained model parameters are saved, serving as the foundation for downstream applications.  
* **Finally**, the data and trained model parameters are provided to scMethCraft functions to perform various downstream analyses.


### Environment
Anndata object is a Python object/container designed to store single-cell data in Python packege [**anndata**](https://anndata.readthedocs.io/en/latest/) which is seamlessly integrated with [**scanpy**](https://scanpy.readthedocs.io/en/stable/), a widely-used Python library for single-cell data analysis. scMethCraft is built upon [PyTorch](https://pytorch.org/), an open-source deep learning framework in Python.


### 1. Training Data Generation

The main workflow for training data generation is as follows:  
First, a **cell-by-region matrix** is constructed from site-level methylation data. This step is primarily adapted from **EpiScanpy**, with extensions and modifications to enhance its applicability.  
Next, the matrix is associated with the corresponding metadata. Subsequently, regions that are unsuitable for training (e.g., with excessive missing values) are filtered out.  
Finally, the corresponding sequence information is extracted from the reference genome, partially leveraging functionality from **basset**.

Key functions involved:

| Function        | Description                                      |
|----------------|--------------------------------------------------|
| `count_matrix` | Generate a matrix from site-level methylation data |
| `assemble_meta`| Associate the matrix with metadata               |
| `filter_region`| Filter regions for training                      |
| `obtain_seq`   | Retrieve the corresponding sequence from the reference genome |

A tutorial is provided to demonstrate these functions to users:

#### [Training Data Generation](./tutorial/tutorial_create_training_datasets.ipynb)
Sample data can be downloaded from:https://drive.google.com/drive/folders/1vaLj5UoJ5wi46ZJAzbwaEkRZU_UEa_k0?usp=drive_link

### 2. Model Training

The scMethCraft model consists of two modules: a **sequence feature extraction module** and a **similarity weighting module**.  
To train scMethCraft, the following steps are required: read the preprocessed data from the previous step, instantiate the scMethCraft model, and specify two iterators to update each module independently.  
After a certain number of epochs, the trained model parameters are saved.

Key functions involved:

| Function               | Description                                  |
|------------------------|----------------------------------------------|
| `load_seq`             | Load sequence data                           |
| `MethyDataset`         | Create the data structure for training       |
| `Sequence_extraction`  | Sequence feature extraction module of scMethCraft |
| `Similarity_weighting` | Similarity weighting module of scMethCraft   |
| `output_model`         | Output the trained model                     |

A tutorial is provided to demonstrate these functions to users:

#### [Model Training](./tutorial/tutorial_model_training.ipynb)

### 3. Downstream Analysis

scMethCraft can be applied to various downstream analyses, including dimensionality reduction, data enhancement, batch integration, cell type annotation, and DMR (differentially methylated region) identification.

The following tutorials are provided for users:

#### [Dimensionality Reduction](./tutorial/tutorial_cell_embedding.ipynb)
#### [Data Enhancement](./tutorial/tutorial_cell_embedding.ipynb)
#### [Batch Integration](./tutorial/tutorial_batch_integration.ipynb)
#### [Cell type annotation](./tutorial/tutorial_annotation.ipynb)


