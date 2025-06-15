import numpy as np
from ..postprecessing.similarity_norm import GCN_norm

def batch_correction_similarity(batch_label: list, alpha = 0.3) -> np.array:
    batch_correction_matrix = np.tile(batch_label,(len(batch_label),1))
    batch_correction_matrix = (batch_correction_matrix == np.array(batch_label.values).reshape(-1,1))
    x = 1/(1-alpha)
    batch_correction_matrix = (batch_correction_matrix.astype(int)*-1+x)/x
    return batch_correction_matrix

def output_batch_integration(raw_adata,scMethCraft_part1,scMethCraft_part2):
    adata = raw_adata.copy()
    correction_matrix = batch_correction_similarity(raw_adata.obs["batch"],0.2)
    similarity_matrix = scMethCraft_part2.state_dict()["SimilarityLayer1.similarity_matrix"].cpu().numpy()
    similarity_matrix = np.abs(similarity_matrix+similarity_matrix.T)/2
    adata.obsm["Similarity_matrix"] = similarity_matrix
    final_embedding = scMethCraft_part1.state_dict()["final.weight"].cpu().numpy()
    adata.obsm["Loading_matrix"] = final_embedding
    adata.obsm["Cell_embedding"] = (GCN_norm(similarity_matrix,0.8)*correction_matrix)@final_embedding
    return adata