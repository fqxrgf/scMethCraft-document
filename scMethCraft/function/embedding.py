import numpy as np
from ..postprecessing.similarity_norm import GCN_norm

def output_embedding(raw_adata,scMethCraft_part1,scMethCraft_part2):
    adata = raw_adata.copy()
    similarity_matrix = scMethCraft_part2.state_dict()["SimilarityLayer1.similarity_matrix"].cpu().numpy()
    similarity_matrix = np.abs(similarity_matrix+similarity_matrix.T)/2
    adata.obsm["Similarity_matrix"] = similarity_matrix

    final_embedding = scMethCraft_part1.state_dict()["final.weight"].cpu().numpy()
    adata.obsm["Loading_matrix"] = final_embedding
    adata.obsm["Cell_embedding"] = GCN_norm(similarity_matrix)@final_embedding
    return adata