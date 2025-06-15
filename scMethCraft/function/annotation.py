import numpy as np
from sklearn.preprocessing import LabelEncoder

def cell_annotation(adata,train_index,test_index,type_column = "MajorType"):
    label_encoder = LabelEncoder()
    
    adata.obs[type_column] = label_encoder.fit_transform(adata.obs[type_column])
    train_data = adata.obsm["Similarity_matrix"][train_index]
    test_data = adata.obsm["Similarity_matrix"][test_index]

    train_label = adata.obs[type_column][train_index]
    test_label = adata.obs[type_column][test_index]

    test_adata = adata[test_index].copy()
    transfer_matrix = test_adata.obsm["Similarity_matrix"][:,train_index]

    k = 10
    similarity_matrix = transfer_matrix
    n_targets, n_sources = similarity_matrix.shape
    knn_indices = np.zeros((n_targets, k), dtype=int)
    knn_distances = np.zeros((n_targets, k))

    for i in range(n_targets):
            similarities = similarity_matrix[i, :]
            partitioned_indices = np.argpartition(-similarities, k)[:k]
            sorted_indices = partitioned_indices[np.argsort(-similarities[partitioned_indices])]
            knn_indices[i, :] = sorted_indices

    predict = np.zeros((n_targets, 1), dtype=int)
    for i in range(n_targets):
        labels = train_label[knn_indices[i]]
        weights = similarities[sorted_indices]
        weighted_votes = np.bincount(labels, weights=weights)
        predict[i, :] = np.argmax(weighted_votes)
    adata.obs[type_column] = label_encoder.inverse_transform(adata.obs[type_column])
    
    return label_encoder.inverse_transform(predict)