import numpy as np

def GCN_norm(similarity_matrix_,alpha = 1):
    similarity_matrix_ = similarity_matrix_*(1-np.eye(similarity_matrix_.shape[0]))
    row_sum = similarity_matrix_.sum(axis = 1)
    D_minus_1_2 = np.diag(np.power(row_sum+1, -0.5))
    similarity_matrix_ = D_minus_1_2@similarity_matrix_@D_minus_1_2
    similarity_matrix_ = similarity_matrix_+alpha*np.eye(similarity_matrix_.shape[0])
    return similarity_matrix_

def Standard_norm(similarity_matrix_,alpha = 1):
    similarity_matrix_ = similarity_matrix_*(1-np.eye(similarity_matrix_.shape[0]))
    row_sum = similarity_matrix_.sum(axis = 1)
    similarity_matrix_ = (similarity_matrix_/row_sum).T
    similarity_matrix_ = similarity_matrix_+alpha*np.eye(similarity_matrix_.shape[0])
    return similarity_matrix_
