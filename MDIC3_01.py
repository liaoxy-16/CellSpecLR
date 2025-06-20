import numpy as np
import pandas as pd

def compute_binary_ccc(GRN_file, expression_file, CCC_output_file, CCC_binaryzation_output_file, threshold=0.25):
    print('Reading expression data...')
    expression_data = pd.read_csv(expression_file, delimiter='\t', index_col=0).values
    E = np.matrix(expression_data)

    print('Reading GRN data...')
    # Each element represents the regulation of the row index gene by the column index gene
    GRN_data = pd.read_csv(GRN_file, delimiter='\t', index_col=0).values
    G = np.matrix(GRN_data)

    print('Computing CCC matrix...')
    T = np.dot(G, E)
    Tp = np.linalg.pinv(T)
    CCC_result = np.dot(Tp, E)
    print('CCC_result:\n',CCC_result)
    CCC_result_df = pd.DataFrame(CCC_result, index=range(CCC_result.shape[0]), columns=range(CCC_result.shape[1]))
    absolute_max = CCC_result_df.abs().max().max()
    CCC_result_scaled = CCC_result_df / absolute_max
    CCC_result_scaled.to_csv(CCC_output_file)

    abs_matrix = np.abs(CCC_result_df.to_numpy(dtype='float32'))
    flat = abs_matrix.flatten()
    #print('len(flat)=',len(flat))
    top_n = int(len(flat) * threshold)
    #print('top_n=',top_n)
    top_indices = np.argsort(flat)[::-1][:top_n]

    bin_matrix = np.zeros_like(abs_matrix)
    np.put(bin_matrix, top_indices, 1)
    #np.fill_diagonal(bin_matrix, 0)

    bin_df = pd.DataFrame(bin_matrix, index=range(bin_matrix.shape[0]), columns=range(bin_matrix.shape[1]))
    bin_df.to_csv(CCC_binaryzation_output_file)

    print('Done.')
    return CCC_result_scaled,bin_df


# Example
compute_binary_ccc(
    GRN_file='./data/example_GRN.txt',
    expression_file='./data/example_data.txt',
    CCC_output_file='./data/CCC_Network.csv',
    CCC_binaryzation_output_file='./data/CCC_Network_binaryzation.csv',
    threshold=0.25
)
