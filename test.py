from CellSpecLR import *
import numpy as np
import pandas as pd


if __name__ == '__main__':
# ======================================================================================
# 1.MDIC3

    # Read the CSV file with the first row as the column index and the first column as the row index
    # index: gene   column: cell
    expression_data = pd.read_csv('./data/data.txt', delimiter='\t', index_col=0)
    print('expression matrix:\n', expression_data)
    expression_data = expression_data.values

    E = np.matrix(expression_data)
    print('E:\n',E)

    # ------------------------------------------------------------------

    # Each element represents the regulation of the row index gene by the column index gene
    GRN_data = pd.read_csv('./data/data_GRN.txt', delimiter='\t', index_col=0)
    print('GRN matrix:\n', GRN_data)
    GRN_data = GRN_data.values

    # Transposed data: the regulation of the row index gene by the column index gene
    GRN_data = np.transpose(GRN_data)
    print('GRN transpose:\n', GRN_data)

    G = np.matrix(GRN_data)
    print('GRN matrix:\n',G)

    CCC_result = MDIC3(G, E)

    # output
    np.savetxt('./data/CCC_result.txt', CCC_result, delimiter='\t')
    print('CCC_result:\n',CCC_result)

    df = pd.DataFrame(CCC_result, index=range(1, CCC_result.shape[0] + 1),
                      columns=range(1, CCC_result.shape[1] + 1))
    df.to_csv('./data/CCC_result.csv')
    print('CCC_result df:\n', df)


# ======================================================================================
# 2.filter_CCC_matrix
# Load the CCC_result matrix from a CSV file
    CCC_result = pd.read_csv('./data/CCC_result.csv', index_col=0).to_numpy(dtype='float32')

    # Set t = 0.25 to process the top 25% largest elements by absolute value
    result_matrix = filter_CCC_matrix(CCC_result, 0.25, './data/filter_CCC_result.csv')
    print("The processed matrix:\n", result_matrix)

# ======================================================================================
# 3.CCC_list
    df = pd.read_csv('./data/filter_CCC_result.csv', index_col=0)
    df_labels = pd.read_csv('./data/data_label.txt', delimiter='\t')

    # Call the function and get the result
    result_df = process_cell_communication(df, df_labels)

    # Save the result to a file
    result_df.to_csv('./data/CCC_list.txt', index=False, sep='\t')
    print('Processing complete, result saved to ./data/CCC_list.txt')

# ======================================================================================
# 4.CCC_strength
    # Example usage
    # Read the label data and communication relationship data
    df1 = pd.read_csv('./data/data_label.txt', delimiter='\t')
    df2 = pd.read_csv('./data/CCC_list.txt', delimiter='\t')

    # Calculate the label counts
    label_counts_df = CellType_Count(df1)
    print('label_counts_df:\n', label_counts_df)

    # Calculate the proportion-adjusted communication matrix using CWCT
    CWCT_matrix = CWCT(label_counts_df, df2)
    print('CWCT_matrix:\n', CWCT_matrix)

    # Save results to CSV files
    label_counts_df.to_csv('./data/CellType_counts.csv', header=True)
    CWCT_matrix.to_csv("./data/CWCT_matrix.csv")

# ======================================================================================
# 5.Rank_sum_test
    # Example usage
    expression_data = pd.read_csv('./data/data.txt', delimiter='\t', index_col=0)
    cell_types = pd.read_csv('./data/data_label.txt', delimiter='\t')

    # Set p_value_threshold to 0.01 for more stringent filtering (optional)
    ranksums_gene_list, ranksums_matrix = Rank_sum_test(expression_data, cell_types, p_value_threshold=0.05)

    ranksums_gene_list.to_csv('./data/ranksums_gene_list.txt', index=False, sep='\t')
    ranksums_matrix.to_csv('./data/ranksums_data.txt', sep='\t')

# ======================================================================================
# 6.association_score
    # Read input files
    print('Reading cell communication label data...')
    CCC_list = pd.read_csv('./data/CCC_list.txt', delimiter='\t')
    print('CCC_list:\n', CCC_list)

    # Read communication result data
    abs_ccc = pd.read_csv('./data/CCC_result.csv', index_col=0)
    print('Reading communication result data...')
    abs_ccc = abs_ccc.abs()
    abs_ccc.index = range(1, len(abs_ccc.index) + 1)
    abs_ccc.columns = range(1, len(abs_ccc.columns) + 1)
    print('abs_ccc:\n', abs_ccc)

    # Read gene expression matrix
    print('Reading gene expression matrix data...')
    df_data = pd.read_csv('./data/ranksums_data.txt', delimiter='\t', index_col=0)
    df_data.columns = range(1, len(df_data.columns) + 1)
    df_data.index.name = None  # Remove index name
    print('expression matrix:\n', df_data)


    print('Calling Association_Score function...')
    abs_score, rel_score, sorted_scores = Association_Score(
        source_cell_type='cDC1',
        target_cell_type='NKT',
        df_data=df_data,
        CCC_list=CCC_list,
        abs_ccc=abs_ccc,
        threshold=0.5
    )

    # Save results to CSV files
    print('Saving Absolute_Score to CSV file...')
    abs_score.to_csv('./data/cDC1-NKT_AAS.csv')

    print('Saving Relative_Score to CSV file...')
    rel_score.to_csv('./data/cDC1-NKT_RAS.csv')

    print('Saving sorted relative scores to CSV file...')
    sorted_scores.to_csv('./data/cDC1-NKT_RAS_sort.txt', index=False, sep='\t')

    print('Saving completed!')