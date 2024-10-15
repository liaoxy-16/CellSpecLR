import pandas as pd
from scipy.stats import ranksums

def Rank_sum_test(expression_data: pd.DataFrame, cell_types: pd.DataFrame, p_value_threshold: float = 0.05):
    # Set the column names of the expression matrix to sample numbers starting from 1
    expression_data.columns = range(1, len(expression_data.columns) + 1)

    # Transpose the expression matrix and merge the cell type information into the gene expression matrix
    merged_data = expression_data.T
    merged_data['label'] = cell_types.set_index('index')['label']

    # Perform Wilcoxon rank-sum test and filter out significantly different genes
    significantly_different_genes = []
    for gene in expression_data.index:
        expression_data_temp = merged_data[gene]
        for cell_type in cell_types['label'].unique():
            current_type_expression = expression_data_temp[merged_data['label'] == cell_type]
            other_types_expression = expression_data_temp[merged_data['label'] != cell_type]

            # Perform Wilcoxon rank-sum test
            stat, p_value = ranksums(current_type_expression, other_types_expression)

            # If p-value is less than the specified threshold, record the significantly different gene
            if p_value < p_value_threshold:
                significantly_different_genes.append((gene, cell_type, p_value))

    # Convert the significantly different genes and their p-values to a DataFrame
    diff_genes_df = pd.DataFrame(significantly_different_genes, columns=['gene', 'label', 'p_value'])

    # Get the unique list of significantly different genes
    unique_genes = diff_genes_df['gene'].drop_duplicates().reset_index(drop=True)
    ranksums_gene_list = pd.DataFrame(unique_genes, columns=['gene'])

    # Generate a new gene expression matrix based on the filtered genes
    ranksums_matrix = expression_data.loc[ranksums_gene_list['gene'], :]

    # Return the list of significantly different genes and the expression matrix
    return ranksums_gene_list, ranksums_matrix

