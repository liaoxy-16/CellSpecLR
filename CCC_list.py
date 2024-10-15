import pandas as pd

def process_cell_communication(df: pd.DataFrame, df_labels: pd.DataFrame) -> pd.DataFrame:
    """
    Process cell communication data, calculate communication pairs, and merge label information.

    Parameters:
    - df: DataFrame, cell communication matrix with values 0 or 1
    - df_labels: DataFrame, contains cell label information, must include 'index' and 'label' columns

    Returns:
    - merged_df: DataFrame, contains cell communication pairs and their corresponding labels
    """
    # Calculate the number of elements with a value of 1
    count_1 = df[df == 1].count().sum()
    print('Number of elements with a value of 1:', count_1)

    # Find the row and column indices of elements with a value of 1
    indices_1 = df[df == 1].stack().reset_index()
    indices_1.columns = ['Row', 'Column', 'Communication']

    # Print the number of rows (i.e., the number of communication pairs)
    num_rows = indices_1.shape[0]
    #print("Number of rows:", num_rows)


    # Convert 'Row' and 'Column' to integer type to ensure consistency for merging
    indices_1['Row'] = indices_1['Row'].astype(int)
    indices_1['Column'] = indices_1['Column'].astype(int)
    df_labels['index'] = df_labels['index'].astype(int)


    # Merge with label data to generate the label for the sender cell (label_row)
    merged_df = pd.merge(indices_1, df_labels, left_on='Row', right_on='index')
    merged_df.rename(columns={'label': 'label_row'}, inplace=True)

    # Continue merging to generate the label for the receiver cell (label_column)
    merged_df = pd.merge(merged_df, df_labels, left_on='Column', right_on='index')
    merged_df.rename(columns={'label': 'label_column'}, inplace=True)
    print("merged_df\n:", merged_df)

    # Return the merged result, keeping only the necessary columns
    return merged_df[['Row', 'Column', 'label_row', 'label_column']]


