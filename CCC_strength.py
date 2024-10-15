import pandas as pd


def CellType_Count(label_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the count of each label.

    :param label_df: A DataFrame containing labels, must include a 'label' column.
    :param communication_df: A DataFrame containing communication relationships.
    :return: A DataFrame containing the count of each label, formatted as ['label', 'count'].
    """
    # Calculate the number of occurrences for each label
    label_counts = label_df['label'].value_counts().reset_index()
    label_counts.columns = ['label', 'count']
    return label_counts


def CWCT(label_counts: pd.DataFrame, communication_df: pd.DataFrame) -> pd.DataFrame:
    """
    Populate the communication matrix and calculate the proportion of each element
    based on the product of label counts.

    :param label_counts: A DataFrame containing the count of each label, 'label' column should be the index.
    :param communication_df: A DataFrame containing communication relationships with 'label_row' and 'label_column'.
    :return: A new communication matrix with communication weights between cell types .
    """
    # Ensure 'label' is set as the index in label_counts
    if 'label' in label_counts.columns:
        label_counts.set_index('label', inplace=True)

    # Get unique labels from label_counts and create a zero matrix for communication
    label_indices = label_counts.index.tolist()
    ccc_matrix = pd.DataFrame(0, index=label_indices, columns=label_indices)

    # Populate the communication matrix based on the relationships in communication_df
    for _, row in communication_df.iterrows():
        a = row["label_row"]
        b = row["label_column"]
        if a in ccc_matrix.index and b in ccc_matrix.columns:
            ccc_matrix.at[a, b] += 1

    # Create a new DataFrame to store the calculated communication weights between cell types
    ccc_matrix_new = pd.DataFrame(0.0, index=ccc_matrix.index, columns=ccc_matrix.columns)

    # Calculate the proportion for each element
    for i, row in ccc_matrix.iterrows():
        for j, value in row.items():
            # Get the product of the corresponding count values
            count_product = label_counts.loc[i, 'count'] * label_counts.loc[j, 'count']
            # Compute the new value
            new_value = value / count_product if count_product != 0 else 0.0
            ccc_matrix_new.loc[i, j] = new_value

    return ccc_matrix_new



