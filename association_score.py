import pandas as pd
import numpy as np


def Association_Score(source_cell_type: str, target_cell_type: str,
                      df_data: pd.DataFrame, CCC_list: pd.DataFrame,
                      abs_ccc: pd.DataFrame, threshold: float):
    # Filter source and target cell types
    filtered_CCC_list = CCC_list[(CCC_list['label_row'] == source_cell_type) &
                                 (CCC_list['label_column'] == target_cell_type)]

    # Calculate communication weights
    filtered_CCC_list['weight'] = filtered_CCC_list.apply(
        lambda row: abs_ccc.loc[row['Row'], row['Column']], axis=1)

    # Reset index
    filtered_CCC_list.index = range(len(filtered_CCC_list.index))

    # Get the corresponding sender and receiver cells
    send_cell_list = filtered_CCC_list['Row'].tolist()
    receive_cell_list = filtered_CCC_list['Column'].tolist()

    # iloc index starts from 0
    send_cell_list = [x - 1 for x in send_cell_list]
    receive_cell_list = [x - 1 for x in receive_cell_list]

    # Extract sender and receiver cell data from df_data
    send_cell_data = df_data.iloc[:, send_cell_list]
    send_cell_data.columns = range(0, len(send_cell_data.columns))
    print('send_cell_data:\n', send_cell_data)

    receive_cell_data = df_data.iloc[:, receive_cell_list]
    receive_cell_data.columns = range(0, len(receive_cell_data.columns))
    print('receive_cell_data:\n', receive_cell_data)

    # Count the number of non-zero elements in each row
    non_zero_counts_df1 = (send_cell_data != 0).sum(axis=1)
    non_zero_counts_df2 = (receive_cell_data != 0).sum(axis=1)

    # Filter rows where the number of zero elements is less than threshold * the number of cell pairs
    ccc_length = len(send_cell_list)
    filtered_indices_send = send_cell_data.index[non_zero_counts_df1 < (ccc_length * threshold)]
    filtered_indices_receive = receive_cell_data.index[non_zero_counts_df2 < (ccc_length * threshold)]

    # Drop these rows
    send_cell_data = send_cell_data.drop(index=filtered_indices_send)
    receive_cell_data = receive_cell_data.drop(index=filtered_indices_receive)

    # Reset index
    row = send_cell_data.index
    column = receive_cell_data.index
    send_cell_data.index = range(len(send_cell_data.index))
    receive_cell_data.index = range(len(receive_cell_data.index))
    print('final send_cell_data:\n', send_cell_data)
    print('final receive_cell_data:\n', receive_cell_data)

    # Create a zero matrix to store scores
    score_matrix = np.zeros((len(send_cell_data.index), len(receive_cell_data.index)))
    print('score_matrix: shape=', score_matrix.shape)

    # Calculate communication scores between cell pairs
    count = 1
    for i in range(score_matrix.shape[0]):
        print('Progress Bars: ',count,'/',score_matrix.shape[0])
        count += 1
        for j in range(score_matrix.shape[1]):
            score = 0
            for k in range(len(filtered_CCC_list.index)):
                score += filtered_CCC_list.loc[k, 'weight'] * np.log(
                    send_cell_data.loc[i, k] * receive_cell_data.loc[j, k] + 1)
            score /= len(filtered_CCC_list.index)
            score_matrix[i][j] = score

    # Generate the score matrix
    Absolute_Score = pd.DataFrame(score_matrix, index=row, columns=column)

    # Calculate weighted quantile means for each row and column
    row_means = 0.5 * Absolute_Score.quantile(0.5, axis=1) + \
                0.25 * (Absolute_Score.quantile(0.25, axis=1) + Absolute_Score.quantile(0.75, axis=1))

    column_means = 0.5 * Absolute_Score.quantile(0.5) + \
                   0.25 * (Absolute_Score.quantile(0.25) + Absolute_Score.quantile(0.75))

    # Calculate relative scores
    Relative_Score = Absolute_Score.copy()
    for i in range(Relative_Score.shape[0]):
        for j in range(Relative_Score.shape[1]):
            Relative_Score.iloc[i, j] = Relative_Score.iloc[i, j] - (row_means[i] + column_means[j])

    # Sort the relative scores and reset the index
    sorted_df = Relative_Score.stack(dropna=False).sort_values(ascending=False).reset_index()
    sorted_df.columns = ['Source_gene', 'Target_gene', 'RAS']

    return Absolute_Score, Relative_Score, sorted_df



