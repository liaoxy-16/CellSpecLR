
import pandas as pd
import numpy as np

def compute_celltype_communication(adj_file, label_file, output_dir1, output_dir2, output_dir3):
    # Read the cell label file and calculate the number of each cell type
    cell_to_label_df = pd.read_csv(label_file, delimiter='\t')
    cell_to_label_df.index = range(len(cell_to_label_df.index))
    label_counts = cell_to_label_df['label'].value_counts().reset_index()
    label_counts.columns = ['label', 'count']
    label_counts.to_csv(output_dir1, header=True, index=False)


    # Read the cell-cell communication adjacency matrix
    adj_matrix = pd.read_csv(adj_file, index_col=0)
    adj_matrix.index = adj_matrix.index.astype(int)
    adj_matrix.columns = adj_matrix.columns.astype(int)

    label_df = label_counts.copy()
    labels = label_df['label'].unique()
    label_count = label_df.set_index('label')['count']

    print('labels:\n', labels)
    print('label_count:\n', label_count)

    cell_to_label = pd.read_csv(label_file, delimiter='\t', index_col=0)['label']
    cell_to_label = cell_to_label.to_dict()
    print('cell_to_label:\n', cell_to_label)

    new_df = pd.DataFrame(0.0, index=labels, columns=labels)

    cell_labels = adj_matrix.index.map(cell_to_label)
    print('cell_labels:\n', cell_labels)

    for i, cell_i in enumerate(adj_matrix.index):
        # print('i=',i,' cell_i=',cell_i)
        for j, cell_j in enumerate(adj_matrix.columns):
            if adj_matrix.iloc[i, j] == 1:
                label_i = cell_labels[cell_i]
                label_j = cell_labels[cell_j]
                new_df.loc[label_i, label_j] += 1

    print('new_df:\n', new_df)

    normalized_df = new_df.copy()
    for label_i in labels:
        for label_j in labels:
            count_i = label_count.get(label_i, 1)
            count_j = label_count.get(label_j, 1)
            normalized_df.loc[label_i, label_j] = new_df.loc[label_i, label_j] / (count_i * count_j)
    print('normalized_df:\n', normalized_df)

    min_val = normalized_df.min().min()
    max_val = normalized_df.max().max()
    normalized = (normalized_df - min_val) / (max_val - min_val)
    normalized.to_csv(output_dir2)
    print('normalized:\n', normalized)

    long_df = normalized.stack().reset_index()
    long_df.columns = ['Out_CellType', 'In_CellType', 'Weight']
    sorted_df = long_df.sort_values(by='Weight', ascending=False)
    sorted_df.to_csv(output_dir3, index=False)

    print('sorted_df :\n', sorted_df )



compute_celltype_communication(
    adj_file='./data/CCC_Network_binaryzation.csv',
    label_file='./data/example_label.txt',
    output_dir1='./data/Cell_Type_count.csv',
    output_dir2='./data/Cell_Type_Comm_Weight.csv',
    output_dir3='./data/Cell_Type_Comm_Weight_sort.csv'
)

