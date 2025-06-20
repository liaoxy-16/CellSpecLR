import pandas as pd
import time
import numpy as np
import os

def compute_penalty(
    expression_file: str,
    cell_label_file: str,
    comm_weight_file: str,
    ccc_pair_file: str,
    output_dir: str,
    threshold: float = 0.3,
    top_n_pairs: int = 5
):
    begin = time.asctime()
    print("Begin time: " + begin)

    expression_df = pd.read_csv(expression_file, index_col=0, sep='\t')
    cell_type_df = pd.read_csv(cell_label_file, sep='\t')
    comm_weight_df = pd.read_csv(comm_weight_file, index_col=0).abs()

    expression_df.columns = expression_df.columns.astype(int)
    cell_type_df['index'] = cell_type_df['index'].astype(int)
    comm_weight_df.index = comm_weight_df.index.astype(int)
    comm_weight_df.columns = comm_weight_df.columns.astype(int)

    print('expression_df:\n', expression_df)
    print('cell_type_df:\n', cell_type_df)
    print('comm_weight_df:\n', comm_weight_df)

    CellType_CCC_list_df = pd.read_csv(ccc_pair_file)
    print('CellType_CCC_list_df:\n', CellType_CCC_list_df)

    CellType_CCC_list_df = CellType_CCC_list_df[
        CellType_CCC_list_df['Out_CellType'] != CellType_CCC_list_df['In_CellType']
    ]
    print('Filtered CellType_CCC_list_df:\n', CellType_CCC_list_df)

    CellType_CCC_list_df = CellType_CCC_list_df.iloc[0:top_n_pairs,:]
    print('Sliced CellType_CCC_list_df:\n', CellType_CCC_list_df)

    os.makedirs(output_dir, exist_ok=True)

    for index, row in CellType_CCC_list_df.iterrows():
        out_cell_type = row['Out_CellType']
        in_cell_type = row['In_CellType']

        print(f"Processing {out_cell_type} -> {in_cell_type}")

        source_cells = cell_type_df[cell_type_df['label'] == out_cell_type]['index'].values
        target_cells = cell_type_df[cell_type_df['label'] == in_cell_type]['index'].values

        ccc_list = []
        for source_cell in source_cells:
            for target_cell in target_cells:
                weight = comm_weight_df.at[source_cell, target_cell]
                ccc_list.append({'source': source_cell, 'target': target_cell, 'weight': weight})

        ccc_df = pd.DataFrame(ccc_list)
        #print('ccc_df:\n', ccc_df)

        ccc_num = ccc_df.shape[0]
        #print('ccc_num:', ccc_num)

        ccc_df = ccc_df[ccc_df['weight'] != 0]
        #print('Filtered ccc_df:\n', ccc_df)

        ligand_df = pd.DataFrame(index=expression_df.index, columns=ccc_df['source'].values)
        receptor_df = pd.DataFrame(index=expression_df.index, columns=ccc_df['target'].values)

        ligand_values = expression_df[ccc_df['source'].values].values
        receptor_values = expression_df[ccc_df['target'].values].values

        ligand_df.iloc[:, :] = ligand_values
        receptor_df.iloc[:, :] = receptor_values

        #print('ligand_df:\n', ligand_df)
        #print('receptor_df:\n', receptor_df)

        ligand_non_zero_count = ligand_df.apply(lambda row: (row != 0).sum(), axis=1)
        ligand_row_threshold = int(ligand_df.shape[1] * threshold)
        ligand_df_filtered = ligand_df[ligand_non_zero_count > ligand_row_threshold]

        receptor_non_zero_count = receptor_df.apply(lambda row: (row != 0).sum(), axis=1)
        receptor_row_threshold = int(receptor_df.shape[1] * threshold)
        receptor_df_filtered = receptor_df[receptor_non_zero_count > receptor_row_threshold]

        #print('ligand_df_filtered:\n', ligand_df_filtered)
        #print('receptor_df_filtered:\n', receptor_df_filtered)

        penalty_df = pd.DataFrame(index=ligand_df_filtered.index, columns=receptor_df_filtered.index)

        for i, ligand_gene in enumerate(ligand_df_filtered.index):
            #print('i=', i)
            for j, receptor_gene in enumerate(receptor_df_filtered.index):
                total_weighted_sum1 = 0
                total_weighted_sum2 = 0

                ligand_values = ligand_df_filtered.loc[ligand_gene].values.astype(float)
                receptor_values = receptor_df_filtered.loc[receptor_gene].values.astype(float)

                if np.any(ligand_values == 0) and np.any(receptor_values != 0):
                    condition = (ligand_values == 0) & (receptor_values != 0)
                    penalty_values1 = receptor_values[condition] * ccc_df['weight'].values[condition]
                    total_weighted_sum1 = penalty_values1.sum()

                if np.any(ligand_values != 0) and np.any(receptor_values == 0):
                    condition = (ligand_values != 0) & (receptor_values == 0)
                    penalty_values2 = ligand_values[condition] * ccc_df['weight'].values[condition]
                    total_weighted_sum2 = penalty_values2.sum()

                penalty = (total_weighted_sum1 + total_weighted_sum2) / ccc_num
                penalty_df.loc[ligand_gene, receptor_gene] = penalty

        print("penalty_df:\n", penalty_df)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        output_filename = f'{output_dir}{out_cell_type}_{in_cell_type}_Penalty.csv'
        penalty_df.to_csv(output_filename)
        output_filename_txt = f'{output_dir}{out_cell_type}_{in_cell_type}_Penalty.txt'
        penalty_df.to_csv(output_filename_txt, sep='\t')

        print(f"Results saved as {output_filename}")

    print("End time: " + time.asctime())


compute_penalty(
    expression_file='./data/example_data.txt',
    cell_label_file='./data/example_label.txt',
    comm_weight_file='./data/CCC_Network_binaryzation.csv',
    ccc_pair_file='./data/Cell_Type_Comm_Weight_sort.csv',
    output_dir='./Results_Penalty/',
    threshold=0.3,
    top_n_pairs=5
)
