import pandas as pd
import numpy as np
import time
import gc
import os

def compute_Reward(
    expression_file,
    cell_label_file,
    comm_weight_file,
    ccc_pair_file,
    output_dir,
    threshold=0.2,
    top_n_pairs=5
):
    print("Begin time: " + time.asctime())

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

    #CellType_CCC_list_df = CellType_CCC_list_df[CellType_CCC_list_df['Out_CellType'] != CellType_CCC_list_df['In_CellType']]
    CellType_CCC_list_df = CellType_CCC_list_df.iloc[0:top_n_pairs, :]
    print('Filtered CellType_CCC_list_df:\n', CellType_CCC_list_df)

    for index, row in CellType_CCC_list_df.iterrows():
        out_cell_type = row['Out_CellType']
        in_cell_type = row['In_CellType']

        print(f"Processing {out_cell_type} -> {in_cell_type}")

        source_cells = cell_type_df[cell_type_df['label'] == out_cell_type]['index'].values
        target_cells = cell_type_df[cell_type_df['label'] == in_cell_type]['index'].values

        ccc_list = [
            {'source': sc, 'target': tc, 'weight': comm_weight_df.at[sc, tc]}
            for sc in source_cells for tc in target_cells
        ]
        ccc_df = pd.DataFrame(ccc_list)
        #print('ccc_df:\n', ccc_df)

        ccc_num = ccc_df.shape[0]
        #print('ccc_num:', ccc_num)

        ccc_df = ccc_df[ccc_df['weight'] != 0]
        #print('Filtered ccc_df:\n', ccc_df)

        ligand_df = pd.DataFrame(index=expression_df.index, columns=ccc_df['source'].values)
        receptor_df = pd.DataFrame(index=expression_df.index, columns=ccc_df['target'].values)

        ligand_df.iloc[:, :] = expression_df[ccc_df['source'].values].values
        receptor_df.iloc[:, :] = expression_df[ccc_df['target'].values].values

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

        Reward_df = pd.DataFrame(index=ligand_df_filtered.index, columns=receptor_df_filtered.index)

        for i, ligand_gene in enumerate(ligand_df_filtered.index):
            #print('i=', i)
            for j, receptor_gene in enumerate(receptor_df_filtered.index):
                ligand_values = ligand_df_filtered.loc[ligand_gene].values.astype(float)
                receptor_values = receptor_df_filtered.loc[receptor_gene].values.astype(float)
                product_values = ligand_values * receptor_values

                if np.all(product_values == 0):
                    Reward = 0
                else:
                    log_values = np.log(product_values + 1)
                    weighted_values = log_values * ccc_df['weight'].values
                    Reward = weighted_values.sum() / ccc_num if weighted_values.sum() != 0 else 0

                Reward_df.loc[ligand_gene, receptor_gene] = Reward

        print("Reward_df:\n", Reward_df)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        output_filename = f'{output_dir}{out_cell_type}_{in_cell_type}_Reward.csv'
        Reward_df.to_csv(output_filename)
        Reward_df.to_csv(output_filename.replace('.csv', '.txt'), sep='\t')

        # long_df = Reward_df.stack().reset_index()
        # long_df.columns = ['ligand', 'receptor', 'Reward']
        # sorted_df = long_df.sort_values(by='Reward', ascending=False)
        #
        # print("sorted_df:\n", sorted_df)
        #
        # sorted_df.to_csv(output_filename.replace('_Reward.csv', '_Reward_sort.csv'), index=False)
        # sorted_df.to_csv(output_filename.replace('_Reward.csv', '_Reward_sort.txt'), sep='\t', index=False)

        print(f"Results saved as {output_filename}.")

        del ccc_list, ccc_df, ligand_df, receptor_df, ligand_df_filtered, receptor_df_filtered, Reward_df
        gc.collect()

    print("End time: " + time.asctime())

compute_Reward(
    expression_file='./data/example_data.txt',
    cell_label_file='./data/example_label.txt',
    comm_weight_file='./data/CCC_Network_binaryzation.csv',
    ccc_pair_file='./data/Cell_Type_Comm_Weight_sort.csv',
    output_dir='./Results_Reward/',
    threshold=0.3,
    top_n_pairs=5
)

