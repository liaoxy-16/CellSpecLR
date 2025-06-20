import pandas as pd
import time
import numpy as np
import scipy.stats as stats
from scipy.stats import norm
import os

def compute_Dependency(
    expression_file: str,
    cell_label_file: str,
    comm_weight_file: str,
    ccc_pair_file: str,
    threshold: float = 0.3,
    top_n_pairs: int = 5
):
    begin = time.asctime()
    print("Begin time: " + begin)

    expression_df = pd.read_csv(expression_file, index_col=0, sep='\t')  # Gene expression matrix
    cell_type_df = pd.read_csv(cell_label_file, sep='\t') # Cell type information, including 'index' and 'label'
    comm_weight_df = pd.read_csv(comm_weight_file, index_col=0).abs()  # Intercellular communication weight

    expression_df.columns = expression_df.columns.astype(int)
    cell_type_df['index'] = cell_type_df['index'].astype(int)
    comm_weight_df.index = comm_weight_df.index.astype(int)
    comm_weight_df.columns = comm_weight_df.columns.astype(int)

    print('expression_df:\n', expression_df)
    print('cell_type_df:\n', cell_type_df)
    print('comm_weight_df:\n', comm_weight_df)

    CellType_CCC_list_df = pd.read_csv(ccc_pair_file)
    #CellType_CCC_list_df = CellType_CCC_list_df[CellType_CCC_list_df['Out_CellType'] != CellType_CCC_list_df['In_CellType']]
    print('CellType_CCC_list_df:\n', CellType_CCC_list_df)

    CellType_CCC_list_df = CellType_CCC_list_df[
        CellType_CCC_list_df['Out_CellType'] != CellType_CCC_list_df['In_CellType']
        ]
    print('Filtered CellType_CCC_list_df:\n', CellType_CCC_list_df)

    CellType_CCC_list_df = CellType_CCC_list_df.iloc[0:top_n_pairs, :]
    print('Sliced CellType_CCC_list_df:\n', CellType_CCC_list_df)

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
        print('ccc_df:\n', ccc_df)

        ccc_df = ccc_df[ccc_df['weight'] != 0]
        print('Filtered ccc_df:\n', ccc_df)

        ligand_df = pd.DataFrame(index=expression_df.index, columns=ccc_df['source'].values)
        receptor_df = pd.DataFrame(index=expression_df.index, columns=ccc_df['target'].values)

        ligand_values = expression_df[ccc_df['source'].values].values
        receptor_values = expression_df[ccc_df['target'].values].values

        ligand_df.iloc[:, :] = ligand_values
        receptor_df.iloc[:, :] = receptor_values

        print("Ligand DataFrame:\n", ligand_df)
        print("Receptor DataFrame:\n", receptor_df)

        threshold_val = threshold

        ligand_non_zero_count = ligand_df.apply(lambda row: (row != 0).sum(), axis=1)
        ligand_row_threshold = int(ligand_df.shape[1] * threshold_val)
        ligand_df_filtered = ligand_df[ligand_non_zero_count > ligand_row_threshold]

        receptor_non_zero_count = receptor_df.apply(lambda row: (row != 0).sum(), axis=1)
        receptor_row_threshold = int(receptor_df.shape[1] * threshold_val)
        receptor_df_filtered = receptor_df[receptor_non_zero_count > receptor_row_threshold]

        print("Filtered Ligand DataFrame:\n", ligand_df_filtered)
        print("Filtered Receptor DataFrame:\n", receptor_df_filtered)


        ligand_df_filtered.columns = range(len(ligand_df_filtered.columns))
        receptor_df_filtered.columns = range(len(receptor_df_filtered.columns))
        print("Filtered Ligand DataFrame:\n", ligand_df_filtered)
        print("Filtered Receptor DataFrame:\n", receptor_df_filtered)

        print('ligand_num=', len(ligand_df_filtered.index), ' receptor_num=', len(receptor_df_filtered.index))

        ligand_array = ligand_df_filtered.values.astype(float)
        receptor_array = receptor_df_filtered.values.astype(float)

        ligand_mean = ligand_array.mean(axis=1, keepdims=True)
        ligand_std = ligand_array.std(axis=1, keepdims=True)
        ligand_normalized = (ligand_array - ligand_mean) / ligand_std

        receptor_mean = receptor_array.mean(axis=1, keepdims=True)
        receptor_std = receptor_array.std(axis=1, keepdims=True)
        receptor_normalized = (receptor_array - receptor_mean) / receptor_std

        n_samples = ligand_array.shape[1]

        correlation_matrix = np.dot(ligand_normalized, receptor_normalized.T) / n_samples

        t_values = correlation_matrix * np.sqrt((n_samples - 2) / (1 - correlation_matrix ** 2))
        pvalue_matrix = 2 * (1 - norm.cdf(np.abs(t_values)))

        log_pvalue_matrix = -np.log10(pvalue_matrix)
        finite_vals = log_pvalue_matrix[np.isfinite(log_pvalue_matrix)]
        max_val = finite_vals.max() if finite_vals.size > 0 else 0
        log_pvalue_matrix[np.isinf(log_pvalue_matrix)] = max_val


        logp_df = pd.DataFrame(
            log_pvalue_matrix,
            index=ligand_df_filtered.index,
            columns=receptor_df_filtered.index
        )


        print("logp_df:\n", logp_df)

        #os.makedirs('./Results_Dependency/', exist_ok=True)
        if not os.path.exists('./Results_Dependency/'):
            os.makedirs('./Results_Dependency/')


        logp_df.to_csv(f'./Results_Dependency/{out_cell_type}_{in_cell_type}_Dependency.csv')
        logp_df.to_csv(f'./Results_Dependency/{out_cell_type}_{in_cell_type}_Dependency.txt', sep='\t')


    print("End time: " + time.asctime())

compute_Dependency(
    expression_file='./data/example_data.txt',
    cell_label_file='./data/example_label.txt',
    comm_weight_file='./data/CCC_Network_binaryzation.csv',
    ccc_pair_file='./data/Cell_Type_Comm_Weight_sort.csv',
    threshold=0.3,
    top_n_pairs=5
)
