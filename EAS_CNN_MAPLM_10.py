import os
import pandas as pd


results_folder = './Results'
delete_folder = './Results_filter'

os.makedirs(delete_folder, exist_ok=True)

# receptor: label=1
df_labels = pd.read_csv('./data/sequence/example_protein_label.csv')
filtered_genes = df_labels[df_labels['label'] == 1]['Gene Names']
print('df_labels:\n',df_labels)
print('filtered_genes:\n',filtered_genes)
genes_to_keep = set(filtered_genes)

print(f'sum(label=1)：{len(genes_to_keep)}')

for file_name in os.listdir(results_folder):
    if file_name.endswith('_EAS_sort.txt'):
        input_file_path = os.path.join(results_folder, file_name)
        print(f'{input_file_path}')

        df = pd.read_csv(input_file_path, sep='\t')

        filtered_df = df[df['receptor'].isin(genes_to_keep)]
        print('filtered_df:\n', filtered_df)

        output_file_path = os.path.join(delete_folder, file_name)

        filtered_df.to_csv(output_file_path, sep='\t', index=False)

        print(f'finish: {file_name} -> {output_file_path}')

print(f'All files have been processed and the results are saved in {delete_folder}')
