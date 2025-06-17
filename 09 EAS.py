import os
import pandas as pd

def compute_EAS_scores(
    reward_folder: str = "./Results_Reward",
    penalty_folder: str = "./Results_Penalty",
    dependency_folder: str = "./Results_Dependency",
    output_folder: str = "./Results"
):
    os.makedirs(output_folder, exist_ok=True)

    for reward_file in os.listdir(reward_folder):
        if reward_file.endswith("_Reward.txt"):
            cell_type_pair = reward_file.replace("_Reward.txt", "")

            penalty_file = f"{cell_type_pair}_Penalty.txt"
            dependency_file = f"{cell_type_pair}_Dependency.txt"

            reward_file_path = os.path.join(reward_folder, reward_file)
            penalty_file_path = os.path.join(penalty_folder, penalty_file)
            dependency_file_path = os.path.join(dependency_folder, dependency_file)

            print(f"{cell_type_pair}:")

            if not os.path.exists(penalty_file_path):
                print(f"Penalty file is missing：{penalty_file}")
                continue
            if not os.path.exists(dependency_file_path):
                print(f"Dependency file is missing：{dependency_file}")
                continue

            df_reward = pd.read_csv(reward_file_path, sep="\t", index_col=0)
            df_penalty = pd.read_csv(penalty_file_path, sep="\t", index_col=0)
            df_dependency = pd.read_csv(dependency_file_path, sep="\t", index_col=0)
            # print(f"df_Reward shape: {df_reward.shape}")
            # print(f"df_Penalty shape: {df_penalty.shape}")
            # print(f"df_Dependency shape: {df_dependency.shape}")

            common_genes = df_reward.index.intersection(df_penalty.index).intersection(df_dependency.index)
            common_columns = df_reward.columns.intersection(df_penalty.columns).intersection(df_dependency.columns)
            df_reward = df_reward.loc[common_genes, common_columns]
            df_penalty = df_penalty.loc[common_genes, common_columns]
            df_dependency = df_dependency.loc[common_genes, common_columns]

            # print(f"df_Reward shape: {df_reward.shape}")
            # print(f"df_Penalty shape: {df_penalty.shape}")
            # print(f"df_Dependency shape: {df_dependency.shape}")

            max1 = df_reward.max().max()
            max2 = df_penalty.max().max()
            #print('Reward.max=', max1)
            #print('Penalty.max=', max2)
            scale_factor = max1 / max2
            df_penalty *= scale_factor

            df_reward_penalty = df_reward - df_penalty

            max3 = df_reward_penalty.max().max()
            #print('Reward_Penalty.max=', max3)

            max4 = df_dependency.max().max()
            #print('Dependency.max=', max4)

            scale_factor_dependency = max3 / max4
            df_dependency *= scale_factor_dependency

            df_EAS = df_reward_penalty + df_dependency

            long_df = df_EAS.stack().reset_index()
            long_df.columns = ['ligand', 'receptor', 'EAS']
            sorted_df = long_df.sort_values(by='EAS', ascending=False)

            final_csv = os.path.join(output_folder, f"{cell_type_pair}_EAS_sort.csv")
            final_txt = os.path.join(output_folder, f"{cell_type_pair}_EAS_sort.txt")
            sorted_df.to_csv(final_csv, index=False)
            sorted_df.to_csv(final_txt, index=False, sep='\t')

            print(f"finish：{cell_type_pair}")

    print(f"\nAll files have been processed and the results are saved in {output_folder}")


compute_EAS_scores(
    reward_folder="./Results_Reward",
    penalty_folder="./Results_Penalty",
    dependency_folder="./Results_Dependency",
    output_folder="./Results"
)
