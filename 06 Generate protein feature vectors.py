import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

def extract_features(protein_sequence):
    """ Compute AAC + SOC feature vector for a protein sequence """
    # Compute AAC features (20 dimensions)
    aac_vector = np.array([protein_sequence.count(aa) / len(protein_sequence) for aa in amino_acids])

    # Compute SOC features (180 dimensions)
    soc_vector = []
    for prop in df_normalized.columns.tolist():  # Traverse 9 physicochemical properties
        prop_values = df_normalized.loc[amino_acids, prop].to_dict()  # Get amino acid property values
        seq_values = np.array([prop_values[aa] for aa in protein_sequence if aa in prop_values])

        # Compute sequence-order correlation
        for lambda_val in lambda_values:
            if len(seq_values) > lambda_val:
                soc_value = np.mean(seq_values[:-lambda_val] * seq_values[lambda_val:])
            else:
                soc_value = 0  # Set to 0 if sequence is too short
            soc_vector.append(soc_value)

    # Combine AAC + SOC to form a 200-dimensional feature
    feature_vector = np.concatenate([aac_vector, soc_vector])
    return feature_vector

# Read physicochemical property feature data
feature_file = "./data/aaindex/SL_aaindex_feature.csv"  # Physicochemical properties file
protein_file = "./data/sequence/example_sequence.csv"  # Protein sequence file

# Read physicochemical property data, set amino acid abbreviations as row index
df = pd.read_csv(feature_file, index_col=0)
print('df:\n',df)

# Normalize 9 physicochemical properties to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
df_normalized = pd.DataFrame(scaler.fit_transform(df), index=df.index, columns=df.columns)
print('df_normalized:\n', df_normalized)

# Define 20 standard amino acids
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
print('amino_acids=',amino_acids)

# Set lambda range
lambda_values = list(range(1, 21))

# Read protein sequence data
data = pd.read_csv(protein_file)
print('data:\n',data)

# Remove rows with empty Sequence and reset index
data = data.dropna(subset=["Sequence"]).reset_index(drop=True)
print('filtered data:\n', data)

# Batch compute feature vectors for all protein sequences
feature_matrix = data["Sequence"].apply(extract_features)
print('feature_matrix.shape=',feature_matrix.shape)
print('feature_matrix:\n',feature_matrix)

# Create feature DataFrame
feature_columns = [f"feature{i + 1}" for i in range(200)]
feature_df = pd.DataFrame(feature_matrix.tolist(), columns=feature_columns)

# Combine all columns in data except Sequence with the feature vectors
metadata_columns = [col for col in data.columns if col != "Sequence"]
output_df = pd.concat([data[metadata_columns], feature_df], axis=1)
print('output_df:\n', output_df)

# Save to CSV file
output_df.to_csv("./data/sequence/protein_feature_vectors.csv", index=False)
output_df.to_csv("./data/sequence/protein_feature_vectors.txt", index=False, sep='\t')
