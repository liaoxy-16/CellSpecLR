import torch
import torch.nn as nn
import pandas as pd
import numpy as np

# Load data
file_path = "./data/sequence/protein_feature_vectors.txt"
data = pd.read_csv(file_path, sep='\t')
print('data:\n',data)

# Extract 200-dimensional feature columns
feature_cols = [f"feature{i}" for i in range(1, 201)]
X = torch.tensor(data[feature_cols].values, dtype=torch.float32)

# Define 1D CNN structure (must match the training structure)
class ProteinCNN(nn.Module):
    def __init__(self):
        super(ProteinCNN, self).__init__()
        self.conv1 = nn.Conv1d(in_channels=1, out_channels=32, kernel_size=3, padding=1)
        self.conv2 = nn.Conv1d(in_channels=32, out_channels=64, kernel_size=3, padding=1)
        self.pool = nn.MaxPool1d(kernel_size=2)
        self.fc1 = nn.Linear(64 * 50, 128)
        self.fc2 = nn.Linear(128, 1)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = x.unsqueeze(1)  # (batch_size, 1, 200)
        x = self.pool(self.relu(self.conv1(x)))  # (batch_size, 32, 100)
        x = self.pool(self.relu(self.conv2(x)))  # (batch_size, 64, 50)
        x = x.view(x.size(0), -1)  # Flatten
        x = self.relu(self.fc1(x))
        x = self.sigmoid(self.fc2(x))
        return x

# Load 10 models
num_models = 10
model_paths = [f"./model/best_protein_cnn_fold{i+1}.pth" for i in range(num_models)]

# Store predictions from 10 models
all_predictions = np.zeros((len(X), num_models))

# Load each model and make predictions
for i, model_path in enumerate(model_paths):
    model = ProteinCNN()
    model.load_state_dict(torch.load(model_path))
    model.eval()  # Set to evaluation mode

    # Predict
    with torch.no_grad():
        all_predictions[:, i] = model(X).squeeze().numpy()  # Store predictions from each model

# Compute average prediction values across 10 models
avg_prediction = np.mean(all_predictions, axis=1)

# Generate binary classification result (ensemble prediction)
binary_pred = (avg_prediction >= 0.5).astype(int)

# Create DataFrame to save average predictions and final prediction results
non_feature_cols = [col for col in data.columns if col not in feature_cols]
output_df = data[non_feature_cols].copy()
output_df["Prediction"] = avg_prediction
output_df["Binary Prediction"] = binary_pred

# Save DataFrame
output_df.to_csv("./data/sequence/protein_prediction_label.csv", index=False)
output_df.to_csv("./data/sequence/protein_prediction_label.txt", index=False, sep='\t')
