import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
import random
from sklearn.model_selection import KFold
from torch.utils.data import DataLoader, TensorDataset
import os

def seed_everything(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # If using GPU
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

seed_everything(42)  # Set global seed

# Load data
file_path = "./data/sequence/protein_feature_vectors.txt"
data = pd.read_csv(file_path, sep='\t')

# Extract features and labels
X = torch.tensor(data.loc[:, 'feature1':'feature200'].values, dtype=torch.float32)  # 200-dimensional features
y = torch.tensor(data["label"].values, dtype=torch.float32).unsqueeze(1)  # (batch_size, 1)

dataset = TensorDataset(X, y)

# Define 1D CNN model
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

# 10-fold cross-validation
kf = KFold(n_splits=10, shuffle=True, random_state=42)

# Create model directory
os.makedirs("./model", exist_ok=True)

# Record results
results = []

for fold, (train_idx, test_idx) in enumerate(kf.split(X)):
    print(f"\n========= Fold {fold + 1} =========")

    train_dataset = TensorDataset(X[train_idx], y[train_idx])
    test_dataset = TensorDataset(X[test_idx], y[test_idx])

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True, num_workers=0)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False, num_workers=0)

    # Reset random seed for reproducibility
    seed_everything(42)

    model = ProteinCNN()
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-4)

    # Early stopping parameters
    num_epochs = 500
    patience = 50  # Max number of epochs with no improvement
    best_val_loss = float('inf')
    best_train_loss = float('inf')
    best_epoch = 0  # Record best epoch
    early_stop_counter = 0
    best_model_state = None  # Save best model parameters

    for epoch in range(1, num_epochs + 1):
        model.train()
        total_train_loss = 0
        for batch_X, batch_y in train_loader:
            optimizer.zero_grad()
            output = model(batch_X)
            loss = criterion(output, batch_y)
            loss.backward()
            optimizer.step()
            total_train_loss += loss.item()

        avg_train_loss = total_train_loss / len(train_loader)

        # Evaluate model
        model.eval()
        total_val_loss = 0
        with torch.no_grad():
            for batch_X, batch_y in test_loader:
                output = model(batch_X)
                loss = criterion(output, batch_y)
                total_val_loss += loss.item()

        avg_val_loss = total_val_loss / len(test_loader)
        print(f"Epoch [{epoch}/{num_epochs}] - Training Loss: {avg_train_loss:.4f} - Validation Loss: {avg_val_loss:.4f}")

        # Early stopping logic
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            best_train_loss = avg_train_loss
            best_epoch = epoch  # Record best epoch
            early_stop_counter = 0  # Reset counter
            best_model_state = model.state_dict()  # Save best model
        else:
            early_stop_counter += 1

        if early_stop_counter >= patience:
            print(f"Early stopping at epoch {epoch}, best validation loss: {best_val_loss:.4f}")
            break

    print(f"Fold {fold + 1} - Best Training Loss: {best_train_loss:.4f} - Best Validation Loss: {best_val_loss:.4f} - Best Epoch: {best_epoch}")

    # Save best model for this fold
    model_path = f"./model/best_protein_cnn_fold{fold+1}.pth"
    torch.save(best_model_state, model_path)
    print(f"Model for Fold {fold + 1} saved to {model_path}")

    # Load best model for test prediction
    best_model = ProteinCNN()
    best_model.load_state_dict(torch.load(model_path))
    best_model.eval()

    y_true, y_pred = [], []
    with torch.no_grad():
        for batch_X, batch_y in test_loader:
            output = best_model(batch_X)
            pred_labels = (output >= 0.5).float()  # Binarize predictions
            y_true.extend(batch_y.numpy().flatten())
            y_pred.extend(pred_labels.numpy().flatten())

    y_true = np.array(y_true)
    y_pred = np.array(y_pred)

    # Calculate TP, TN, FP, FN
    TP = np.sum((y_pred == 1) & (y_true == 1))  # True Positive
    FP = np.sum((y_pred == 1) & (y_true == 0))  # False Positive
    FN = np.sum((y_pred == 0) & (y_true == 1))  # False Negative
    TN = np.sum((y_pred == 0) & (y_true == 0))  # True Negative

    # Calculate Precision, Recall, F1-score
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    # Record detailed results for each fold
    results.append({
        "Fold": fold + 1,
        "Best Epoch": best_epoch,
        "Train Loss": best_train_loss,
        "Val Loss": best_val_loss,
        "Precision": precision,
        "Recall": recall,
        "F1 Score": f1_score,
        "TP": TP,
        "TN": TN,
        "FP": FP,
        "FN": FN
    })

# Generate results DataFrame and save
results_df = pd.DataFrame(results)
results_df.to_csv("./model/results.csv", index=False)
print("\nResults saved to ./model/results.csv")
