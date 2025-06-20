# CellSpecLR
![image](https://github.com/liaoxy-16/CellSpecLR/blob/master/fig1.png)

# Abstract
Ligand-receptor (L-R) interactions play a critical role in cell-cell communication, signal transduction and immune response, making their accurate prediction essential for understanding biological processes and advancing drug development. In this study, we propose a novel algorithm, the Cell-Specific Ligand-Receptor Pair Prediction Algorithm (CellSpecLR), designed to predict cell specific L-R pairs from single-cell data in tissue microenvironment. Unlike existing methods that rely on established L-R databases, our algorithm operates independently of pre-existing data, aiming to uncover a broader range of multispecies L-R interactions. This innovative approach can enhance our understanding of biological processes and disease mechanisms, offering significant potential for new drug development and personalized medicine. To further support and extend the utility of our approach, we have also developed a user-friendly, multispecies L-R database website (http://compbiol.top:2023/CSMSLRdb).

# Requirements
python>=3.7  
numpy>=1.21.6  
pandas>=1.3.5  
torch>=1.10.1  
scipy>=1.7.3  

# Code Usage Instructions

## Code Usage Example
Below is the standard workflow and description of each script in this project:
### 1. Calculate the Cell-Cell Communication Matrix
Run the script python MDIC3_01.py to compute the cell-cell communication matrix.

#### Input files:
Gene expression matrix: ./data/example_data.txt<br>
Gene regulatory network matrix: ./data/example_GRN.txt<br>

#### Output files:
Communication matrix: ./data/CCC_Network.csv<br>
Binarized communication matrix: ./data/CCC_Network_binaryzation.csv<br>

### 2. Calculate Communication Weights Between Cell Types
Run the script python Cell_Type_Comm_Weight_02.py to compute communication weights among different cell types.

#### Input files:
Binarized communication matrix: ./data/CCC_Network_binaryzation.csv<br>
Cell type annotation file: ./data/example_label.txt<br>

#### Output files:
Cell type counts: Cell_Type_count.csv<br>
Cell-type communication weight matrix: Cell_Type_Comm_Weight.csv<br>
Sorted communication weights: Cell_Type_Comm_Weight_sort.csv<br>


### 3. Compute the Reward Metric
Run the script python Reward_03.py to calculate the Reward metric for each communicating cell-type pair.

#### Input files:
Gene expression matrix: ./data/example_data.txt<br>
Cell type annotation: ./data/example_label.txt<br>
Binarized communication matrix: ./data/CCC_Network_binaryzation.csv<br>
Sorted communication weight file: Cell_Type_Comm_Weight_sort.csv<br>

#### Output directory:
All results are saved in ./Results_Reward/<br>
File naming format: {out_cell_type}_{in_cell_type}_Reward.txt<br>

### 4. Compute the Penalty Metric
Run the script python Penalty_04.py to calculate the Penalty metric.

#### Input files: Same as for the Reward step

#### Output directory: ./Results_Penalty/
File naming format: {out_cell_type}_{in_cell_type}_Penalty.txt

### 5. Compute the Dependency Metric
Run the script python Dependency_05.py to calculate the Dependency metric.

#### Input files: Same as for the Reward step

#### Output directory: ./Results_Dependency/
File naming format: {out_cell_type}_{in_cell_type}_Dependency.txt

### 6. Generate Protein Feature Vectors (for Receptor Gene Selection)
Run the script python Generate_protein_feature_vectors_06.py

#### Input files:
Physicochemical properties: ./data/aaindex/SL_aaindex_feature.csv<br>
Protein sequence data: ./data/sequence/example_sequence.csv<br>

#### Output file:
Feature vectors: ./data/sequence/protein_feature_vectors.txt

### 7. Train Membrane-Associated Protein Prediction Model
Run the script python CNN_MAPLM_train_07.py

#### Input file:
Protein feature vectors: ./data/sequence/protein_feature_vectors.txt

#### Output:
Optimal model obtained via 10-fold cross-validation

### 8. Test Membrane-Associated Protein Prediction Model
Run the script python CNN_MAPLM_test_08.py to predict membrane association for unlabeled proteins.

#### Input file:
Unlabeled protein feature vectors: ./data/sequence/protein_feature_vectors.txt

#### Output:
Predicted labels for each protein: values close to 1 indicate membrane-associated proteins, and values close to 0 indicate non-membrane proteins

### 9. Compute EAS Scores
Run the script python EAS_09.py to integrate Reward, Penalty, and Dependency metrics into a final EAS score.

#### Input directories:
./Results_Reward/<br>
./Results_Penalty/<br>
./Results_Dependency/<br>

#### Output directory: ./Results/
File naming format: {out_cell_type}_{in_cell_type}_EAS_sort.txt

### 10. Compute Final Filtered EAS Scores
Run the script python EAS_CNN_MAPLM_10.py

#### Input files:
Protein label file: ./data/sequence/example_protein_label.csv (label=1 indicates membrane-associated proteins)<br>
Original EAS result folder: ./Results/<br>

#### Output directory: ./Results_filter/
File naming format: {out_cell_type}_{in_cell_type}_EAS_sort.txt, representing the filtered final EAS scores after receptor selection


## Descriptions for Each Script
### MDIC3_01.py: Compute Cell-Cell Communication Matrix
#### Input Files:
1 Gene Expression File (.txt)<br>
Rows represent gene names. <br>
Columns represent individual cells, indexed starting from 0. (e.g., 0, 1, 2, ...) <br>
|       | 0   | 1   | ... |  n  |
|-------|-----|-----|-----|-----|
|gene 1 |     |     |     |     |
|gene 2 |     |     |     |     |
|  ...  |     |     |     |     |
|gene m |     |     |     |     |

2 Gene Regulatory Network File (.txt)
Both rows and columns are gene names. <br>
Values represent the regulatory weight of the column gene on the row gene. <br>
|      |gene 1|gene 2| ...  |gene m|
|------|------|------|------|------|
|gene 1|      |      |      |      |
|gene 2|      |      |      |      |
|  ... |      |      |      |      |
|gene m|      |      |      |      |

3 Threshold Parameter:<br>
A float between 0 and 1. <br>
Determines the proportion of cell pairs with the highest communication weights to be retained.<br>
Default: 0.25 (i.e., the top 25% of cell pairs with the highest absolute communication weights are considered as having communication).<br>

#### Output Files
1 CCC_Network.csv<br>
A communication weight matrix.<br>
Row and column indices represent cell indices.<br>
Values indicate the communication weight from one cell to another.<br>

2 CCC_Network_binaryzation.csv<br>
A binarized version of the communication matrix.<br>
Cell pairs with absolute weights in the top percentage (defined by threshold) are marked as 1 (indicating communication), others as 0.<br>

### Cell_Type_Comm_Weight_02.py: Compute Communication Weights Between Cell Types
#### Input Files
1 CCC_Network_binaryzation.csv<br>
A binary cell-to-cell communication matrix.<br>
Rows and columns represent cell indices.<br>
Value = 1 indicates communication exists between the two cells; 0 indicates no communication.<br>

2 Cell Type File (.txt) <br>
Contains two columns: index and label. <br>
index: cell indices starting from 0. <br>
label: the corresponding cell type for each cell. <br>
|index|label|
|-----|-----|
|  0  |     |  
|  1  |     |   
| ... |     |     
|  n  |     | 

#### Output Files
1 Cell_Type_count.csv<br>
A table listing each cell type along with the number of cells.<br>

2 Cell_Type_Comm_Weight.csv<br>
A matrix with cell types as both rows and columns.<br>
Each value represents the total communication weight from one cell type to another.<br>

3 Cell_Type_Comm_Weight_sort.csv<br>
A ranked list of cell-type pairs based on their communication weights.<br>
Contains three columns: source cell type, target cell type, and communication weight.<br>



### Reward_03.py: Compute Reward Score Between Genes
#### Input Files
1 Gene Expression File (.txt)<br>
Rows represent genes.<br>
Columns represent cells indexed from 0.<br>

2 Cell Type File (.txt)<br>
Two columns: index (cell index) and label (cell type).<br>

3 CCC_Network_binaryzation.csv<br>
A binary matrix indicating cell-cell communication.<br>

4 Cell_Type_Comm_Weight_sort.csv<br>
Ranked communication weights between cell types.<br>
Includes three columns: source type, target type, weight.<br>

5 Parameters:<br>
Threshold (default = 0.3):<br>
Retain genes with non-zero expression in more than 30% of the target cells involved in communication.<br>
top_n_pairs (default = 5):<br>
Only the top N cell-type communication pairs are considered for computing gene associations.<br>

#### Output File
Reward Matrix:<br>
A square matrix where rows and columns are gene names.<br>
Each value represents the Reward score between a pair of genes.<br>

### Penalty_04.py: Compute Penalty Score Between Genes
#### Input Files
1 Gene Expression File (.txt)<br>
Rows = genes, columns = cell indices (starting from 0).<br>

2 Cell Type File (.txt)<br>
Includes index (cell ID) and label (cell type).<br>

3 CCC_Network_binaryzation.csv<br>
A binary matrix indicating whether a cell pair has significant communication (1 or 0).<br>

4 Cell_Type_Comm_Weight_sort.csv<br>
Ranked communication weights between cell types.<br>

5 Parameters:<br>
Threshold (default = 0.3):<br>
Retain genes expressed in more than 30% of the target cells involved in communication.<br>
top_n_pairs (default = 5):<br>
Use only the top N communicating cell-type pairs to evaluate gene relationships.<br>

#### Output File
Penalty Matrix:<br>
Rows and columns are gene names.<br>
Values represent the Penalty score between each gene pair.<br>

### 05 Dependency.py: Compute Dependency Score Between Genes
#### Input Files
1 Gene Expression File (.txt)<br>
Rows = genes, columns = cell indices (starting from 0).<br>

2 Cell Type File (.txt)<br>
Two columns: index (cell ID) and label (cell type).<br>

3 CCC_Network_binaryzation.csv<br>
Binary communication matrix (1 for significant communication, 0 otherwise).<br>

4 Cell_Type_Comm_Weight_sort.csv<br>
Ranked communication weights between cell types.<br>

5 Parameter:<br>
Threshold (default = 0.3):<br>
Retain genes expressed in more than 30% of the target cells involved in communication.<br>
top_n_pairs (default = 5):<br>
Use only the top N communicating cell-type pairs to evaluate gene relationships.<br>

#### Output File
Dependency Matrix:<br>
A gene-gene matrix with Dependency scores.<br>
Row and column names are gene identifiers.<br>


### Generate_protein_feature_vectors_06.py: Generate Protein Feature Vectors
#### Input Files
1 SL_aaindex_feature.csv<br>
Rows: amino acid symbols.<br>
Columns: physicochemical properties associated with each amino acid.<br>

2 Protein Sequence File<br>
Must contain a column named Sequence that holds amino acid sequences.<br>

#### Output File
1 Protein Feature Vector File<br>
The original file with the Sequence column removed.<br>
Replaced by 200 new columns: feature1 to feature200, representing numerical features for each protein.<br>

### CNN_MAPLM_train_07.py: Train the CNN-MAPLM Model
#### Input File
1 Protein Feature File<br>
Must include:<br>
feature1 ~ feature200: a 200-dimensional feature vector for each protein.<br>
label: a binary indicator
1: membrane-associated protein (potential receptor gene).<br>
0: non-membrane-associated protein.<br>

#### Output Files
1 Cross-Validation Models<br>
Best-performing model for each fold in 10-fold cross-validation.<br>

2 results.csv<br>
Performance metrics recorded for each test fold (e.g., accuracy, precision, recall).<br>

### CNN_MAPLM_test_08.py: Test the CNN-MAPLM Model
#### Input File
1 Protein Feature File<br>
Must contain:<br>
feature1 ~ feature200: numerical descriptors of proteins.<br>

#### Output File
1 Predicted Labels<br>
The model outputs a label (1 or 0) indicating whether each protein is predicted to be membrane-associated.<br>

### EAS_09.py: Calculate Expression Association Strength (EAS)
#### Input Folders
1 ./Results_Reward/<br>
Contains Reward matrices for each cell-type communication pair.<br>
File naming: {out_cell_type}_{in_cell_type}_Reward.txt.<br>

2 ./Results_Penalty/<br>
Contains Penalty matrices for each cell-type communication pair.<br>
File naming: {out_cell_type}_{in_cell_type}_Penalty.txt.<br>

3 ./Results_Dependency/<br>
Contains Dependency matrices for each cell-type communication pair.<br>
File naming: {out_cell_type}_{in_cell_type}_Dependency.txt.<br>

#### Output Folder
1 ./Results/<br>
Stores sorted gene-gene EAS scores for each cell-type communication pair.<br>
File naming: {out_cell_type}_{in_cell_type}_EAS_sort.txt.<br>

### EAS_CNN_MAPLM_10.py: Filter Potential Receptor Genes Using CNN-MAPLM Predictions
#### Input Folder
./Results/<br>
Contains EAS-sorted gene-gene files for each communicating cell-type pair.<br>

#### Output Folder
./Results_filter/<br>
Stores filtered EAS files for each cell-type pair, containing only gene pairs associated with CNN-MAPLM-predicted receptor genes.<br>






