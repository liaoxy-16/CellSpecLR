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

## 01 MDIC3.py: Compute Cell-Cell Communication Matrix
### Input Files:
1 Gene Expression File (.txt)<br>
Rows represent gene names. <br>
Columns represent individual cells, indexed starting from 0. (e.g., 0, 1, 2, ...) <br>

2 Cell Type Annotation File (.txt) <br>
Contains two columns: index and label. <br>
index: cell indices starting from 0. <br>
label: the corresponding cell type for each cell. <br>

3 Gene Regulatory Network File (.txt)
Both rows and columns are gene names. <br>
Values represent the regulatory weight of the column gene on the row gene. <br>

4 Threshold Parameter:<br>
A float between 0 and 1. <br>
Determines the proportion of cell pairs with the highest communication weights to be retained.<br>
Default: 0.25 (i.e., the top 25% of cell pairs with the highest absolute communication weights are considered as having communication).<br>

### Output Files
1 CCC_Network.csv<br>
A communication weight matrix.<br>
Row and column indices represent cell indices.<br>
Values indicate the communication weight from one cell to another.<br>

2 CCC_Network_binaryzation.csv<br>
A binarized version of the communication matrix.<br>
Cell pairs with absolute weights in the top percentage (defined by threshold) are marked as 1 (indicating communication), others as 0.<br>

## 02 Cell_Type_Comm_Weight.py: Compute Communication Weights Between Cell Types
### Input Files
1 CCC_Network_binaryzation.csv<br>
A binary cell-to-cell communication matrix.<br>
Rows and columns represent cell indices.<br>
Value = 1 indicates communication exists between the two cells; 0 indicates no communication.<br>

2 Cell Type File (.txt)<br>
A two-column table with:<br>
index: cell index (starting from 0).<br>
label: corresponding cell type.<br>

### Output Files
1 Cell_Type_count.csv<br>
A table listing each cell type along with the number of cells.<br>

2 Cell_Type_Comm_Weight.csv<br>
A matrix with cell types as both rows and columns.<br>
Each value represents the total communication weight from one cell type to another.<br>

3 Cell_Type_Comm_Weight_sort.csv<br>
A ranked list of cell-type pairs based on their communication weights.<br>
Contains three columns: source cell type, target cell type, and communication weight.<br>



## 03 Reward.py: Compute Reward Score Between Genes
### Input Files
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

### Output File
Reward Matrix:<br>
A square matrix where rows and columns are gene names.<br>
Each value represents the Reward score between a pair of genes.<br>

## 04 Penalty.py: Compute Penalty Score Between Genes
### Input Files
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

### Output File
Penalty Matrix:<br>
Rows and columns are gene names.<br>
Values represent the Penalty score between each gene pair.<br>

## 05 Dependency.py: Compute Dependency Score Between Genes
### Input Files
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

### Output File
Dependency Matrix:<br>
A gene-gene matrix with Dependency scores.<br>
Row and column names are gene identifiers.<br>


## 06 Generate_protein_feature_vectors.py: Generate Protein Feature Vectors
### Input Files
1 SL_aaindex_feature.csv<br>
Rows: amino acid symbols.<br>
Columns: physicochemical properties associated with each amino acid.<br>

2 Protein Sequence File<br>
Must contain a column named Sequence that holds amino acid sequences.<br>

### Output File
1 Protein Feature Vector File<br>
The original file with the Sequence column removed.<br>
Replaced by 200 new columns: feature1 to feature200, representing numerical features for each protein.<br>

## 07 CNN-MAPLM_train.py: Train the CNN-MAPLM Model
### Input File
1 Protein Feature File<br>
Must include:<br>
feature1 ~ feature200: a 200-dimensional feature vector for each protein.<br>
label: a binary indicator
1: membrane-associated protein (potential receptor gene).<br>
0: non-membrane-associated protein.<br>

### Output Files
1 Cross-Validation Models<br>
Best-performing model for each fold in 10-fold cross-validation.<br>

2 results.csv<br>
Performance metrics recorded for each test fold (e.g., accuracy, precision, recall).<br>

## 08 CNN-MAPLM_test.py: Test the CNN-MAPLM Model
### Input File
1 Protein Feature File<br>
Must contain:<br>
feature1 ~ feature200: numerical descriptors of proteins.<br>

### Output File
1 Predicted Labels<br>
The model outputs a label (1 or 0) indicating whether each protein is predicted to be membrane-associated.<br>

## 09 EAS.py: Calculate Expression Association Strength (EAS)
### Input Folders
1 ./Results_Reward/<br>
Contains Reward matrices for each cell-type communication pair.<br>
File naming: {out_cell_type}_{in_cell_type}_Reward.txt.<br>

2 ./Results_Penalty/<br>
Contains Penalty matrices for each cell-type communication pair.<br>
File naming: {out_cell_type}_{in_cell_type}_Penalty.txt.<br>

3 ./Results_Dependency/<br>
Contains Dependency matrices for each cell-type communication pair.<br>
File naming: {out_cell_type}_{in_cell_type}_Dependency.txt.<br>

### Output Folder
1 ./Results/<br>
Stores sorted gene-gene EAS scores for each cell-type communication pair.<br>
File naming: {out_cell_type}_{in_cell_type}_EAS_sort.txt.<br>

## 10 EAS_CNN-MAPLM.py: Filter Potential Receptor Genes Using CNN-MAPLM Predictions
### Input Folder
./Results/<br>
Contains EAS-sorted gene-gene files for each communicating cell-type pair.<br>

### Output Folder
./Results_filter/<br>
Stores filtered EAS files for each cell-type pair, containing only gene pairs associated with CNN-MAPLM-predicted receptor genes.<br>






