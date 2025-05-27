# CellSpecLR
![image](https://github.com/liaoxy-16/CellSpecLR/blob/master/fig1.png)

# Abstract
Ligand-receptor (L-R) interactions play a critical role in cell-cell communication, signal transduction and immune response, making their accurate prediction essential for understanding biological processes and 
advancing drug development. In this study, we propose a novel algorithm, the Cell-Specific Ligand-Receptor Pair Prediction Algorithm (CellSpecLR), designed to predict cell specific L-R pairs from single-cell data 
in tissue microenvironment. Unlike existing methods that rely on established L-R databases, our algorithm operates independently of pre-existing data, aiming to uncover a broader range of multispecies L-R interactions.
This innovative approach can enhance our understanding of biological processes and disease mechanisms, offering significant potential for new drug development and personalized medicine. To further support and extend the 
utility of our approach, we have also developed a user-friendly, multispecies L-R database website (http://compbiol.top:2023/CSMSLRdb).

# Requirements
python==3.7  
numpy==1.21.6  
pandas==1.3.5  
torch==1.10.1  
scipy==1.7.3  

# Code Usage Instructions

## 01 MDIC3.py: Compute Cell-Cell Communication Matrix
### Input Files:
1 Gene Expression File (.txt)
Rows represent gene names <br>
Columns represent individual cells, indexed starting from 0 (e.g., 0, 1, 2, ...) <br>

2 Cell Type Annotation File (.txt) <br>
Contains two columns: index and label <br>
index: cell indices starting from 0 <br>
label: the corresponding cell type for each cell <br>

3 Gene Regulatory Network File (.txt)
Both rows and columns are gene names
Values represent the regulatory weight of the column gene on the row gene

4 Threshold Parameter:
A float between 0 and 1
Determines the proportion of cell pairs with the highest communication weights to be retained
Default: 0.25 (i.e., the top 25% of cell pairs with the highest absolute communication weights are considered as having communication)

### Output Files
1 CCC_Network.csv
A communication weight matrix
Row and column indices represent cell indices
Values indicate the communication weight from one cell to another

2 CCC_Network_binaryzation.csv
A binarized version of the communication matrix
Cell pairs with absolute weights in the top percentage (defined by threshold) are marked as 1 (indicating communication), others as 0

## 02 Cell_Type_Comm_Weight.py: Compute Communication Weights Between Cell Types
### Input Files
1 CCC_Network_binaryzation.csv
A binary cell-to-cell communication matrix
Rows and columns represent cell indices
Value = 1 indicates communication exists between the two cells; 0 indicates no communication
2 Cell Type File (.txt)
A two-column table with:
index: cell index (starting from 0)
label: corresponding cell type

### Output Files
1 Cell_Type_count.csv
A table listing each cell type along with the number of cells
2 Cell_Type_Comm_Weight.csv
A matrix with cell types as both rows and columns
Each value represents the total communication weight from one cell type to another
3 Cell_Type_Comm_Weight_sort.csv
A ranked list of cell-type pairs based on their communication weights
Contains three columns: source cell type, target cell type, and communication weight



## 03 Reward.py: Compute Reward Score Between Genes
### Input Files
1 Gene Expression File (.txt)
Rows represent genes
Columns represent cells indexed from 0
2 Cell Type File (.txt)
Two columns: index (cell index) and label (cell type)
3 CCC_Network_binaryzation.csv
A binary matrix indicating cell-cell communication
4 Cell_Type_Comm_Weight_sort.csv
Ranked communication weights between cell types
Includes three columns: source type, target type, weight
5 Parameters:
Threshold (default = 0.3):
Retain genes with non-zero expression in more than 30% of the target cells involved in communication
top_n_pairs (default = 5):
Only the top N cell-type communication pairs are considered for computing gene associations

### Output File
Reward Matrix
A square matrix where rows and columns are gene names
Each value represents the Reward score between a pair of genes

## 04 Penalty.py: Compute Penalty Score Between Genes
### Input Files
1 Gene Expression File (.txt)
Rows = genes, columns = cell indices (starting from 0)
2 Cell Type File (.txt)
Includes index (cell ID) and label (cell type)
3 CCC_Network_binaryzation.csv
A binary matrix indicating whether a cell pair has significant communication (1 or 0)
4 Cell_Type_Comm_Weight_sort.csv
Ranked communication weights between cell types
5 Parameters:
Threshold (default = 0.3):
Retain genes expressed in more than 30% of the target cells involved in communication
top_n_pairs (default = 5):
Use only the top N communicating cell-type pairs to evaluate gene relationships

### Output File
Penalty Matrix
Rows and columns are gene names
Values represent the Penalty score between each gene pair

## 05 Dependency.py: Compute Dependency Score Between Genes
### Input Files
1 Gene Expression File (.txt)
Rows = genes, columns = cell indices (starting from 0)

2 Cell Type File (.txt)
Two columns: index (cell ID) and label (cell type)

3 CCC_Network_binaryzation.csv
Binary communication matrix (1 for significant communication, 0 otherwise)

4 Cell_Type_Comm_Weight_sort.csv
Ranked communication weights between cell types

5 Parameter:
Threshold (default = 0.3):
Retain genes expressed in more than 30% of the target cells involved in communication

### Output File
Dependency Matrix
A gene-gene matrix with Dependency scores
Row and column names are gene identifiers

