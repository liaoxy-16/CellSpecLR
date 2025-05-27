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

## MDIC3.py: Compute Cell-to-Cell Communication Matrix
### Input Files:
1 Gene Expression File (.txt)
Rows represent gene names
Columns represent individual cells, indexed starting from 0 (e.g., 0, 1, 2, ...)

2 Cell Type Annotation File (.txt)
Contains two columns: index and label
index: cell indices starting from 0
label: the corresponding cell type for each cell

3 Gene Regulatory Network File (.txt)
Both rows and columns are gene names
Values represent the regulatory weight of the column gene on the row gene

4 Threshold Parameter:
A float between 0 and 1
Determines the proportion of cell pairs with the highest communication weights to be retained
Default: 0.25 (i.e., the top 25% of cell pairs with the highest absolute communication weights are considered as having communication)

