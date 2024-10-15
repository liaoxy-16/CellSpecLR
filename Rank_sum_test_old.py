import pandas as pd
from scipy.stats import ranksums

expression_data = pd.read_csv('./data/data.txt', delimiter='\t',index_col=0)
print('expression_data:\n',expression_data)
expression_data.columns = range(1,len(expression_data.columns)+1)

cell_types = pd.read_csv('./data/data_label.txt',delimiter='\t')
print('cell_types:\n',cell_types)


# 将细胞类型信息合并到基因表达矩阵中
merged_data = expression_data .T
print('merged_data:\n',merged_data)
merged_data['label'] = cell_types.set_index('index')['label']
print('merged_data:\n',merged_data)


# 进行Wilcoxon秩和检验
significantly_different_genes = []

for gene in expression_data.index:
    expression_data_temp = merged_data[gene]
    for cell_type in cell_types['label'].unique():
        # 将当前细胞类型与其他类型的细胞进行比较
        current_type_expression = expression_data_temp[merged_data['label'] == cell_type]
        other_types_expression = expression_data_temp[merged_data['label'] != cell_type]

        # 执行Wilcoxon秩和检验
        stat, p_value = ranksums(current_type_expression, other_types_expression)

        if p_value < 0.05:  # 显著性水平0.05
            significantly_different_genes.append((gene, cell_type, p_value))

# 输出差异表达的基因及其在哪个细胞类型中普遍表达
# for gene, cell_type, p_value in significantly_different_genes:
#     print(f"基因 {gene} 在细胞类型 {cell_type} 中差异表达 (p值={p_value})")

# 将diff_genes列表转换为DataFrame
diff_genes_df = pd.DataFrame(significantly_different_genes, columns=['gene', 'label', 'p_value'])
print('diff_genes_df:\n',diff_genes_df)
# 保存DataFrame为CSV文件
diff_genes_df.to_csv('./data/ranksums_gene_label.csv', index=False)

# 将 DataFrame 保存为 txt 文件
diff_genes_df.to_csv('./data/ranksums_gene_label.txt', sep='\t', index=False)

# 读取CSV文件的第一列数据
data = pd.read_csv('./data/ranksums_gene_label.txt', delimiter='\t')
print('data:\n',data)
column_data = data.iloc[:, 0]
print('column_data:\n',column_data)

# 删除重复的元素
unique_data = column_data.drop_duplicates()
print('unique_data:\n',unique_data)

# 将结果输出到新的CSV文件
#unique_data.to_csv('./data/filter_gene_list.txt', index=False, header=['gene'])
unique_data.to_csv('./data/ranksums_gene_list.txt', index=False, sep='\t')

#---------------------------------------------------
# 将通过秩和检验删选出的基因重新生成基因表达矩阵
# 读取第一个CSV文件，获取基因名列表
gene_list = pd.read_csv('./data/ranksums_gene_list.txt', delimiter='\t')
#gene_list = gene_list.iloc[:, 0]
print('gene_list:\n',gene_list)
gene_list = gene_list.iloc[:, 0]
print('gene_list:\n',gene_list)

# 读取第二个CSV文件，基因作为行索引
expression_matrix = pd.read_csv('./data/data.txt', delimiter='\t',index_col=0)
print('expression_matrix:\n',expression_matrix)

expression_matrix.columns = range(1,len(expression_matrix.columns)+1)
print('expression_matrix:\n',expression_matrix)

# 根据基因名进行行筛选
ranksums_matrix = expression_matrix.loc[gene_list, :]
print('ranksums_matrix:\n',ranksums_matrix)

# 将结果写入新的CSV文件
#filtered_matrix.to_csv('./data/filter_fish.txt', index=False, sep='\t')
ranksums_matrix.to_csv('./data/ranksums_data.txt',sep='\t')