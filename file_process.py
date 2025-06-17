import pandas as pd

#---------------------------------------
# 读取文件为DataFrame
df = pd.read_csv('./data/filtered_NL.txt', sep='\t',index_col=0)
print('df:\n',df)

# 提取前20行和前100列
df_subset = df.iloc[:20, :100]
print('df_subset:\n',df_subset)

# 保存为新的文件
df_subset.to_csv('./data/example_data.txt', sep='\t')
df_subset.to_csv('./data/example_data.csv')
#---------------------------------------
# 读取文件为DataFrame
df2 = pd.read_csv('./data/filtered_NL_label.txt', sep='\t')
print('df2:\n',df2)

# 提取前100行
df2_subset = df2.iloc[:100]
print('df2_subset:\n',df2_subset)

# 保存为新的文件
df2_subset.to_csv('./data/example_label.txt', sep='\t',index=False)
df2_subset.to_csv('./data/example_label.csv',index=False)


#---------------------------------------
# 读取文件为DataFrame
df3 = pd.read_csv('./data/filtered_NL_GRN_matrix.txt', sep='\t',index_col=0)
print('df3:\n',df3)
row_index_df = df_subset.index

# 在 df3 中筛选：只保留行索引和列索引都在 df 的行索引中的部分
df3_filtered = df3.loc[
    df3.index.intersection(row_index_df),       # 行索引筛选
    df3.columns.intersection(row_index_df)      # 列索引筛选
]
print('df3_filtered:\n',df3_filtered)
df3_filtered.to_csv('./data/example_GRN.txt', sep='\t')
df3_filtered.to_csv('./data/example_GRN.csv')
