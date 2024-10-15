import pandas as pd
import numpy as np

ccc_label_list = pd.read_csv('./data/CCC_list.txt', delimiter='\t')
print('ccc_label_list:\n', ccc_label_list)

ccc_label_list = ccc_label_list[(ccc_label_list['label_row'] == 'T cell') & (ccc_label_list['label_column'] == 'T cell')]
print('ccc_label_list(T cell->T cell):\n', ccc_label_list)

abs_ccc = pd.read_csv('./data/CCC_result.csv', index_col=0)
print('ccc:\n', abs_ccc)
abs_ccc = abs_ccc.abs()
abs_ccc.index = range(1,len(abs_ccc.index)+1)
abs_ccc.columns = range(1,len(abs_ccc.columns)+1)
print('abs_ccc:\n', abs_ccc)

#print(abs_ccc.at[1,2])
ccc_label_list['weight'] = ccc_label_list.apply(lambda row: abs_ccc.loc[row['Row'], row['Column']], axis=1)
print('ccc_label_list(T cell->T cell):\n', ccc_label_list)

ccc_label_list.index = range(len(ccc_label_list.index))
print('final ccc_label_list(T cell->T cell):\n', ccc_label_list)

print('The expression file is being read...')
# 读取 CSV 文件，将第一行作为列索引，第一列作为行索引
df_fish = pd.read_csv('./data/data.txt', delimiter='\t', index_col=0)
#df_fish = pd.read_csv('./data_preprocess/filter_fish.txt', delimiter='\t', index_col=0)
df_fish.columns = range(1,len(df_fish.columns)+1)

# Removing index names
df_fish.index.name = None

print('expression matrix:\n', df_fish)
# df_fish = df_fish.iloc[:,0:100]
# print('expression matrix 100:\n', df_fish)

#-----------------------------------------------------
print('获得基因对应的出度细胞和入度细胞')
send_cell_list = ccc_label_list['Row'].tolist()
receive_cell_list = ccc_label_list['Column'].tolist()
print('send_cell_list:\n',send_cell_list)
print('receive_cell_list:\n',receive_cell_list)

# iloc 索引是从0开始
send_cell_list = [x-1 for x in send_cell_list]
receive_cell_list = [x-1 for x in receive_cell_list]
print('send_cell_list:\n',send_cell_list)
print('receive_cell_list:\n',receive_cell_list)

print('len(send_cell_list)=',len(send_cell_list))
print('len(receive_cell_list)=',len(receive_cell_list))


send_cell_data = df_fish.iloc[:, send_cell_list]
send_cell_data.columns = range(0,len(send_cell_data.columns))
print('send_cell_data:\n',send_cell_data)

receive_cell_data = df_fish.iloc[:, receive_cell_list]
receive_cell_data.columns = range(0,len(receive_cell_data.columns))
print('receive_cell_data:\n',receive_cell_data)

#----------------------------------------------------------------------
#计算每一行非零元素的数量
non_zero_counts_df1 = (send_cell_data != 0).sum(axis=1)
non_zero_counts_df2 = (receive_cell_data != 0).sum(axis=1)

# 分别筛选出两个 DataFrame 中一行的零元素数量小于 通讯细胞对数量的%50 的行索引
ccc_length = len(send_cell_list)
print('CCC num=',ccc_length)
filtered_indices_send = send_cell_data.index[non_zero_counts_df1 < (ccc_length*0.5)]
filtered_indices_receive = receive_cell_data.index[non_zero_counts_df2 < (ccc_length*0.5)]
# 删除这些行
send_cell_data = send_cell_data.drop(index=filtered_indices_send)
receive_cell_data = receive_cell_data.drop(index=filtered_indices_receive)
print('delete 0 send_cell_data:\n',send_cell_data)
print('delete 0 receive_cell_data:\n',receive_cell_data)

#send_cell_data.to_csv('./data/T cell-T cell/send_cell_data.csv')
#receive_cell_data.to_csv('./data/T cell-T cell/receive_cell_data.csv')


row = send_cell_data.index
column = receive_cell_data.index

send_cell_data.index = range(len(send_cell_data.index))
receive_cell_data.index = range(len(receive_cell_data.index))
print('final send_cell_data:\n',send_cell_data)
print('final receive_cell_data:\n',receive_cell_data)

#---------------------------------------------------------------------
# 创建一个全零矩阵
score_matrix = np.zeros((len(send_cell_data.index), len(receive_cell_data.index)))
print('shape=',score_matrix.shape)

count = 0
for i in range(score_matrix.shape[0]):
    print(count)
    count += 1
    for j in range(score_matrix.shape[1]):
        score = 0
        for k in range(len(ccc_label_list.index)):
            score = score + ccc_label_list.loc[k,'weight']*\
                    np.log(send_cell_data.loc[i,k]*receive_cell_data.loc[j,k]+1)
        score = score/len(ccc_label_list.index)
        score_matrix[i][j] = score

Absolute_Score = pd.DataFrame(score_matrix,index=row,columns=column)
print('Absolute_Score:\n',Absolute_Score)

Absolute_Score.to_csv('./data/T cell-T cell/T cell-T cell_AAS.csv')
print('-----------------------------------------------')
Relative_Score = Absolute_Score
# 计算每行和每列的平均值
#row_means = Absolute_Score.mean(axis=1)
#column_means = Absolute_Score.mean(axis=0)

# 计算每行和每列的平均值
row_means = 0.5 * Absolute_Score.quantile(0.5, axis=1) + \
            0.25 * (Absolute_Score.quantile(0.25, axis=1) + Absolute_Score.quantile(0.75, axis=1))

column_means = 0.5 * Absolute_Score.quantile(0.5) + \
               0.25 * (Absolute_Score.quantile(0.25) + Absolute_Score.quantile(0.75))

print(row_means)
print(column_means)
# 对每个元素进行操作
for i in range(Relative_Score.shape[0]):
    for j in range(Relative_Score.shape[1]):
        Relative_Score.iloc[i, j] = Relative_Score.iloc[i, j]-(row_means[i]+column_means[j])

print('Relative_Score:\n', Relative_Score)
Relative_Score.to_csv('./data/T cell-T cell/T cell-T cell_RAS.csv')
print('-----------------------------------------------')
# 从大到小排序 dropna=Fafishe：是NAN也保留
sorted_df = Relative_Score.stack(dropna=False).sort_values(ascending=False)
print('sorted_df:\n', sorted_df)
#sorted_df.to_csv('./data/score_uniprot_fish_send-receive_sort.txt', sep='\t', header=['gene1', 'gene2', 'score'])

# 重置索引，将多级索引转换为普通列
sorted_df = sorted_df.reset_index()
new_column_names = ['gene1', 'gene2', 'RAS']
sorted_df.columns = new_column_names
print('sorted_df:\n', sorted_df)

#保存前100000行
# sorted_df_subset = sorted_df.head(100000)
# print('sorted_df_subset:\n', sorted_df_subset)
sorted_df.to_csv('./data/T cell-T cell/T cell-T cell_RAS_sort.txt', index=False, sep='\t')
#sorted_df_subset.to_csv('./data/ALL_fish_APOE+FIB-FBN1+FIB_score_sort100000_100.txt', sep='\t', header=['gene1', 'gene2', 'score'])
