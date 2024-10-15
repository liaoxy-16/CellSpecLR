# CCC_list.csv：记录有细胞通讯的两个细胞的索引和类型，有四列，第一列是发送细胞的索引，
# 第二列是接收细胞的索引，第三列是发送细胞的类型，第四列是接收细胞的类型。
#——————————————————————————————————————————————————————————————————————————————————

import pandas as pd

# 读取CSV文件
df = pd.read_csv('./data/filter_CCC_result.csv',index_col=0)
print('df:\n',df)

# 计算每个值的数量
count_0 = df[df == 0].count().sum()
count_1 = df[df == 1].count().sum()
count_minus_1 = df[df == -1].count().sum()

# 打印结果
print('数量为0的元素个数:', count_0)
print('数量为1的元素个数:', count_1)
print('数量为-1的元素个数:', count_minus_1)

#------------------------------------------------------------
# 查找元素为1的行列索引
indices_1 = df[df == 1].stack().reset_index()
indices_1.columns = ['Row', 'Column','Communication']
print('indices_1:\n',indices_1)

# 获取行数
num_rows = indices_1.shape[0]
# 打印行数
print("行数:", num_rows)

# 保存元素为1的行列索引到新的CSV文件
indices_1.to_csv('./data/cell_communication_list.txt', index=False,sep='\t')

#-----------------------------------------------------------------------------
# LS_cell_communication_label_list.csv：

#df1 = pd.read_csv('./data_preprocess/filter_fish_label.txt', delimiter='\t',names=['index', 'label'])
df1 = pd.read_csv('./data_preprocess/data_label.txt', delimiter='\t')
print('df1:\n',df1)
# df1 = df1.iloc[0:100,:]
# print('df1:\n', df1)

# 读取第二个CSV文件的'Row'列和'Column'列
df2 = pd.read_csv('./data/cell_communication_list.txt', usecols=['Row', 'Column'],delimiter='\t')
print('df2:\n',df2)

# 将第二个文件的'Row'列与第一个文件进行合并，生成对应的'label'列
merged_df = pd.merge(df2, df1, left_on='Row', right_on='index')
merged_df.rename(columns={'label': 'label_row'}, inplace=True)
print('merged_df:\n',merged_df)

# 将第二个文件的'Column'列与第一个文件进行合并，生成对应的'label'列
merged_df = pd.merge(merged_df, df1, left_on='Column', right_on='index')
merged_df.rename(columns={'label': 'label_column'}, inplace=True)
print('merged_df:\n',merged_df)

# 保存合并后的结果到新的CSV文件
merged_df[['Row','Column','label_row','label_column']].\
    to_csv('./data/CCC_list.txt', index=False,sep='\t')


