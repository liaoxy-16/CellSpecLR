import pandas as pd


#df1 = pd.read_csv('./data_preprocess/filter_fish_label.txt', delimiter='\t', names=['number', 'label'])
df1 = pd.read_csv('./data_preprocess/filter_fish_label.txt', delimiter='\t')
df1.index = range(1, len(df1.index) + 1)
print('df1:\n', df1)
# df1 = df1.iloc[0:100,:]
# print('df1:\n', df1)

# 计算每个label对应的行数
label_counts = df1['label'].value_counts().reset_index()
label_counts.columns = ['label', 'count']
print('label_counts:\n', label_counts)


# 保存结果到新的CSV文件
label_counts.to_csv('./data/filter_fish_Cell_label_count.csv', header=True)



# 细胞组间通讯
# 读入第一个CSV文件，获取label列的信息
# df1 = pd.read_csv("./data/LS_Cell_label_count.csv",index_col=0)
# print('df1:\n',df1)
label_indices = label_counts["label"].unique().tolist()
print('label_indices:',label_indices)
# 创建全0矩阵，行列索引都是label列
# labels = sorted(list(set(label_indices)))  # 获取唯一的label值并排序
# print('labels:',labels)
ccc_matrix = pd.DataFrame(0, index=label_indices, columns=label_indices)  # cell cell communication

# 读入第二个CSV文件，遍历ab两列，对应位置加1
df2 = pd.read_csv("./data/filter_fish_cell_communication_label_list.txt",delimiter='\t')
print('df2:\n',df2)
for index, row in df2.iterrows():
    #print('index=',index)
    #print('row=\n',row)
    a = row["label_row"]
    b = row["label_column"]
    #print('a=',a,' b=',b)
    #if a in label_indices and b in label_indices:
    ccc_matrix.at[a, b] += 1

print(ccc_matrix)

# 将新矩阵保存到CSV文件
ccc_matrix.to_csv("./data/filter_fish_cell_label_communication.csv")

#--------------------------------------------------------------------------
label_counts = pd.read_csv('./data/filter_fish_Cell_label_count.csv', index_col=0)
print('label_counts:\n',label_counts)
# 将label_counts中的'label'列设置为索引
label_counts.set_index('label', inplace=True)
#label_counts.index = label_counts.index.map(float)
print('label_counts set_index:\n',label_counts)

ccc_matrix = pd.read_csv('./data/filter_fish_cell_label_communication.csv', index_col=0)
print('ccc_matrix:\n',ccc_matrix)

# 对ccc_matrix 中的每个元素进行操作，除以行列索引对应的count值相加后的值
#ccc_matrix_new = ccc_matrix.copy()  # 创建一个新的DataFrame，用于存储操作后的结果
ccc_matrix_new = pd.DataFrame(0.0, index=ccc_matrix.index, columns=ccc_matrix.columns)
print('ccc_matrix_new:\n',ccc_matrix_new)

for i, row in ccc_matrix.iterrows():
    print('i:\n', i)
    print('row:\n', row)
    for j, value in row.items():
        print('j:\n', j)
        print('value:\n', value)
        count_sum = label_counts.loc[float(i), 'count'] * label_counts.loc[float(j), 'count']  # 获取对应的count值相加的结果
        new_value = value / count_sum  # 计算新的值
        ccc_matrix_new.loc[i, j] = new_value  # 将新的值赋值给新的DataFrame


print('ccc_matrix_new:\n',ccc_matrix_new)


ccc_matrix_new.to_csv("./data/filter_fish_cell_label_communication_proportion.csv")