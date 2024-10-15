#根据元素的大小将前10%的元素置为1，后10%的元素置为-1，其余元素置为0

import pandas as pd
import numpy as np

# 读取CSV文件
matrix = pd.read_csv('./data/CCC_result.csv', index_col=0)
matrix = matrix.to_numpy(dtype='float32')

# 将矩阵的所有元素取绝对值
abs_matrix = np.abs(matrix)

# 展平矩阵并获取需要的百分位数的索引位置
flattened = abs_matrix.flatten()
n_elements = len(flattened)
top_10_percent_index = int(n_elements * 0.1)
bottom_10_percent_index = int(n_elements * 0.9)
print('n_elements=',n_elements)
print('top_10_percent_index=',top_10_percent_index)
print('bottom_10_percent_index=',bottom_10_percent_index)

# 对展平的矩阵排序
sorted_indices = np.argsort(flattened)
#print('sorted_indices:\n',sorted_indices)

# 获取前10%和后10%的元素的索引
top_indices = sorted_indices[:top_10_percent_index]
bottom_indices = sorted_indices[bottom_10_percent_index:]

# 创建一个全0的数组，与原矩阵同形状
new_matrix = np.zeros_like(abs_matrix)

# 设置前10%的元素为1，后10%的元素为-1
np.put(new_matrix, top_indices, 1)
np.put(new_matrix, bottom_indices, -1)

# 将索引从0开始转换为从1开始
row_index = np.arange(1, new_matrix.shape[0] + 1)
col_index = np.arange(1, new_matrix.shape[1] + 1)

# 创建DataFrame对象
new_matrix_df = pd.DataFrame(new_matrix, index=row_index, columns=col_index)
print('new_matrix_df:\n',new_matrix_df)


#new_matrix_df.to_csv('./data/filter_CCC_result.csv')

np.fill_diagonal(new_matrix_df.values, 0)

print("新的DataFrame:\n",new_matrix_df)
new_matrix_df.to_csv('./data/filter_CCC_result.csv')

