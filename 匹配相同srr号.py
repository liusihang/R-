import os
import pandas as pd

# 定义文件夹路径
#folder_path = r'G:\ARG输出结果\acs.est.2c08116\nofilterdata'

# 读取填空的 Excel 文件
#fill_in_file = pd.read_excel(r'G:\ARG输出结果\acs.est.2c08116\nofilterdata\PRJNA758994数据说明.xlsx', index_col=0)

# 定义文件夹路径
folder_path = r'G:\ARG输出结果\acs.est.2c08116\nofilterdata\add'

# 读取excel文件
lookup_df = pd.read_excel(r'G:\ARG输出结果\acs.est.2c08116\nofilterdata\PRJNA758994数据说明.xlsx\add', index_col=0)

# 遍历文件夹中的所有csv文件
for file_name in os.listdir(folder_path):
    if file_name.endswith('_counts.csv'):
        # 读取csv文件为DataFrame
        df = pd.read_csv(os.path.join(folder_path, file_name))
        
        # 将第一列设置为索引
        df = df.set_index(df.columns[0])
        #lookup_df = lookup_df.set_index(lookup_df.columns[0])
        
        # 根据索引进行匹配
        matched_df = lookup_df.loc[df.index]
        
        # 将csv文件和匹配到的数据合并
        #result_df = pd.merge(df, matched_df, left_index=True, right_index=True)
        
        # 将结果写入csv文件
        matched_df.to_csv(os.path.join(folder_path, file_name.split('.')[0] + '_matched.csv'))