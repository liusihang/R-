import os
import pandas as pd

# 定义文件夹路径
folder_path = r'G:\ARG输出结果\acs.est.2c08116\nofilterdata'

# 遍历文件夹中的所有文件
for file_name in os.listdir(folder_path):
    # 判断文件名是否以 _counts.csv 结尾
    if file_name.endswith('_counts.csv'):
        # 读取 CSV 文件为 DataFrame
        df = pd.read_csv(os.path.join(folder_path, file_name))

        # 删除第一列中的字符串
        df.iloc[:, 0] = df.iloc[:, 0].str.replace('_1.fastq.gz_paired.fq.gz', '')

        # 将 DataFrame 写入 CSV 文件
        df.to_csv(os.path.join(folder_path, file_name), index=False)
