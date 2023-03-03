import os
import zipfile

# 获取当前路径
path = os.getcwd()

# 压缩文件名
zip_filename = "data.zip"

# 遍历当前路径下的所有子文件夹
for root, dirs, files in os.walk(path):
    for file in files:
        # 提取txt和csv文件
        if file.endswith(".txt") or file.endswith(".csv"):
            # 获取文件绝对路径
            file_path = os.path.join(root, file)
            # 添加到zip文件中
            with zipfile.ZipFile(zip_filename, "a") as zip:
                zip.write(file_path, arcname=file)


# 关闭压缩文件
zip_file.close()



