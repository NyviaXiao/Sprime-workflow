import pandas
import pandas as pd
import os
import sys

pop = sys.argv[1]
datapath = "/public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/{}".format(pop)
f1 = os.listdir(datapath)
for filename in f1:
    outname = os.path.join(datapath,filename)
    data = pd.read_csv(outname,sep="\t")
    chrseq = filename[7]
    print(chrseq)
    #data = pd.read_csv("out.chr1.compare_update",sep="\t")
    print(data.head())

    # 按照第一列的相同值进行分组
    grouped_data = data.groupby(data.columns[0])
    count=0
    # 打印每个分组的第一行数据
    for name, group in grouped_data:
        # 重置索引，使得每组的索引从0开始
        group = group.reset_index(drop=True)
        #print(f"Group {name}:")
        #print(type(group))#dataframe
        #break
        # 遍历第12列到最后一列
        chrseq = group.iloc[0, 2]
        for col in group.columns[11:]:
            # 初始化连续0的计数器和记录最先出现0的行和最后连续0的行的索引
            start_index = 0
            end_index = 0
            index = 0
            first_zero_row = None
            last_zero_row = None


            first_index = group[col].index[0]
            last_index = group[col].index[-1]

            index = first_index
            while index < last_index:
                if group[col][index]==0:
                    start_index = index

                    while  index < group[col].index[-1] and group[col][index + 1] == 0:
                        index += 1
                        end_index = index
                    if end_index>start_index:
                        with open(os.path.join("/public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/{}".format(pop),str(chrseq)+".txt"),"a") as f:
                            f.write(f"{group.iloc[start_index,2] } {group.iloc[start_index,3] } {group.iloc[end_index,3] } {end_index-start_index+1 } {col} {group.iloc[start_index,0]}")
                            f.write('\n')
                        # print(f"{group.iloc[start_index,2] } {group.iloc[start_index,3] } {group.iloc[end_index,3] } {end_index-start_index+1 } {col} {group.iloc[start_index,0]}")
                    index+=1
                else:
                    index+=1
#
#
#

# import pandas as pd
#
# data = pd.read_csv("out.chr1.compare_update", sep="\t")
# print(data.head())
#
# # 按照第一列的相同值进行分组
# grouped_data = data.groupby(data.columns[0])
#
# # 打印每个分组的第一行数据
# for name, group in grouped_data:
#     print(f"Group {name}:")
#
#     # 遍历第12列到最后一列
#     for col in group.columns[11:]:
#         # 初始化连续0的计数器和记录最先出现0的行和最后连续0的行的索引
#         start_index = 0
#         end_index = 0
#         index = group.index[0]
#         first_zero_row = None
#         last_zero_row = None
#
#         while index < len(group) - 1:
#             if group[col].iloc[index] == 0:
#                 start_index = index
#
#                 while index < len(group) - 1 and group[col].iloc[index + 1] == 0:
#                     index += 1
#                     end_index = index
#
#                 if end_index > start_index:
#                     print(f"chorm:{group.iloc[start_index, 2]}, start:{group.iloc[start_index, 3]}, end:{group.iloc[end_index, 3]}, length:{end_index - start_index + 1}, pop{col}")
#                 index += 1
#             else:
#                 index += 1
#
#     print("\n")

