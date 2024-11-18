import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 加载数据
file_path = 'class_count.csv'
df = pd.read_csv(file_path)

# 获取基因列表（排除 'name' 和 'stage' 列）
genes = df.columns[1:-1]
custom_colors = ['#DEECF6', '#AFC8E2', '#E2F2CD', '#B6DAA7', '#E8E0EF',  '#C2B1D7', '#F9D5D5', '#E89DA0']
for gene in genes:
    plt.figure(figsize=(8, 6))

    # 小提琴图
    sns.violinplot(x='stage', y=gene, data=df, palette=custom_colors)

    # 叠加平均值并用直线连接
    sns.pointplot(x='stage', y=gene, data=df, color='grey', markers='o', linestyles='-', ci=None)

    # 去掉上边框和右边框
    sns.despine(top=True, right=True)

    # 设置标题和轴标签的字体大小
    plt.title(f'Expression of {gene} across time series', fontsize=26)
    plt.ylabel('Expression Level', fontsize=24)
    plt.xlabel('Pathological time series (n)', fontsize=24)

    # 设置坐标轴刻度字体大小
    plt.xticks(rotation=0, fontsize=24)
    plt.yticks(fontsize=24)

    plt.tight_layout()

    # 保存每个基因的小提琴图
    plt.savefig(f'violin_plot_{gene}_with_mean.pdf')
    plt.close()