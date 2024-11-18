import csv
import math
import statistics

def calculate_degree(csv_file):
    gene_degrees = {}
    total_degree = 0

    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)

        for row in reader:
            gene1 = row['Node1']
            gene2 = row['Node2']

            if gene1 not in gene_degrees:
                gene_degrees[gene1] = 0
            if gene2 not in gene_degrees:
                gene_degrees[gene2] = 0

            gene_degrees[gene1] += 1
            gene_degrees[gene2] += 1

            total_degree += 2  # Each edge contributes to 2 degrees (one for each node)

    return gene_degrees, total_degree


def load_gene_network(csv_file):
    gene_network = {}
    total_gene_degree = 0

    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)

        for row in reader:
            gene1 = row['Node1']
            gene2 = row['Node2']
            weight = float(row['Weight'])

            if gene1 not in gene_network:
                gene_network[gene1] = {}
            if gene2 not in gene_network:
                gene_network[gene2] = {}

            gene_network[gene1][gene2] = weight
            gene_network[gene2][gene1] = weight

            total_gene_degree += 2  # Each edge contributes to 2 degrees (one for each node)

    return gene_network, total_gene_degree


def calculate_pcc_avg_subset(gene_network, gene_subset):
    pcc_sum = 0
    count = 0

    for gene1 in gene_subset:
        for gene2 in gene_subset:
            if gene1 != gene2 and gene2 in gene_network.get(gene1, {}):
                weight = gene_network[gene1][gene2]
                pcc_sum += abs(weight)
                count += 1

    return pcc_sum / count if count != 0 else 0

def calculate_pcc_avg_neighborhood(gene_network, gene_subset):
    pcc_avg_list = []  # List to store PCC absolute average for each gene in the neighborhood

    for gene in gene_subset:
        pcc_sum = 0
        count = 0
        neighborhood = set(gene_network.get(gene, {}).keys()) - set(gene_subset)

        for neighbor in neighborhood:
            weight = gene_network[gene][neighbor]
            pcc_sum += abs(weight)
            count += 1

        # Calculate the average PCC absolute value for the current gene's neighborhood
        pcc_avg = pcc_sum / count if count != 0 else 0
        pcc_avg_list.append(pcc_avg)

    # Calculate the average of PCC absolute averages for all genes in the subset
    avg_pcc_avg = statistics.mean(pcc_avg_list) if pcc_avg_list else 0
    return avg_pcc_avg

def calculate_pcc_sum_neighborhood(gene_network, gene_subset):
    pcc_sum = 0

    for gene in gene_subset:
        neighborhood = set(gene_network.get(gene, {}).keys()) - set(gene_subset)
        for neighbor in neighborhood:
            weight = gene_network[gene][neighbor]
            pcc_sum += abs(weight)

    return pcc_sum


def calculate_pcc_avg_edges(gene_network, gene_subset):
    pcc_sum = 0
    count = 0

    for gene1 in gene_subset:
        for gene2 in gene_subset:
            if gene1 < gene2 and gene2 in gene_network.get(gene1, {}):
                weight = gene_network[gene1][gene2]
                pcc_sum += abs(weight)
                count += 1

    return pcc_sum / count if count != 0 else 0


def calculate_subset_degree(gene_network, gene_subset):
    subset_degrees = {gene: 0 for gene in gene_subset}

    for gene1 in gene_subset:
        for gene2 in gene_subset:
            if gene1 != gene2 and gene2 in gene_network.get(gene1, {}):
                subset_degrees[gene1] += 1

    return sum(subset_degrees.values())


# 计算每个基因的度
csv_file = 'data/stage3.csv'  # 你的CSV文件路径
gene_degrees, total_degree = calculate_degree(csv_file)

# 你提供的子集
gene_subset = ['PLCG1', 'FOS', 'BRCA1', 'IMPDH2', 'ATF2', 'SMAD3', 'E2F4', 'TCF3', 'GRB2', 'ATF5', 'MET', 'E2F5', 'PTGS2', 'MYC', 'ACTL6B', 'ELK1', 'GAB1', 'FAF1', 'EGR1', 'SMARCA4', 'E2F2', 'TPMT', 'SHC3', 'NFYA', 'PIK3CA', 'TP53', 'ARID1B', 'GMPS', 'MZT1', 'ARID2', 'SMARCA2', 'CTNNB1', 'FOSB']


# 加载基因网络
gene_network, total_gene_degree = load_gene_network(csv_file)

# 计算子集内的PCC绝对值的平均值
pcc_avg_subset = calculate_pcc_avg_subset(gene_network, gene_subset)
print("\nPCC Absolute Average for Gene Subset:")
print(f"{pcc_avg_subset}")

# 计算子网络与其一阶边的PCC绝对值的和
pcc_sum_neighborhood = calculate_pcc_sum_neighborhood(gene_network, gene_subset)
print("\nPCC Absolute Sum for Neighborhood of Genes in Subset:")
print(f"{pcc_sum_neighborhood}")

# 计算子集内边的PCC绝对值的平均值
pcc_avg_edges = calculate_pcc_avg_edges(gene_network, gene_subset)
print("\nPCC Absolute Average for Edges in Gene Subset:")
print(f"{pcc_avg_edges}")

# 计算子集内节点的度的总和
#subset_degree_sum = sum(gene_degrees.get(gene, 0) for gene in gene_subset)
subset_degree_sum = calculate_subset_degree(gene_network, gene_subset)
print("\nDegree Sum of Nodes in the Subset:", subset_degree_sum)

# 计算整个图
print("\nTotal Degree of All Nodes in the Graph:", total_gene_degree)

def calculate_y(a, b, c, d):
    if c == 0 or c + d == 0:
        return None  # Avoid division by zero or negative log value
    return -((d/9836) * math.log2(a/9836))

# 测试
a = subset_degree_sum
b = total_gene_degree
c = pcc_avg_subset
d = pcc_sum_neighborhood
e = pcc_avg_edges
y = calculate_y(a, b, c, d)
z = c/d
entropy = math.exp(-y)
print("y =", y)
print("z =", z)
print("entropy =", entropy)
