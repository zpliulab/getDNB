import networkx as nx
from control import read_graph_from_csv
import matplotlib.pyplot as plt
# 创建一个70节点，108条边的连通网络（这里假设你已经有了你的网络数据）
file_path = 'minist_HCC _001.csv'  # 替换为你的 CSV 文件路径
edges = read_graph_from_csv(file_path)

G_original = nx.Graph(edges)
selected_nodes = ['NFYA', 'FOS', 'PLCG1', 'SMARCA4', 'MET', 'PTGS2', 'ARID1B', 'GMPS', 'SMAD3', 'CTNNB1', 'ELK1', 'SHC3', 'E2F4', 'TP53', 'MYC', 'FOSB', 'E2F2', 'E2F5', 'BRCA1', 'TCF3', 'EGR1', 'ATF2', 'SMARCA2']
all_shortest_path_nodes = set(selected_nodes)

# 找到每对节点之间的最短路径，并将路径上的节点添加到all_shortest_path_nodes中

for i in range(len(selected_nodes)):
    for j in range(i+1, len(selected_nodes)):
        shortest_path = nx.shortest_path(G_original, source=selected_nodes[i], target=selected_nodes[j])
        all_shortest_path_nodes.update(shortest_path)

# 创建一个新图，将原始网络中的节点和最短路径上的节点添加到新图中
G_combined = nx.Graph()
G_combined.add_nodes_from(all_shortest_path_nodes)
for u, v in G_original.edges():
    if u in all_shortest_path_nodes and v in all_shortest_path_nodes:
        G_combined.add_edge(u, v)
nx.draw(G_combined, with_labels=True)
plt.show()
# 输出新图的节点和边
print("Nodes in combined graph:", G_combined.nodes())
print("Edges in combined graph:", G_combined.edges())