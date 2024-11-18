import networkx as nx
import csv
import matplotlib.pyplot as plt


def read_graph_from_csv(file_path):
    edges = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            source, target = row
            edges.append((source, target))
    return edges


def modified_greedy_dominating_set(graph, degree_threshold=5):
    dominating_set = set()
    nodes = set(graph.nodes())

    # Add nodes with degree greater than the threshold to the dominating set
    high_degree_nodes = {node for node, degree in graph.degree if degree > degree_threshold}
    dominating_set.update(high_degree_nodes)
    nodes -= high_degree_nodes

    # Remove the neighbors of these high degree nodes from the remaining nodes
    for node in high_degree_nodes:
        nodes -= set(graph.neighbors(node))

    # Apply the greedy algorithm to cover the remaining nodes
    while nodes:
        node = max(nodes, key=lambda n: len(nodes & set(graph.neighbors(n))))
        dominating_set.add(node)
        nodes -= set([node] + list(graph.neighbors(node)))

    return dominating_set


def visualize_graph(graph, dominating_set):
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos, with_labels=True, font_weight='bold')
    nx.draw_networkx_nodes(graph, pos, nodelist=dominating_set, node_color='r')
    plt.show()


if __name__ == "__main__":
    file_path = 'minist_HCC_001.csv'  # Update the file path
    edges = read_graph_from_csv(file_path)
    G = nx.Graph(edges)
    dominating_set = modified_greedy_dominating_set(G, degree_threshold=5)  # Nodes with degree > 5
    print("最小支配集:", dominating_set)
    visualize_graph(G, dominating_set)
