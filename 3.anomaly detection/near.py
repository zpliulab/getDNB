import numpy as np

def load_embeddings(file_path):
    embeddings = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            node = parts[0]
            embedding = np.array([float(x) for x in parts[1:]])
            embeddings[node] = embedding
    return embeddings

def find_closest_nodes(target_node, embeddings):
    target_embedding = embeddings[target_node]

    distances = {}
    for node, embedding in embeddings.items():
        if node != target_node:
            distance = np.linalg.norm(target_embedding - embedding)
            distances[node] = distance

    # Find the two nodes with the smallest distances
    closest_nodes = sorted(distances, key=distances.get)[:8]

    return closest_nodes

def main():
    embedding_file = 'GCN_embeddings_HCC_1.txt'
    target_node = 'GDF2'

    embeddings = load_embeddings(embedding_file)
    closest_nodes = find_closest_nodes(target_node, embeddings)

    print(f"The two closest nodes to {target_node} are: {closest_nodes}")

if __name__ == "__main__":
    main()
