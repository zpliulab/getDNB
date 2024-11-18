import numpy as np
from sklearn.metrics import pairwise_distances

def load_embeddings(file_path):
    embeddings = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            node = parts[0]
            embedding = np.array([float(x) for x in parts[1:]])
            embeddings[node] = embedding
    return embeddings

def load_cluster_centers(cluster_center_file):
    with open(cluster_center_file, 'r') as file:
        cluster_centers = file.read().splitlines()
    return cluster_centers

def calculate_distance(embedding, cluster_center):
    return np.linalg.norm(embedding - cluster_center)

def calculate_anomaly_scores(embeddings, cluster_centers):
    anomaly_scores = {}

    for node, embedding in embeddings.items():
        distances = [calculate_distance(embedding, embeddings[center]) for center in cluster_centers]
        min_distance = min(distances)
        anomaly_scores[node] = min_distance

    return anomaly_scores

def main():
    embedding_file = 'GCN_embeddings_HCC_2.txt'
    cluster_center_file = 'center_HCC_1.txt'

    embeddings = load_embeddings(embedding_file)
    cluster_centers = load_cluster_centers(cluster_center_file)

    anomaly_scores = calculate_anomaly_scores(embeddings, cluster_centers)
    # Print or save the anomaly scores
    for node, score in anomaly_scores.items():
        print(f"Node: {node}, Anomaly Score: {score}")

if __name__ == "__main__":
    main()
