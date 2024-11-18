import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score


with open('GCN_embeddings_HCC_1.txt', 'r') as file:
    lines = file.readlines()

data = [line.strip().split(' ') for line in lines]

embeddings = np.array([list(map(float, entry[1:])) for entry in data])

range_n_clusters = range(2, 30)

# 计算肘部法的SSD
#ssd = []
#for num_clusters in range_n_clusters:
#    kmeans = KMeans(n_clusters=num_clusters, n_init=1000, random_state=42)
#    kmeans.fit(embeddings)
 #   ssd.append(kmeans.inertia_)

# 绘制SSD曲线
#plt.figure(figsize=(10, 6))
#plt.plot(range_n_clusters, ssd, marker='o')
#plt.xlabel('Number of clusters')
#plt.ylabel('Sum of Squared Distances (SSD)')
#plt.title('Elbow Method For Optimal k')
#plt.show()

# 计算轮廓系数
#silhouette_avg = []
#for num_clusters in range_n_clusters:
#   kmeans = KMeans(n_clusters=num_clusters, n_init=1000, random_state=42)
#   cluster_labels = kmeans.fit_predict(embeddings)
#   silhouette_avg.append(silhouette_score(embeddings, cluster_labels))

# 绘制轮廓系数曲线
#plt.figure(figsize=(10, 6))
#plt.plot(range_n_clusters, silhouette_avg, marker='o')
#plt.xlabel('Number of clusters')
#plt.ylabel('Average Silhouette Score')
#plt.title('Silhouette Method For Optimal k')
#plt.show()

num_clusters = 24
kmeans = KMeans(n_clusters=num_clusters, n_init=1000, random_state=0)
kmeans.fit(embeddings)

cluster_centers = kmeans.cluster_centers_

for center in cluster_centers:
    idx = np.argmin(np.linalg.norm(embeddings - center, axis=1))
    node_name = data[idx][0]
    print(f"Node name: {node_name}")
