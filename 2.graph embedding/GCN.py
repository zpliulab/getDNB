import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
from sklearn.model_selection import train_test_split
from torch_geometric.data import Data
from torch_geometric.nn import GraphConv

# 设置随机种子
seed = 42
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
import numpy as np
np.random.seed(seed)

# Read CSV file
df = pd.read_csv('data/stage8.csv')

# Map node names to integers
node_mapping = {node: idx for idx, node in enumerate(set(df['Node1']).union(set(df['Node2'])))}
df['Node1'] = df['Node1'].map(node_mapping)
df['Node2'] = df['Node2'].map(node_mapping)

# Split into train and test sets
train_df, test_df = train_test_split(df, test_size=0.2, random_state=seed)


# Build GCN model
class GCN(nn.Module):
    def __init__(self, num_nodes, embedding_dim):
        super(GCN, self).__init__()
        self.embedding = nn.Embedding(num_nodes, embedding_dim)
        self.gc1 = GraphConv(embedding_dim, embedding_dim)
        self.relu = nn.ReLU()

    def forward(self, data):
        x = self.embedding.weight
        edge_index = data.edge_index
        edge_weight = data.edge_attr
        x = self.gc1(x, edge_index, edge_weight)
        x = self.relu(x)
        return x


# Define training function
def train(model, train_data, num_epochs=100, lr=0.01):
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss()

    for epoch in range(num_epochs):
        model.train()
        optimizer.zero_grad()

        # Create a PyTorch Geometric Data object
        edge_index = torch.tensor([train_data['Node1'].values, train_data['Node2'].values], dtype=torch.long)
        edge_weight = torch.tensor(train_data['Weight'].values, dtype=torch.float32)
        data = Data(x=model.embedding.weight, edge_index=edge_index, edge_attr=edge_weight)

        output = model(data)
        target = model.embedding.weight
        loss = criterion(output, target)
        loss.backward()
        optimizer.step()

        if epoch % 10 == 0:
            print(f'Epoch {epoch}/{num_epochs}, Loss: {loss.item()}')


# Create model and train
num_nodes = len(node_mapping)
embedding_dim = 16
model = GCN(num_nodes, embedding_dim)
train(model, train_df)

# Save node embeddings to a text file
embeddings_file = 'GCN_embeddings_HCC_8.txt'
with open(embeddings_file, 'w') as f:
    for node, idx in node_mapping.items():
        embedding = model.embedding.weight[idx].detach().numpy()
        f.write(f'{node} ' + ' '.join(map(str, embedding)) + '\n')

print(f'Node embeddings saved to {embeddings_file}')
