# getDNB
**getDNB: Identifying Dynamic Network Biomarkers of Hepatocellular Carcinoma from Time-Varying Gene Regulations Utilizing Graph Embedding Techniques for Anomaly Detection**

We proposed a computational strategy named getDNB to discover dynamic network biomarker (DNB) by using graph embedding technology (get) for anomaly detection
## Requirements
* python>=3.9
* R version 4.3.3
## Files:
* 1. data process: The data process contains the original data of data sets GSE6764, GSE89377 and TCGA-LIHC, data processing programs and output results. In data processing programs, it includes the normalization of data, the acquisition of network weights, and the construction of specific temporal networks.

*** 2. graph embedding: **The graph embedding contains the raw data, programs, and outputs used by the graph convolution neural network.

*** 3. anomaly detection: **The anomaly detection contains the input data, anomaly detection procedures and output results. Here K-means.py is used to cluster by using node vector representation to obtain cluster center. anomaly detection.py is used to calculate the node deviation score. near.py is used to get multiple nodes that are closest to a node.

*** 4. crucial point: **The crucial point contains the input abnormal gene set network, the key point acquisition program, the output DNBs and the connected network of DNBs.
