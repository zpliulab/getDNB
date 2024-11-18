import pandas as pd
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score, f1_score, precision_score, recall_score, precision_recall_curve
from imblearn.over_sampling import SMOTE
import matplotlib.pyplot as plt
import numpy as np

# File names of the CSV files
file_names = ['count1.csv', 'count2.csv', 'count3.csv', 'count4.csv', 'count5.csv', 'count6.csv']

# Initialize lists to store evaluation metrics
auc_values = []
auprc_values = []
accuracy_values = []
f1_values = []
precision_values = []
recall_values = []

for file_name in file_names:
    # Load data
    data = pd.read_csv(file_name, index_col=0)

    # Separate features and target
    X = data.iloc[:, :-1].values  # Gene count data
    y = data.iloc[:, -1].values  # Pathological stage data

    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=0)

    # Apply SMOTE to the training data
    smote = SMOTE(random_state=42)
    X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)

    # Train SVM classifier
    classifier = svm.SVC(kernel='linear', probability=True)
    classifier.fit(X_train_resampled, y_train_resampled)

    # Predictions on test data
    y_pred = classifier.predict(X_test)
    y_prob = classifier.predict_proba(X_test)

    # Compute evaluation metrics
    fpr, tpr, _ = roc_curve(y_test, y_prob[:, 1], pos_label=1)
    roc_auc = auc(fpr, tpr)
    precision, recall, _ = precision_recall_curve(y_test, y_prob[:, 1], pos_label=1)
    pr_auc = auc(recall, precision)
    acc = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, pos_label=1)
    prec = precision_score(y_test, y_pred, pos_label=1)
    rec = recall_score(y_test, y_pred, pos_label=1)

    # Store evaluation metrics
    auc_values.append(roc_auc)
    auprc_values.append(pr_auc)
    accuracy_values.append(acc)
    f1_values.append(f1)
    precision_values.append(prec)
    recall_values.append(rec)

    # Print evaluation metrics
    print(f"File: {file_name}")
    print(f"Accuracy: {acc:.4f}")
    print(f"F1-score: {f1:.4f}")
    print(f"Precision: {prec:.4f}")
    print(f"Recall: {rec:.4f}")
    print(f"AUC: {roc_auc:.4f}")
    print(f"AUPRC: {pr_auc:.4f}")


metrics = {
    'AUC': auc_values,
    'AUPRC': auprc_values,
    'Accuracy': accuracy_values,
    'F1-score': f1_values,
    'Precision': precision_values,
    'Recall': recall_values
}

means = {metric: np.mean(values) for metric, values in metrics.items()}
stds = {metric: np.var(values) for metric, values in metrics.items()}

# Print mean and std for each metric
for metric in metrics:
    mean_value = means[metric]
    std_value = stds[metric]
    print(f"{metric}: Mean = {mean_value:.4f}, Std = {std_value:.4f}")

means = {metric: np.mean(values) for metric, values in metrics.items()}
stds = {metric: np.var(values) for metric, values in metrics.items()}

# Append mean values to the lists
for metric, values in metrics.items():
    values.append(means[metric])
    stds[metric] = [stds[metric]] * len(file_names) + [stds[metric]]  # Use std for each file's std, and mean's std

# Define new x-axis labels
x_labels = [f'Time {i + 1}' for i in range(len(file_names))]
x_labels.append('Means')

# Define colors for the bars
colors = ['#DEECF6', '#AFC8E2', '#E2F2CD', '#B6DAA7', '#E8E0EF',  '#C2B1D7', '#E89DA0']

# Function to plot each metric
def plot_metric(metric_name, values, stds):
    plt.figure(figsize=(8, 6))
    bars = plt.bar(range(len(x_labels)), values, yerr=stds, capsize=5, color=colors, width=0.6)
    plt.xticks(range(len(x_labels)), x_labels, rotation=0, ha='center', fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(f'{metric_name} Scores for time series', fontsize=18, pad=60)  # Adjust the pad parameter to move the title higher
    plt.ylabel(metric_name, fontsize=18)
    plt.ylim([0, 1])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    for bar, val, std in zip(bars, values, stds):
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 0.08, f'{round(yval, 4)}', ha='center', fontsize=13)
    plt.tight_layout()
    plt.savefig(f'{metric_name}_scores.pdf')
    plt.show()

# Plot and save each metric
for metric, values in metrics.items():
    plot_metric(metric, values, stds[metric])
