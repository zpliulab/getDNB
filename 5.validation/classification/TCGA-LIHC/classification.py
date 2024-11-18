import pandas as pd
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score, f1_score, precision_score, recall_score, precision_recall_curve
from imblearn.over_sampling import SMOTE
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# File names of the CSV files
file_names = ['count1.csv', 'count2.csv', 'count3.csv', 'count4.csv']

# Colors for each plot
colors = ['#DEECF6', '#AFC8E2', '#E2F2CD', '#B6DAA7',  '#E89DA0']

# Initialize lists to store AUC and AUPRC values
auc_values = []
auprc_values = []

# Initialize plots
plt.figure(figsize=(7, 6))
plt.title('AUROC Curve of Time Series', fontsize=16)
plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontsize=16)
plt.ylabel('True Positive Rate', fontsize=16)

for i, file_name in enumerate(file_names):
    # Load data
    data = pd.read_csv(file_name, index_col=0)

    # Separate features and target
    X = data.iloc[:, :-1].values  # Gene count data
    y = data.iloc[:, -1].values  # Pathological stage data

    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.35, random_state=0)

    # Train SVM classifier
    classifier = svm.SVC(kernel='linear', probability=True)
    classifier.fit(X_train, y_train)

    # Predictions on test data
    y_pred = classifier.predict(X_test)
    y_prob = classifier.predict_proba(X_test)

    # Compute evaluation metrics
    fpr, tpr, _ = roc_curve(y_test, y_prob[:, 1], pos_label=1)
    roc_auc = auc(fpr, tpr)

    precision, recall, _ = precision_recall_curve(y_test, y_prob[:, 1], pos_label=1)
    pr_auc = auc(recall, precision)

    # Store AUC and AUPRC values
    auc_values.append(roc_auc)
    auprc_values.append(pr_auc)

    # Print evaluation metrics
    print(f"File: {file_name}")
    print(f"Accuracy: {accuracy_score(y_test, y_pred):.4f}")
    print(f"F1-score: {f1_score(y_test, y_pred, pos_label=1):.4f}")
    print(f"Precision: {precision_score(y_test, y_pred, pos_label=1):.4f}")
    print(f"Recall: {recall_score(y_test, y_pred, pos_label=1):.4f}")
    print(f"AUC: {roc_auc:.4f}")
    print(f"AUPRC: {pr_auc:.4f}")

    # Interpolate the ROC curve for smoothing
    fpr_smooth = np.linspace(0, 1, 1000)  # More points for smoother curve
    tpr_interp = interp1d(fpr, tpr, kind='linear')(fpr_smooth)

    # Plot ROC curve
    plt.subplot(1, 1, 1)
    plt.plot(fpr_smooth, tpr_interp, color=colors[i], lw=3, label=f'AUROC for time {i+1} (AUC = {roc_auc:.4f})')

# Calculate the mean AUC and AUPRC
mean_auc = np.mean(auc_values)
mean_auprc = np.mean(auprc_values)

# Calculate the standard deviation of AUC and AUPRC
std_auc = np.var(auc_values)
std_auprc = np.var(auprc_values)

# Compute mean ROC curve
mean_fpr = np.linspace(0, 1, 1000)
tprs = []
for i in range(len(file_names)):
    data = pd.read_csv(file_names[i], index_col=0)
    X = data.iloc[:, :-1].values
    y = data.iloc[:, -1].values
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.35, random_state=0)
    smote = SMOTE(random_state=42)
    #X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)
    classifier = svm.SVC(kernel='linear', probability=True)
    classifier.fit(X_train, y_train)
    y_prob = classifier.predict_proba(X_test)
    fpr, tpr, _ = roc_curve(y_test, y_prob[:, 1], pos_label=1, drop_intermediate=False)
    tpr_interp = interp1d(fpr, tpr, kind='linear')(mean_fpr)
    tprs.append(tpr_interp)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[0] = 0.0
mean_tpr[-1] = 1.0

# Plot mean ROC curve
plt.subplot(1, 1, 1)
plt.plot(mean_fpr, mean_tpr, color='#E89DA0', lw=3, label=f'Mean AUROC (AUC = {mean_auc:.4f} ± {std_auc:.4f})')
plt.legend(loc="lower right", fontsize=16)

# Show plots
plt.tight_layout()
plt.savefig("roc_pr_curves.pdf", format="pdf")
plt.show()
