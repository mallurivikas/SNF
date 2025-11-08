import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import SpectralClustering
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics import adjusted_rand_score, silhouette_score, normalized_mutual_info_score, accuracy_score
from sklearn.preprocessing import LabelEncoder
import os
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("PHASE 1: DATA PREPARATION & GROUND TRUTH DISCOVERY")
print("="*70)

# Load the Excel files to find cancer labels
data_dir = "Dataset"
excel_files = [
    "Gene_expression_final (1).xlsx",
    "miRNA_target.xlsx", 
    "TCGA_meth_717samples.xlsx"
]

print("Looking for cancer subtype labels in the data...")

# Check each file for potential label columns
cancer_labels = None
label_source = None

for file in excel_files:
    file_path = os.path.join(data_dir, file)
    if os.path.exists(file_path):
        print(f"\nChecking {file}...")
        
        # Load the file
        df = pd.read_excel(file_path)
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {list(df.columns[:10])}...")  # Show first 10 columns
        
        # Look for columns that might contain cancer labels
        potential_label_cols = []
        for col in df.columns:
            col_lower = str(col).lower()
            if any(word in col_lower for word in ['cancer', 'subtype', 'type', 'class', 'label', 'group', 'cluster']):
                potential_label_cols.append(col)
        
        if potential_label_cols:
            print(f"  Potential label columns found: {potential_label_cols}")
            for col in potential_label_cols:
                unique_vals = df[col].nunique()
                sample_vals = df[col].value_counts().head()
                print(f"    {col}: {unique_vals} unique values")
                print(f"    Sample values: {dict(sample_vals)}")
        else:
            print(f"  No obvious label columns found")

# If no labels found in main files, check if there's a separate labels file
print(f"\nLooking for separate label files in Dataset directory...")
all_files = os.listdir(data_dir) if os.path.exists(data_dir) else []
for file in all_files:
    if not file in excel_files:
        print(f"  Found additional file: {file}")

# Extract REAL cancer subtype labels from the Subtype-3 column
print(f"\n" + "="*50)
print("EXTRACTING REAL CANCER SUBTYPE LABELS")
print("="*50)

# Load the first dataset to get the actual cancer labels
gene_expr_file = os.path.join(data_dir, "Gene_expression_final (1).xlsx")
if os.path.exists(gene_expr_file):
    gene_df = pd.read_excel(gene_expr_file)
    n_samples = gene_df.shape[0]
    print(f"Total samples: {n_samples}")
    
    # Extract the actual cancer subtype labels from 'Subtype-3' column
    if 'Subtype-3' in gene_df.columns:
        subtype_col = gene_df['Subtype-3']
        print(f"\n✓ Found 'Subtype-3' column!")
        print(f"Unique cancer subtypes: {subtype_col.unique()}")
        print(f"Subtype distribution:")
        for subtype, count in subtype_col.value_counts().items():
            print(f"  {subtype}: {count} samples")
        
        # Convert to numeric labels
        from sklearn.preprocessing import LabelEncoder
        label_encoder = LabelEncoder()
        true_labels = label_encoder.fit_transform(subtype_col)
        
        # Store the mapping
        label_mapping = dict(zip(label_encoder.classes_, range(len(label_encoder.classes_))))
        print(f"\nLabel mapping: {label_mapping}")
        
        cancer_labels = true_labels
        label_source = "Subtype-3 column (real cancer subtypes)"
        print(f"✓ Using REAL cancer subtype labels from data")
        print(f"Numeric label distribution: {np.bincount(true_labels)}")
        
    else:
        print("❌ 'Subtype-3' column not found!")
        cancer_labels = None
        label_source = None
else:
    print("❌ Gene expression file not found!")
    cancer_labels = None
    label_source = None

# Save the REAL cancer subtype labels
if cancer_labels is not None:
    # Create the Actual Labels folder
    actual_labels_folder = "Actual SubType"
    os.makedirs(actual_labels_folder, exist_ok=True)
    
    # Also save the original subtype names for reference
    labels_df = pd.DataFrame({
        'sample_index': range(len(cancer_labels)),
        'subtype_name': gene_df['Subtype-3'].values if 'Subtype-3' in gene_df.columns else ['Unknown']*len(cancer_labels),
        'subtype_numeric': cancer_labels
    })
    
    # Save to the Actual Labels folder
    output_path = os.path.join(actual_labels_folder, 'real_cancer_subtype_labels.csv')
    labels_df.to_csv(output_path, index=False)
    print(f"\n✓ Saved REAL cancer subtype labels to: {output_path}")
    print(f"  Source: {label_source}")
    print(f"  Number of subtypes: {len(np.unique(cancer_labels))}")
    
    # Show the mapping between names and numbers
    if 'Subtype-3' in gene_df.columns:
        for name, num in zip(gene_df['Subtype-3'].unique(), 
                           [np.where(cancer_labels == i)[0][0] for i in np.unique(cancer_labels)]):
            num_samples = np.sum(cancer_labels == label_encoder.transform([name])[0])
            print(f"  {name} -> {label_encoder.transform([name])[0]} ({num_samples} samples)")
else:
    print("❌ Could not extract cancer subtype labels!")

print(f"\n" + "="*70)
print("PHASE 1 COMPLETE - REAL CANCER LABELS EXTRACTED")
print("="*70)