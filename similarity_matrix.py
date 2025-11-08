import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
from sklearn.preprocessing import StandardScaler
import os
import warnings
warnings.filterwarnings('ignore')

# Set up the data directory
data_dir = "Dataset"
excel_files = [
    "Gene_expression_final (1).xlsx",
    "miRNA_target.xlsx", 
    "TCGA_meth_717samples.xlsx"
]

def load_excel_data(file_path):
    """Load Excel file and return dataframe"""
    try:
        df = pd.read_excel(file_path)
        print(f"Loaded {file_path}: Shape {df.shape}")
        return df
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

def preprocess_data(df, data_name):
    """Preprocess data for similarity calculation"""
    print(f"\nPreprocessing {data_name}...")
    print(f"Original shape: {df.shape}")
    
    # Remove non-numeric columns if any
    numeric_df = df.select_dtypes(include=[np.number])
    print(f"Numeric columns only: {numeric_df.shape}")
    
    # Handle missing values
    if numeric_df.isnull().sum().sum() > 0:
        print(f"Missing values found: {numeric_df.isnull().sum().sum()}")
        numeric_df = numeric_df.fillna(numeric_df.mean())
    
    # Standardize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(numeric_df)
    
    return scaled_data, numeric_df.columns

def compute_similarity_matrix(data, method='cosine'):
    """Compute similarity matrix using specified method"""
    if method == 'cosine':
        similarity = cosine_similarity(data)
    elif method == 'euclidean':
        # Convert distances to similarities
        distances = euclidean_distances(data)
        # Normalize distances to [0,1] range and convert to similarity
        max_dist = np.max(distances)
        similarity = 1 - (distances / max_dist)
    else:
        # Correlation similarity
        similarity = np.corrcoef(data)
    
    return similarity

def plot_similarity_matrix(similarity_matrix, title, ax):
    """Plot similarity matrix as heatmap"""
    im = ax.imshow(similarity_matrix, cmap='viridis', aspect='auto')
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('Samples')
    ax.set_ylabel('Samples')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Similarity', rotation=270, labelpad=15)
    
    return im

# Main execution
print("Loading Excel files and computing similarity matrices...\n")

# Load all datasets
datasets = {}
for file in excel_files:
    file_path = os.path.join(data_dir, file)
    if os.path.exists(file_path):
        data = load_excel_data(file_path)
        if data is not None:
            datasets[file.replace('.xlsx', '')] = data
    else:
        print(f"File not found: {file_path}")

# Process each dataset and compute similarity matrices
similarity_matrices = {}
processed_data = {}

for data_name, df in datasets.items():
    scaled_data, columns = preprocess_data(df, data_name)
    processed_data[data_name] = {'data': scaled_data, 'columns': columns}
    
    # Compute similarity matrices using different methods
    cosine_sim = compute_similarity_matrix(scaled_data, 'cosine')
    euclidean_sim = compute_similarity_matrix(scaled_data, 'euclidean')
    correlation_sim = compute_similarity_matrix(scaled_data, 'correlation')
    
    similarity_matrices[data_name] = {
        'cosine': cosine_sim,
        'euclidean': euclidean_sim,
        'correlation': correlation_sim
    }
    
    print(f"Computed similarity matrices for {data_name}")
    print(f"Similarity matrix shape: {cosine_sim.shape}")

# Create comprehensive plots
num_datasets = len(datasets)
similarity_methods = ['cosine', 'euclidean', 'correlation']

# Create output directory for images inside similarity_matrices folder
image_output_dir = os.path.join("similarity_matrices", "images")
if not os.path.exists(image_output_dir):
    os.makedirs(image_output_dir)

# Plot 1: All similarity matrices in a grid
fig, axes = plt.subplots(num_datasets, len(similarity_methods), 
                        figsize=(15, 5 * num_datasets))
fig.suptitle('Similarity Matrices for Different Data Types', fontsize=16, fontweight='bold')

if num_datasets == 1:
    axes = axes.reshape(1, -1)

row = 0
for data_name in datasets.keys():
    for col, method in enumerate(similarity_methods):
        sim_matrix = similarity_matrices[data_name][method]
        title = f"{data_name}\n{method.capitalize()} Similarity"
        plot_similarity_matrix(sim_matrix, title, axes[row, col])
    row += 1

plt.tight_layout()
# Save the comprehensive plot
plt.savefig(os.path.join(image_output_dir, 'all_similarity_matrices_overview.png'), 
            dpi=300, bbox_inches='tight')
plt.show()

# Plot 2: Individual detailed plots for each dataset
for data_name in datasets.keys():
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'Similarity Matrices for {data_name}', fontsize=14, fontweight='bold')
    
    for i, method in enumerate(similarity_methods):
        sim_matrix = similarity_matrices[data_name][method]
        
        # Create detailed heatmap
        sns.heatmap(sim_matrix, 
                   ax=axes[i],
                   cmap='viridis',
                   center=0 if method == 'correlation' else None,
                   square=True,
                   cbar_kws={'label': f'{method.capitalize()} Similarity'})
        
        axes[i].set_title(f'{method.capitalize()} Similarity')
        axes[i].set_xlabel('Samples')
        axes[i].set_ylabel('Samples')
    
    plt.tight_layout()
    # Save individual detailed plot for each dataset
    clean_name = data_name.replace(' ', '_').replace('(', '').replace(')', '')
    filename = f'{clean_name}_detailed_similarity_matrices.png'
    plt.savefig(os.path.join(image_output_dir, filename), 
                dpi=300, bbox_inches='tight')
    plt.show()

# Summary statistics for similarity matrices
print("\n" + "="*60)
print("SIMILARITY MATRIX STATISTICS")
print("="*60)

for data_name in datasets.keys():
    print(f"\n{data_name.upper()}:")
    print("-" * 40)
    
    for method in similarity_methods:
        sim_matrix = similarity_matrices[data_name][method]
        
        # Remove diagonal elements for statistics (self-similarity = 1)
        off_diagonal = sim_matrix[np.triu_indices_from(sim_matrix, k=1)]
        
        print(f"{method.capitalize()} Similarity:")
        print(f"  Mean: {np.mean(off_diagonal):.4f}")
        print(f"  Std:  {np.std(off_diagonal):.4f}")
        print(f"  Min:  {np.min(off_diagonal):.4f}")
        print(f"  Max:  {np.max(off_diagonal):.4f}")
        print()

# Save similarity matrices as CSV files
output_dir = "similarity_matrices"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for data_name in datasets.keys():
    for method in similarity_methods:
        sim_matrix = similarity_matrices[data_name][method]
        filename = f"{data_name}_{method}_similarity.csv"
        filepath = os.path.join(output_dir, filename)
        pd.DataFrame(sim_matrix).to_csv(filepath, index=False)
        print(f"Saved: {filepath}")

print(f"\n✓ Images saved in 'similarity_matrices/images/' folder:")
print(f"  - all_similarity_matrices_overview.png (comprehensive overview)")
for data_name in datasets.keys():
    clean_name = data_name.replace(' ', '_').replace('(', '').replace(')', '')
    print(f"  - {clean_name}_detailed_similarity_matrices.png")

print("\nAnalysis complete! Similarity matrices have been computed and visualized.")
