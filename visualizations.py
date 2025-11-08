import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.neighbors import kneighbors_graph
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("3D CLUSTER VISUALIZATION FOR CANCER SUBTYPES")
print("="*70)

# Load the best results from our analysis
results_path = os.path.join('clustering results', 'clustering_results_summary.csv')
results_df = pd.read_csv(results_path)
labels_path = os.path.join('Actual SubType', 'real_cancer_subtype_labels.csv')
true_labels_df = pd.read_csv(labels_path)

# Get the best method for each data type
best_methods = {}
for data_type in ['Gene_Expression', 'miRNA', 'Methylation', 'Fused']:
    best = results_df[results_df['data_type'] == data_type].iloc[0]
    best_methods[data_type] = {
        'method': best['method'],
        'network_type': best['network_type'],
        'k': int(best['best_k'])
    }

print("🔍 Best methods for each data type:")
for data_type, info in best_methods.items():
    print(f"  {data_type}: {info['method']} ({info['network_type']}, k={info['k']})")

# Load similarity matrices
matrices_data = {
    'Gene_Expression': 'similarity_matrices/Gene_expression_final (1)_euclidean_similarity.csv',
    'miRNA': 'similarity_matrices/miRNA_target_euclidean_similarity.csv',
    'Methylation': 'similarity_matrices/TCGA_meth_717samples_euclidean_similarity.csv',
    'Fused': 'fused_matrices/euclidean_fused_average.csv'
}

similarity_matrices = {}
for name, filepath in matrices_data.items():
    if os.path.exists(filepath):
        matrix = pd.read_csv(filepath).values
        similarity_matrices[name] = matrix
        print(f"✓ Loaded {name}: {matrix.shape}")

# Network construction functions (copied from phase2)
def create_knn_network(similarity_matrix, k):
    """Create KNN network from similarity matrix"""
    distance_matrix = 1 - similarity_matrix
    np.fill_diagonal(distance_matrix, 0)
    knn_graph = kneighbors_graph(distance_matrix, n_neighbors=k, 
                                metric='precomputed', mode='connectivity')
    knn_graph = (knn_graph + knn_graph.T) / 2
    return knn_graph.toarray()

def create_percentile_network(similarity_matrix, percentile):
    """Create network using percentile threshold"""
    threshold = np.percentile(similarity_matrix, percentile)
    network = (similarity_matrix >= threshold).astype(float)
    np.fill_diagonal(network, 0)
    return network

def create_consensus_network(similarity_matrix):
    """Create consensus network using geometric mean approach"""
    sim_no_diag = similarity_matrix.copy()
    np.fill_diagonal(sim_no_diag, 0)
    geometric_mean = np.exp(np.mean(np.log(sim_no_diag[sim_no_diag > 0])))
    network = (similarity_matrix >= geometric_mean).astype(float)
    np.fill_diagonal(network, 0)
    return network

def make_weighted_network(binary_network, similarity_matrix):
    """Convert binary network to weighted using original similarities"""
    weighted = binary_network * similarity_matrix
    np.fill_diagonal(weighted, 0)
    return weighted

def get_network_for_method(similarity_matrix, method_name, network_type):
    """Create network based on method name and type"""
    if 'KNN' in method_name:
        k = int(method_name.split('-')[1])
        binary_network = create_knn_network(similarity_matrix, k)
    elif 'Percentile' in method_name:
        percentile = int(method_name.split('-')[1])
        binary_network = create_percentile_network(similarity_matrix, percentile)
    elif 'Consensus' in method_name:
        binary_network = create_consensus_network(similarity_matrix)
    else:
        raise ValueError(f"Unknown method: {method_name}")
    
    if network_type == 'Weighted':
        return make_weighted_network(binary_network, similarity_matrix)
    else:
        return binary_network

# Function to perform clustering and get labels
def get_cluster_labels(data_type):
    """Get cluster labels for a specific data type using best method"""
    matrix = similarity_matrices[data_type]
    method_info = best_methods[data_type]
    
    # Create network
    network = get_network_for_method(matrix, method_info['method'], method_info['network_type'])
    
    # Perform spectral clustering
    spectral = SpectralClustering(n_clusters=method_info['k'], affinity='precomputed', random_state=42)
    cluster_labels = spectral.fit_predict(network)
    
    return cluster_labels

# Get all cluster labels
print(f"\n🔄 Generating cluster labels for each data type...")
cluster_labels = {}
for data_type in ['Gene_Expression', 'miRNA', 'Methylation', 'Fused']:
    labels = get_cluster_labels(data_type)
    cluster_labels[data_type] = labels
    print(f"  {data_type}: {len(np.unique(labels))} clusters")

# Get true labels
true_labels = true_labels_df['subtype_numeric'].values
print(f"  True Labels: {len(np.unique(true_labels))} cancer subtypes")

# Use PCA to get 3D representation for visualization
print(f"\n🎨 Creating 3D representations using PCA...")

# Use the fused matrix as the base for 3D coordinates (most comprehensive)
pca_3d = PCA(n_components=3, random_state=42)
coords_3d = pca_3d.fit_transform(similarity_matrices['Fused'])

explained_variance = pca_3d.explained_variance_ratio_
print(f"PCA Explained Variance: {explained_variance.sum():.1%} total")
print(f"  PC1: {explained_variance[0]:.1%}, PC2: {explained_variance[1]:.1%}, PC3: {explained_variance[2]:.1%}")

# Define colors for clusters and subtypes
colors_clusters = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7']
colors_true = {'HER2': '#FF6B6B', 'LUM': '#4ECDC4', 'TNBC': '#45B7D1'}

# Create 3D visualizations
fig = plt.figure(figsize=(20, 16))

# Plot configurations
plots_config = [
    ('Gene Expression', cluster_labels['Gene_Expression'], 1),
    ('miRNA', cluster_labels['miRNA'], 2),
    ('Methylation', cluster_labels['Methylation'], 3),
    ('Fused Matrix', cluster_labels['Fused'], 4),
    ('True Cancer Subtypes', true_labels, 5)
]

for i, (title, labels, subplot_num) in enumerate(plots_config):
    ax = fig.add_subplot(2, 3, subplot_num, projection='3d')
    
    # Get unique labels and assign colors
    unique_labels = np.unique(labels)
    
    if title == 'True Cancer Subtypes':
        # Use cancer subtype names and specific colors
        subtype_names = true_labels_df['subtype_name'].values
        for label in unique_labels:
            mask = labels == label
            subtype_name = subtype_names[mask][0]  # Get the subtype name
            color = colors_true.get(subtype_name, colors_clusters[label])
            
            ax.scatter(coords_3d[mask, 0], coords_3d[mask, 1], coords_3d[mask, 2],
                      c=color, alpha=0.7, s=30, label=f'{subtype_name}')
    else:
        # Use cluster numbers
        for j, label in enumerate(unique_labels):
            mask = labels == label
            ax.scatter(coords_3d[mask, 0], coords_3d[mask, 1], coords_3d[mask, 2],
                      c=colors_clusters[j], alpha=0.7, s=30, label=f'Cluster {label}')
    
    # Customize the plot
    ax.set_xlabel(f'PC1 ({explained_variance[0]:.1%})', fontsize=10)
    ax.set_ylabel(f'PC2 ({explained_variance[1]:.1%})', fontsize=10)
    ax.set_zlabel(f'PC3 ({explained_variance[2]:.1%})', fontsize=10)
    ax.set_title(f'{title}\n3D Cluster Visualization', fontsize=12, fontweight='bold', pad=20)
    
    # Add legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    
    # Set viewing angle for better visualization
    ax.view_init(elev=20, azim=45)

# Add an overall title
fig.suptitle('Cancer Subtype Discovery: 3D Cluster Visualizations\nUsing PCA Projection of Similarity Networks', 
             fontsize=16, fontweight='bold', y=0.95)

# Create output directory for cluster visualization images
viz_output_dir = "cluster visualizations"
if not os.path.exists(viz_output_dir):
    os.makedirs(viz_output_dir)

plt.tight_layout()
plt.subplots_adjust(top=0.88)
plt.savefig(os.path.join(viz_output_dir, '3D_cluster_visualizations.png'), dpi=300, bbox_inches='tight')
plt.show()

# Create individual detailed 3D plots for better viewing
print(f"\n📊 Creating individual detailed 3D plots...")

# Individual plots with better detail
for title, labels, _ in plots_config:
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    unique_labels = np.unique(labels)
    
    if title == 'True Cancer Subtypes':
        subtype_names = true_labels_df['subtype_name'].values
        for label in unique_labels:
            mask = labels == label
            subtype_name = subtype_names[mask][0]
            color = colors_true.get(subtype_name, colors_clusters[label])
            
            ax.scatter(coords_3d[mask, 0], coords_3d[mask, 1], coords_3d[mask, 2],
                      c=color, alpha=0.8, s=50, label=f'{subtype_name} (n={np.sum(mask)})', 
                      edgecolors='black', linewidth=0.5)
    else:
        for j, label in enumerate(unique_labels):
            mask = labels == label
            ax.scatter(coords_3d[mask, 0], coords_3d[mask, 1], coords_3d[mask, 2],
                      c=colors_clusters[j], alpha=0.8, s=50, label=f'Cluster {label} (n={np.sum(mask)})', 
                      edgecolors='black', linewidth=0.5)
    
    # Enhanced styling
    ax.set_xlabel(f'Principal Component 1 ({explained_variance[0]:.1%} variance)', fontsize=12)
    ax.set_ylabel(f'Principal Component 2 ({explained_variance[1]:.1%} variance)', fontsize=12)
    ax.set_zlabel(f'Principal Component 3 ({explained_variance[2]:.1%} variance)', fontsize=12)
    
    # Method info for clustering results
    if title != 'True Cancer Subtypes':
        # Map display title to data type key
        title_to_key = {
            'Gene Expression': 'Gene_Expression',
            'miRNA': 'miRNA', 
            'Methylation': 'Methylation',
            'Fused Matrix': 'Fused'
        }
        data_type_key = title_to_key.get(title, title)
        method_info = best_methods[data_type_key]
        title_extended = f"{title} Clustering\nMethod: {method_info['method']} ({method_info['network_type']} Network)"
    else:
        title_extended = f"{title}\nGround Truth Labels"
    
    ax.set_title(title_extended, fontsize=14, fontweight='bold', pad=20)
    
    # Legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    
    # Better viewing angle
    ax.view_init(elev=15, azim=45)
    
    # Grid
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save individual plots in cluster visualizations folder
    safe_title = title.replace(' ', '_').replace('/', '_')
    plt.savefig(os.path.join(viz_output_dir, f'3D_{safe_title}_clusters.png'), dpi=300, bbox_inches='tight')
    plt.show()

print(f"\n✅ 3D VISUALIZATION COMPLETE!")
print("="*40)
print(f"📁 Files created in 'cluster visualizations/' folder:")
print("  ✓ 3D_cluster_visualizations.png (all 5 plots together)")
print("  ✓ 3D_Gene_Expression_clusters.png")
print("  ✓ 3D_miRNA_clusters.png") 
print("  ✓ 3D_Methylation_clusters.png")
print("  ✓ 3D_Fused_Matrix_clusters.png")
print("  ✓ 3D_True_Cancer_Subtypes_clusters.png")

print(f"\n🎯 SUMMARY:")
print("Generated 3D visualizations for all clustering approaches:")
print("• Gene Expression clustering (best method)")
print("• miRNA clustering (best method)")  
print("• Methylation clustering (best method)")
print("• Fused Matrix clustering (best method)")
print("• True cancer subtypes (ground truth)")

print(f"\nAll clusters are projected into 3D space using PCA for visualization.")
print(f"Total variance explained by 3 components: {explained_variance.sum():.1%}")
print("Ready for your next instructions! 🚀")