import pandas as pd
import numpy as np
from sklearn.cluster import SpectralClustering
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics import adjusted_rand_score, silhouette_score, normalized_mutual_info_score
from sklearn.preprocessing import LabelEncoder
import os
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("PHASE 2: NETWORK CONSTRUCTION & SPECTRAL CLUSTERING")
print("="*70)

# Load REAL cancer subtype labels from Phase 1
print("Loading REAL cancer subtype labels...")
labels_file = os.path.join('Actual SubType', 'real_cancer_subtype_labels.csv')

if os.path.exists(labels_file):
    labels_df = pd.read_csv(labels_file)
    true_labels = labels_df['subtype_numeric'].values
    subtype_names = labels_df['subtype_name'].values
    
    print(f"✓ Loaded REAL cancer subtypes:")
    unique_subtypes, counts = np.unique(subtype_names, return_counts=True)
    for subtype, count in zip(unique_subtypes, counts):
        numeric_label = labels_df[labels_df['subtype_name'] == subtype]['subtype_numeric'].iloc[0]
        print(f"  {subtype}: {count} samples (encoded as {numeric_label})")
    
    print(f"\nTotal samples with labels: {len(true_labels)}")
else:
    print(f"❌ Real cancer labels file not found at: {labels_file}")
    print("   Run extract_labels.py first to create the file in 'Actual Labels' folder.")
    exit()

# Load similarity matrices
print(f"\nLoading similarity matrices...")
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

# Network construction functions
def create_knn_network(similarity_matrix, k):
    """Create KNN network from similarity matrix"""
    # Convert similarity to distance for KNN (higher similarity = lower distance)
    distance_matrix = 1 - similarity_matrix
    np.fill_diagonal(distance_matrix, 0)
    
    # Create KNN graph
    knn_graph = kneighbors_graph(distance_matrix, n_neighbors=k, 
                                metric='precomputed', mode='connectivity')
    
    # Make symmetric
    knn_graph = (knn_graph + knn_graph.T) / 2
    return knn_graph.toarray()

def create_percentile_network(similarity_matrix, percentile):
    """Create network using percentile threshold"""
    threshold = np.percentile(similarity_matrix, percentile)
    network = (similarity_matrix >= threshold).astype(float)
    np.fill_diagonal(network, 0)  # Remove self-connections
    return network

def create_consensus_network(similarity_matrix):
    """Create consensus network using geometric mean approach"""
    # Remove diagonal
    sim_no_diag = similarity_matrix.copy()
    np.fill_diagonal(sim_no_diag, 0)
    
    # Use geometric mean as adaptive threshold
    geometric_mean = np.exp(np.mean(np.log(sim_no_diag[sim_no_diag > 0])))
    network = (similarity_matrix >= geometric_mean).astype(float)
    np.fill_diagonal(network, 0)
    return network

def make_weighted_network(binary_network, similarity_matrix):
    """Convert binary network to weighted using original similarities"""
    weighted = binary_network * similarity_matrix
    np.fill_diagonal(weighted, 0)
    return weighted

def evaluate_clustering(true_labels, predicted_labels, similarity_matrix):
    """Compute all evaluation metrics"""
    # Remove any samples that might have been filtered out
    min_len = min(len(true_labels), len(predicted_labels))
    true_labels = true_labels[:min_len]
    predicted_labels = predicted_labels[:min_len]
    
    metrics = {}
    metrics['accuracy'] = adjusted_rand_score(true_labels, predicted_labels)  # ARI is better than accuracy for clustering
    metrics['silhouette'] = silhouette_score(similarity_matrix[:min_len, :min_len], predicted_labels, metric='precomputed')
    metrics['ari'] = adjusted_rand_score(true_labels, predicted_labels)
    metrics['nmi'] = normalized_mutual_info_score(true_labels, predicted_labels)
    
    return metrics

# Main analysis function
def analyze_data_type(data_name, similarity_matrix, true_labels):
    """Analyze one data type systematically"""
    print(f"\n{'='*50}")
    print(f"ANALYZING: {data_name}")
    print(f"{'='*50}")
    
    results = []
    
    # Test different network construction methods
    network_methods = {
        'KNN-5': lambda sim: create_knn_network(sim, 5),
        'KNN-10': lambda sim: create_knn_network(sim, 10),
        'KNN-15': lambda sim: create_knn_network(sim, 15),
        'KNN-20': lambda sim: create_knn_network(sim, 20),
        'Percentile-85': lambda sim: create_percentile_network(sim, 85),
        'Percentile-90': lambda sim: create_percentile_network(sim, 90),
        'Percentile-95': lambda sim: create_percentile_network(sim, 95),
        'Consensus': lambda sim: create_consensus_network(sim)
    }
    
    for method_name, method_func in network_methods.items():
        print(f"\nTesting {method_name}...")
        
        try:
            # Create binary network
            binary_network = method_func(similarity_matrix)
            
            # Test both binary and weighted versions
            for network_type in ['Binary', 'Weighted']:
                if network_type == 'Binary':
                    network = binary_network
                else:
                    network = make_weighted_network(binary_network, similarity_matrix)
                
                # Test different cluster numbers
                best_k = 3  # We know there are 3 true subtypes
                best_score = -1
                best_metrics = None
                
                for k in range(2, 8):  # Test k=2 to 7
                    try:
                        # Apply spectral clustering
                        spectral = SpectralClustering(n_clusters=k, affinity='precomputed', 
                                                    random_state=42)
                        predicted_labels = spectral.fit_predict(network)
                        
                        # Evaluate
                        metrics = evaluate_clustering(true_labels, predicted_labels, 1-similarity_matrix)
                        
                        # Use ARI as primary metric for best k selection
                        if metrics['ari'] > best_score:
                            best_score = metrics['ari']
                            best_k = k
                            best_metrics = metrics.copy()
                            
                    except Exception as e:
                        continue
                
                if best_metrics:
                    result = {
                        'data_type': data_name,
                        'method': method_name,
                        'network_type': network_type,
                        'best_k': best_k,
                        'accuracy': best_metrics['accuracy'],
                        'silhouette': best_metrics['silhouette'],
                        'ari': best_metrics['ari'],
                        'nmi': best_metrics['nmi']
                    }
                    results.append(result)
                    
                    print(f"  {network_type}: k={best_k}, ARI={best_metrics['ari']:.3f}, Sil={best_metrics['silhouette']:.3f}")
        
        except Exception as e:
            print(f"  Error with {method_name}: {str(e)[:50]}...")
            continue
    
    return results

# Run analysis for all data types
print(f"\nStarting systematic analysis...")
all_results = []

for data_name, matrix in similarity_matrices.items():
    results = analyze_data_type(data_name, matrix, true_labels)
    all_results.extend(results)

# Create results summary
if all_results:
    results_df = pd.DataFrame(all_results)
    results_df = results_df.sort_values(['data_type', 'ari'], ascending=[True, False])
    
    # Create clustering results folder and save results
    clustering_folder = "clustering results"
    os.makedirs(clustering_folder, exist_ok=True)
    output_path = os.path.join(clustering_folder, 'clustering_results_summary.csv')
    results_df.to_csv(output_path, index=False)
    
    print(f"\n{'='*70}")
    print("RESULTS SUMMARY")
    print(f"{'='*70}")
    
    # Show best result for each data type
    print(f"\nBEST RESULT FOR EACH DATA TYPE:")
    print("-" * 50)
    
    for data_type in results_df['data_type'].unique():
        best_result = results_df[results_df['data_type'] == data_type].iloc[0]
        print(f"{data_type}:")
        print(f"  Method: {best_result['method']} ({best_result['network_type']})")
        print(f"  Best k: {best_result['best_k']}")
        print(f"  ARI: {best_result['ari']:.3f}")
        print(f"  Silhouette: {best_result['silhouette']:.3f}")
        print(f"  NMI: {best_result['nmi']:.3f}")
        print()
    
    # Overall best
    overall_best = results_df.iloc[0]
    print(f"🏆 OVERALL BEST PERFORMER:")
    print(f"  Data Type: {overall_best['data_type']}")
    print(f"  Method: {overall_best['method']} ({overall_best['network_type']})")
    print(f"  ARI Score: {overall_best['ari']:.3f}")
    
    print(f"\n✓ Saved detailed results to: {output_path}")
else:
    print("❌ No results generated!")

print(f"\n{'='*70}")
print("PHASE 2 COMPLETE - NETWORK ANALYSIS DONE")
print("{'='*70}")