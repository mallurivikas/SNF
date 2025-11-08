# Cancer Subtype Discovery using Similarity Network Fusion (SNF)

## 🎯 Project Overview

This project implements a comprehensive cancer subtype discovery pipeline using Similarity Network Fusion (SNF) with multi-omics data. It analyzes gene expression, miRNA, and methylation data to identify cancer subtypes using network-based clustering approaches.

## 🧠 Methodological Approach

Our cancer subtype discovery follows a systematic, evidence-based approach with comprehensive method evaluation:

### Phase 1: Similarity Matrix Generation & Evaluation
**Objective**: Determine optimal similarity measurement approach
- **Multiple Similarity Metrics**: We computed similarity matrices using three different approaches:
  - **Cosine Similarity**: Measures angular similarity between sample vectors
  - **Euclidean Similarity**: Based on geometric distance in feature space  
  - **Correlation Similarity**: Pearson correlation between sample profiles
- **Comparative Analysis**: Evaluated all three methods across all data types
- **Key Decision**: Selected **Euclidean similarity** as optimal based on clustering performance and interpretability

### Phase 2: Multi-omics Data Fusion
**Objective**: Integrate information from different molecular data types
- **Individual Analysis**: First analyzed each data type independently:
  - Gene Expression (717 samples × 311 features)
  - miRNA Expression (717 samples × 74 features)  
  - Methylation (717 samples × 86 features)
- **Fusion Strategies**: Implemented multiple fusion approaches:
  - **Average Fusion**: Simple average of euclidean similarity matrices
  - **Weighted Fusion**: Weighted by feature count (Gene: 54%, miRNA: 13%, Methylation: 15%)
  - **Maximum Fusion**: Element-wise maximum across matrices
- **Recommendation**: Average fusion selected as primary approach for balanced integration

### Phase 3: Comprehensive Network Construction
**Objective**: Systematically test multiple network construction methods
- **K-Nearest Neighbor (KNN) Networks**: Tested with k = 5, 10, 15, 20
  - Creates connections to k most similar samples
  - Evaluated both binary and weighted versions
- **Percentile-based Networks**: Tested with 85th, 90th, 95th percentiles
  - Connects samples above similarity threshold
  - Adaptive to data distribution
- **Consensus Networks**: Geometric mean-based adaptive thresholding
  - Self-adjusting threshold based on data characteristics
- **Network Types**: Each method tested in both configurations:
  - **Binary Networks**: Unweighted connections (0/1)
  - **Weighted Networks**: Connection strength = original similarity value

### Phase 4: Spectral Clustering Optimization
**Objective**: Find optimal clustering parameters for each network
- **Cluster Number Testing**: Systematically tested k = 2 to 7 clusters
- **Method Combinations**: Total evaluation matrix:
  - 4 data types × 8 network methods × 2 network types × 6 k values = **384 combinations**
- **Performance Metrics**: Each combination evaluated using:
  - **ARI (Adjusted Rand Index)**: Agreement with ground truth cancer subtypes
  - **Silhouette Score**: Internal cluster quality
  - **NMI (Normalized Mutual Information)**: Information-theoretic similarity
- **Optimal Selection**: Best method chosen based on highest ARI score for each data type

### Phase 5: 3D Visualization & Interpretation
**Objective**: Create interpretable visualizations of clustering results
- **PCA Projection**: Reduced high-dimensional similarity space to 3D for visualization
- **Comparative Visualization**: Side-by-side comparison of:
  - Individual data type clustering results
  - Fused matrix clustering results  
  - Ground truth cancer subtypes
- **Performance Validation**: Visual confirmation of clustering quality against known cancer biology

## 📊 Key Features

- **Multi-omics Integration**: Combines gene expression, miRNA, and methylation data
- **Similarity Network Fusion**: Creates fused similarity matrices from different data types
- **Network-based Clustering**: Uses spectral clustering on various network construction methods
- **3D Visualization**: Generates interactive 3D cluster visualizations using PCA
- **Performance Evaluation**: Comprehensive evaluation using ARI, Silhouette Score, and NMI
- **Organized Output Structure**: All results saved in logically organized folders

## 🏗️ Project Structure

```
SNF/
├── Dataset/                                    # Input data files
│   ├── Gene_expression_final (1).xlsx        # Gene expression data (717×317)
│   ├── miRNA_target.xlsx                     # miRNA expression data (717×80)
│   └── TCGA_meth_717samples.xlsx            # Methylation data (717×92)
│
├── Actual SubType/                           # Ground truth cancer labels
│   └── real_cancer_subtype_labels.csv       # Real cancer subtype labels (HER2, LUM, TNBC)
│
├── similarity_matrices/                      # Similarity matrices and visualizations
│   ├── images/                              # Similarity matrix visualization plots
│   │   ├── all_similarity_matrices_overview.png
│   │   ├── Gene_expression_final_1_detailed_similarity_matrices.png
│   │   ├── miRNA_target_detailed_similarity_matrices.png
│   │   └── TCGA_meth_717samples_detailed_similarity_matrices.png
│   ├── Gene_expression_final (1)_cosine_similarity.csv
│   ├── Gene_expression_final (1)_euclidean_similarity.csv
│   ├── Gene_expression_final (1)_correlation_similarity.csv
│   ├── miRNA_target_cosine_similarity.csv
│   ├── miRNA_target_euclidean_similarity.csv
│   ├── miRNA_target_correlation_similarity.csv
│   ├── TCGA_meth_717samples_cosine_similarity.csv
│   ├── TCGA_meth_717samples_euclidean_similarity.csv
│   └── TCGA_meth_717samples_correlation_similarity.csv
│
├── fused_matrices/                          # Fused similarity matrices
│   ├── euclidean_fused_average.csv         # Average fusion (recommended)
│   ├── euclidean_fused_weighted.csv        # Weighted fusion by feature count
│   └── euclidean_fused_maximum.csv         # Maximum fusion
│
├── clustering results/                      # Clustering analysis results
│   └── clustering_results_summary.csv      # Performance metrics for all methods
│
├── cluster visualizations/                 # 3D cluster visualizations
│   ├── 3D_cluster_visualizations.png      # Overview of all clustering methods
│   ├── 3D_Gene_Expression_clusters.png    # Individual detailed visualizations
│   ├── 3D_miRNA_clusters.png
│   ├── 3D_Methylation_clusters.png
│   ├── 3D_Fused_Matrix_clusters.png
│   └── 3D_True_Cancer_Subtypes_clusters.png
│
├── Main Scripts/                           # Core analysis pipeline
│   ├── extract_labels.py                  # Step 1: Extract cancer subtype labels
│   ├── similarity_matrix.py               # Step 2: Generate similarity matrices
│   ├── fusion.py                          # Step 3: Network clustering analysis
│   ├── visualizations.py                 # Step 4: Create 3D visualizations
│   └── fused_matrix.py                   # Additional fusion utilities
│
├── Additional Files/                       # Supporting files
│   ├── SNF_similarity_matrices_comparison.png
│   └── SNF_similarity_matrices_detailed.png
│
└── README.md                              # This file
```

## 🚀 Installation & Setup

### Prerequisites

- Python 3.8 or higher
- Git (for cloning the repository)

### Installation Steps

1. **Clone the repository**:
```bash
git clone https://github.com/mallurivikas/SNF
cd SNF
```

2. **Create and activate virtual environment** (recommended):
```bash
# Create virtual environment
python -m venv venv

# Activate virtual environment
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate
```

3. **Install required packages**:
```bash
pip install -r requirements.txt
```

### Required Dependencies

The project uses the following core packages (automatically installed with requirements.txt):

- **Core Data Science**: `numpy`, `pandas`, `scikit-learn`, `scipy`
- **SNF Implementation**: `snfpy` 
- **Visualization**: `matplotlib`, `seaborn`, `plotly`
- **Network Analysis**: `networkx`, `igraph`
- **File Handling**: `openpyxl`, `xlrd`
- **Performance**: `joblib`, `threadpoolctl`

### Data Requirements

Place the following Excel files in the `Dataset/` folder:
- `Gene_expression_final (1).xlsx` - Gene expression data
- `miRNA_target.xlsx` - miRNA expression data  
- `TCGA_meth_717samples.xlsx` - Methylation data

## 📋 Execution Pipeline

Run the following scripts in order to complete the full analysis:

### Step 1: Extract Cancer Subtype Labels
```bash
python extract_labels.py
```
**Purpose**: Extracts real cancer subtype labels (HER2, LUM, TNBC) from the dataset  
**Output**: `Actual SubType/real_cancer_subtype_labels.csv`  
**Duration**: ~10 seconds

### Step 2: Generate Similarity Matrices
```bash
python similarity_matrix.py
```
**Purpose**: Computes similarity matrices using cosine, euclidean, and correlation methods  
**Outputs**: 
- `similarity_matrices/*.csv` - Individual similarity matrices
- `similarity_matrices/images/*.png` - Similarity visualizations  
**Duration**: ~2-3 minutes

### Step 3: Network Clustering Analysis
```bash
python fusion.py
```
**Purpose**: Performs network construction and spectral clustering using various methods  
**Methods Tested**: KNN (k=5,10,15,20), Percentile (85,90,95), Consensus  
**Output**: `clustering results/clustering_results_summary.csv`  
**Duration**: ~15-20 minutes

### Step 4: Create 3D Visualizations
```bash
python visualizations.py
```
**Purpose**: Generates 3D cluster visualizations using PCA projection  
**Outputs**: `cluster visualizations/*.png` - 3D visualization plots  
**Duration**: ~2-3 minutes

## 📊 Expected Results

Based on the analysis, you should expect:

| Data Type | Best Method | ARI Score | Silhouette Score |
|-----------|-------------|-----------|------------------|
| Gene Expression | KNN-10 (Binary) | ~0.774 | ~0.236 |
| miRNA | KNN-20 (Weighted) | ~0.772 | ~0.223 |
| Methylation | Consensus (Binary) | ~0.744 | ~0.364 |
| Fused Matrix | KNN-10 (Weighted) | ~0.741 | ~0.298 |

## 🔬 Technical Details

### Similarity Methods
- **Cosine Similarity**: Measures angular similarity between samples
- **Euclidean Similarity**: Based on geometric distance in feature space
- **Correlation Similarity**: Pearson correlation between sample profiles

### Network Construction
- **KNN Networks**: K-nearest neighbor graphs with varying k values
- **Percentile Thresholding**: Networks based on similarity percentile cutoffs
- **Consensus Networks**: Adaptive thresholding using geometric mean

### Fusion Approaches
- **Average Fusion**: Simple average of similarity matrices
- **Weighted Fusion**: Weighted by number of features in each data type
- **Maximum Fusion**: Element-wise maximum across matrices

### Clustering Evaluation
- **ARI (Adjusted Rand Index)**: Agreement with ground truth cancer subtypes
- **Silhouette Score**: Internal cluster quality measure
- **NMI (Normalized Mutual Information)**: Information-theoretic cluster similarity

## 📈 Performance Interpretation

- **ARI Score > 0.7**: Excellent agreement with true cancer subtypes
- **ARI Score 0.5-0.7**: Good clustering performance
- **ARI Score < 0.5**: Poor clustering performance

## 🛠️ Customization

### Adding New Data Types
1. Place new Excel file in `Dataset/` folder
2. Update file paths in relevant scripts
3. Modify similarity matrix computation as needed

### Trying Different Methods
- Modify network construction parameters in `fusion.py`
- Add new similarity metrics in `similarity_matrix.py`
- Experiment with different clustering algorithms

### Visualization Options
- Adjust PCA components in `visualizations.py`
- Try t-SNE instead of PCA for non-linear dimensionality reduction
- Modify color schemes and plot styling

## 🐛 Troubleshooting

### Common Issues

1. **File Not Found Errors**
   - Ensure all Excel files are in `Dataset/` folder
   - Check file names match exactly (including spaces and parentheses)

2. **Memory Issues**
   - Large similarity matrices may require more RAM
   - Consider reducing sample size for testing

3. **Slow Performance**
   - Spectral clustering can be computationally intensive
   - Consider reducing the range of k values tested

### Dependencies Issues
```bash
# If you encounter import errors, try reinstalling:
pip install -r requirements.txt

# Or install packages individually:
pip install numpy pandas scikit-learn scipy snfpy matplotlib seaborn plotly networkx igraph openpyxl xlrd joblib threadpoolctl
```

### Virtual Environment Setup
If you need to recreate the environment:
```bash
# Deactivate current environment
deactivate

# Remove old environment
rmdir /s venv  # Windows
# rm -rf venv  # macOS/Linux

# Create fresh environment
python -m venv venv
venv\Scripts\activate  # Windows
# source venv/bin/activate  # macOS/Linux

# Install dependencies
pip install -r requirements.txt
```

## 📚 References

- **Similarity Network Fusion**: Wang et al. (2014) "Similarity network fusion for aggregating data types on a genomic scale"
- **Spectral Clustering**: von Luxburg (2007) "A tutorial on spectral clustering"
- **Cancer Subtyping**: Curtis et al. (2012) "The genomic and transcriptomic architecture of 2,000 breast tumours"
