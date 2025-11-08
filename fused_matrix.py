import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set plotting style
plt.style.use('default')
sns.set_palette("viridis")

print("="*70)
print("SNF SIMILARITY MATRICES VISUALIZATION")
print("="*70)

# Load the three original Euclidean similarity matrices
matrices_dir = "similarity_matrices"
datasets = {
    'Gene Expression': 'Gene_expression_final (1)_euclidean_similarity.csv',
    'miRNA Target': 'miRNA_target_euclidean_similarity.csv',
    'Methylation': 'TCGA_meth_717samples_euclidean_similarity.csv'
}

# Load the fused matrix
fused_matrix_path = "fused_matrices/euclidean_fused_average.csv"

print("Loading similarity matrices...")

# Load original matrices
original_matrices = {}
for name, filename in datasets.items():
    filepath = os.path.join(matrices_dir, filename)
    if os.path.exists(filepath):
        matrix = pd.read_csv(filepath).values
        original_matrices[name] = matrix
        print(f"✓ Loaded {name}: {matrix.shape}")
    else:
        print(f"✗ File not found: {filename}")

# Load fused matrix
if os.path.exists(fused_matrix_path):
    fused_matrix = pd.read_csv(fused_matrix_path).values
    print(f"✓ Loaded Fused Matrix: {fused_matrix.shape}")
else:
    print("✗ Fused matrix not found!")

print(f"\nCreating 4-matrix comparison visualization...")

# Create the main visualization: 2x2 grid showing all 4 matrices
fig, axes = plt.subplots(2, 2, figsize=(16, 14))
fig.suptitle('SNF Similarity Matrices: Individual Data Types + Fused Matrix', 
             fontsize=18, fontweight='bold', y=0.95)

# Matrix names and positions
matrix_info = [
    ('Gene Expression\nEuclidean Similarity', original_matrices.get('Gene Expression')),
    ('miRNA Target\nEuclidean Similarity', original_matrices.get('miRNA Target')),
    ('Methylation\nEuclidean Similarity', original_matrices.get('Methylation')),
    ('Fused Matrix\n(Combined)', fused_matrix)
]

positions = [(0,0), (0,1), (1,0), (1,1)]

# Plot each matrix
for (title, matrix), (row, col) in zip(matrix_info, positions):
    if matrix is not None:
        # Show full matrix
        im = axes[row, col].imshow(matrix, cmap='viridis', aspect='auto', 
                                  vmin=0, vmax=1)
        axes[row, col].set_title(title, fontsize=14, fontweight='bold')
        axes[row, col].set_xlabel('Sample Index')
        axes[row, col].set_ylabel('Sample Index')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=axes[row, col], fraction=0.046, pad=0.04)
        cbar.set_label('Similarity', rotation=270, labelpad=15)
        
        # Add statistics as text
        mask = np.triu(np.ones_like(matrix, dtype=bool), k=1)
        off_diagonal = matrix[mask]
        mean_sim = np.mean(off_diagonal)
        std_sim = np.std(off_diagonal)
        
        stats_text = f'Mean: {mean_sim:.3f}\nStd: {std_sim:.3f}'
        axes[row, col].text(0.02, 0.98, stats_text, transform=axes[row, col].transAxes,
                           verticalalignment='top', bbox=dict(boxstyle='round', 
                           facecolor='white', alpha=0.8), fontsize=10)
    else:
        axes[row, col].text(0.5, 0.5, 'Matrix Not Available', 
                           ha='center', va='center', fontsize=14)
        axes[row, col].set_title(title, fontsize=14, fontweight='bold')

plt.tight_layout()
plt.subplots_adjust(top=0.91)
plt.savefig('SNF_similarity_matrices_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# Create a detailed view with sample regions (first 100x100)
fig, axes = plt.subplots(2, 2, figsize=(16, 14))
fig.suptitle('SNF Similarity Matrices: Detailed View (First 100×100 samples)', 
             fontsize=18, fontweight='bold', y=0.95)

sample_size = 100

for (title, matrix), (row, col) in zip(matrix_info, positions):
    if matrix is not None:
        # Show sample region
        sample_matrix = matrix[:sample_size, :sample_size]
        im = axes[row, col].imshow(sample_matrix, cmap='viridis', aspect='auto', 
                                  vmin=0, vmax=1)
        axes[row, col].set_title(title, fontsize=14, fontweight='bold')
        axes[row, col].set_xlabel('Sample Index')
        axes[row, col].set_ylabel('Sample Index')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=axes[row, col], fraction=0.046, pad=0.04)
        cbar.set_label('Similarity', rotation=270, labelpad=15)
        
        # Add statistics
        mean_sim = np.mean(sample_matrix)
        std_sim = np.std(sample_matrix)
        stats_text = f'Mean: {mean_sim:.3f}\nStd: {std_sim:.3f}'
        axes[row, col].text(0.02, 0.98, stats_text, transform=axes[row, col].transAxes,
                           verticalalignment='top', bbox=dict(boxstyle='round', 
                           facecolor='white', alpha=0.8), fontsize=10)
    else:
        axes[row, col].text(0.5, 0.5, 'Matrix Not Available', 
                           ha='center', va='center', fontsize=14)
        axes[row, col].set_title(title, fontsize=14, fontweight='bold')

plt.tight_layout()
plt.subplots_adjust(top=0.91)
plt.savefig('SNF_similarity_matrices_detailed.png', dpi=300, bbox_inches='tight')
plt.show()

# Print summary statistics
print(f"\n{'='*70}")
print("SIMILARITY MATRICES SUMMARY")
print(f"{'='*70}")

print(f"\n📊 MATRIX STATISTICS:")
for title, matrix in zip(['Gene Expression', 'miRNA Target', 'Methylation', 'Fused Matrix'], 
                        [original_matrices.get('Gene Expression'), 
                         original_matrices.get('miRNA Target'),
                         original_matrices.get('Methylation'), 
                         fused_matrix]):
    if matrix is not None:
        mask = np.triu(np.ones_like(matrix, dtype=bool), k=1)
        off_diagonal = matrix[mask]
        print(f"\n{title}:")
        print(f"  Shape: {matrix.shape}")
        print(f"  Mean similarity: {np.mean(off_diagonal):.4f}")
        print(f"  Std deviation: {np.std(off_diagonal):.4f}")
        print(f"  Min similarity: {np.min(off_diagonal):.4f}")
        print(f"  Max similarity: {np.max(off_diagonal):.4f}")

print(f"\n💾 VISUALIZATION FILES CREATED:")
print(f"  ✓ SNF_similarity_matrices_comparison.png (Full matrices)")
print(f"  ✓ SNF_similarity_matrices_detailed.png (Detailed 100×100 view)")

print(f"\n📁 ESSENTIAL FILES RETAINED:")
print(f"  ✓ similarity_matrices/ (Original Euclidean similarity matrices)")
print(f"  ✓ fused_matrices/ (Combined similarity matrix)")
print(f"  ✓ Dataset/ (Original Excel files)")

print(f"\n🗑️  CLEANED UP:")
print(f"  ✓ Removed unnecessary visualization files")
print(f"  ✓ Removed intermediate analysis scripts")
print(f"  ✓ Removed comparison plots folder")

print(f"\n{'='*70}")
print("SNF MATRICES READY FOR NEXT STEPS")
print(f"{'='*70}")
print("You now have the 4 essential similarity matrices visualized:")
print("1. Gene Expression Euclidean Similarity")
print("2. miRNA Target Euclidean Similarity") 
print("3. Methylation Euclidean Similarity")
print("4. Fused Matrix (Combination of all three)")
print("\nReady for your next instructions! 🚀")