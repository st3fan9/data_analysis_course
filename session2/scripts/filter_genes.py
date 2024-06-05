import argparse
from pathlib import Path

import scanpy as sc
import pandas as pd
import numpy as np

# Parse arguments
parser = argparse.ArgumentParser(description='Filter cells based on QC metrics')
parser.add_argument("--input_file", type=Path, required=True, help="Input file with good-quality cells")
parser.add_argument("--output_file", type=Path, required=True, help="Output file without uninformative genes")

args = parser.parse_args()

# Load data
adata = sc.read_h5ad(args.input_file)
print(f"Number of genes before filtering for sample {args.input_file.stem}: {adata.n_vars}")

# Filter mitochondrial genes
# We have annotated these under the "mito" column in adata.var
adata = adata[:, ~adata.var["mito"]].copy()

# Filter genes that are expressed in less than 3 cells
sc.pp.filter_genes(adata, min_cells=3)

# Print the number of genes after filtering
print(f"Number of genes after filtering for sample {args.input_file.stem}: {adata.n_vars}")

# Save the filtered data
adata.write(args.output_file)
