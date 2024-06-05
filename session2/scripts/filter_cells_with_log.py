# Import necessary libraries
import argparse
from pathlib import Path

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from utils.logging import configure_logging

# Parse arguments
parser = argparse.ArgumentParser(description='Filter cells based on QC metrics')
parser.add_argument("--input_file", type=Path, required=True, help="Input file")
parser.add_argument("--output_file", type=Path, required=True, help="Output file with good-quality cells")
parser.add_argument("--qc_violin_plot", type=Path, help="Path to save violin plot of QC metrics")
parser.add_argument("--log", type=Path, help="Path to save log file")

args = parser.parse_args()

# Configure logging
LOGGER.configure_logging(
    log_file=args.log,
    name=Path(args.log).stem
)

# Load input data
adata = sc.read_h5ad(args.input_file)

# Plot QC violin plots
# There is a scanpy violin plot function that does this
sc.pl.violin(adata, ["n_genes", "n_counts", "fraction_mito"], jitter=0.4, multi_panel=True, show=False, save=args.qc_violin_plot)

# Filter cells based on QC metrics
min_n_genes = 200
min_n_counts = 1000
max_fraction_mito = 0.05

# Label cells in unfiltered data based on whether they pass each threshold
adata.obs["pass_n_genes"] = adata.obs["n_genes"] >= min_n_genes
adata.obs["pass_n_counts"] = adata.obs["n_counts"] >= min_n_counts
adata.obs["pass_fraction_mito"] = adata.obs["fraction_mito"] <= max_fraction_mito
adata.obs["pass_all_qc"] = adata.obs["pass_n_genes"] & adata.obs["pass_n_counts"] & adata.obs["pass_fraction_mito"]

LOGGER.info(f"Filtering cells based on QC metrics: n_genes >= {min_n_genes}, n_counts >= {min_n_counts}, fraction_mito <= {max_fraction_mito}")
LOGGER.info(f"Number of cells before filtering: {adata.n_obs}")
LOGGER.info(f"Number of cells passing all QC thresholds: {np.sum(adata.obs['pass_all_qc'])}")
LOGGER.info(f"Fraction of cells passing all QC thresholds: {np.mean(adata.obs['pass_all_qc'])}")
LOGGER.info(f"Fraction of cells passing n_genes threshold: {np.mean(adata.obs['pass_n_genes'])}")
LOGGER.info(f"Fraction of cells passing n_counts threshold: {np.mean(adata.obs['pass_n_counts'])}")
LOGGER.info(f"Fraction of cells passing fraction_mito threshold: {np.mean(adata.obs['pass_fraction_mito'])}")

# Save filtered data
adata_filtered = adata[adata.obs["pass_all_qc"]].copy()
adata_filtered.write_h5ad(args.output_file)
