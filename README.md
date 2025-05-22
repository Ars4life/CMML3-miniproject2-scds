# CMML3-miniproject2-scds
This is the supplementary material of the CMML ICA2
# Single-Cell RNA-Seq Doublet Detection Benchmarking

This repository contains R scripts for benchmarking doublet detection methods in single-cell RNA sequencing (scRNA-seq) data. Below is a detailed description of each script and its purpose.


### 📜 Script Descriptions

#### 1. `02_pretreat.R`
- **Purpose**: Preprocess raw scRNA-seq data for downstream analysis.
- **Steps**:
  - Load raw data (`.rds` format) and extract the count matrix.
  - Create a Seurat object, perform normalization (`LogNormalize`), and scale data.
  - Identify highly variable features and run PCA (20 principal components).
  - Save processed data to `./02_data/03_processed_data/` and log processing metrics (time, cell/gene counts).

#### 2. `03_benchmark_real.R`
- **Purpose**: Evaluate doublet detection methods on real datasets (PBMC and HM-12k).
- **Methods**:
  - **scds**: Compute `cxds`, `bcds`, and `hybrid` doublet scores.
  - **Scrublet** (Python): Predict doublets using thresholds.
  - **DoubletFinder**: Optimize `pK` parameters for doublet detection.
- **Metrics**: Calculate precision, recall, and true negative rate (TNR) at 10%, 20%, and 40% detection rates.

#### 3. `04_benchmark_visualize.R`
- **Purpose**: Generate visualizations for benchmarking results.
- **Plots**:
  - **Grouped bar plots**: Compare precision across methods for PBMC and HM-12k datasets.
  - **Faceted line plots**: Show performance trends across simulated doublet rates.
  - **Composite bar plots**: Summarize precision, recall, and TNR from DE analysis.

#### 4. `05_ratesim.R`
- **Purpose**: Simulate datasets with varying doublet rates (2%, 4%, 8%, 10%, 20%, 40%) and evaluate method performance.
- **Workflow**:
  - Load pre-generated simulated data (`sim_rate.rds`).
  - Compute doublet scores using `scds`.
  - Calculate precision, recall, and TNR at different detection rates. Results are saved to `sim_rate2.csv`.

#### 5. `06_DE.R`
- **Purpose**: Assess the impact of doublet removal on differential expression (DE) analysis.
- **Steps**:
  - Load simulated data with ground-truth DE genes.
  - Detect doublets using `scds`, remove predicted doublets, and re-perform DE analysis (Wilcoxon test).
  - Calculate precision, recall, and TNR to quantify improvements in DE results.

---

### 🛠️ Dependencies
- **R Version**: ≥4.0
- **Key R Packages**: `Seurat`, `scds`, `SingleCellExperiment`, `PRROC`, `reticulate`, `ggplot2`, `dplyr`.
- **Python Modules**: `scrublet`, `scipy`, `numpy` (for Scrublet integration).

---

Adjust file paths as needed for your environment. To reproduce results, execute scripts in numerical order. Ensure all dependencies are installed beforehand.
The dataset mentioned in the project report can all be download from [avialable here](https://zenodo.org/records/4062232#.X6GordD0laQ)
