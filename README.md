# PDO Temporal Drift Characterization Project

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-green.svg)](https://snakemake.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

> **Quantifying genomic, transcriptomic, and functional drift of Patient-Derived Organoids (PDOs) across serial passages and freeze-thaw cycles.**

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Repository Structure](#repository-structure)
3. [Installation & Environment Setup](#installation--environment-setup)
4. [Data Organization](#data-organization)
5. [Running the Pipelines](#running-the-pipelines)
6. [Python Analysis Modules](#python-analysis-modules)
7. [Jupyter Notebooks](#jupyter-notebooks)
8. [Configuration](#configuration)
9. [Testing](#testing)
10. [Contributing](#contributing)

---

## Project Overview

This project systematically characterizes **temporal drift** in Patient-Derived Organoids (PDOs) using a multi-omic longitudinal framework:

| Aim | Focus | Key Tools |
|-----|-------|-----------|
| **Aim 1** | Genomic & epigenomic drift (WES + RRBS) | GATK4, PyClone-VI, MethylKit |
| **Aim 2** | Transcriptomic drift & growth kinetics (RNA-seq, scRNA-seq) | DESeq2, Seurat, Monocle3 |
| **Aim 3** | Freeze-thaw molecular & functional impact | Custom Python MODS framework |

### Multi-Omic Drift Score (MODS)
```
MODS = 0.30×GD + 0.35×TD + 0.15×ED + 0.20×FD
```
Where GD=Genomic Distance, TD=Transcriptomic Distance, ED=Epigenomic Distance, FD=Functional Distance

---

## Repository Structure

```
pdo-drift/
├── .github/
│   └── workflows/          # CI/CD GitHub Actions
│       ├── ci.yml           # Linting, tests
│       └── snakemake.yml    # Pipeline dry-run validation
├── config/
│   ├── config.yaml          # Master project config
│   ├── samples.tsv          # Sample manifest
│   └── resources.yaml       # HPC resource allocation
├── data/
│   ├── raw/                 # Raw FASTQ, BAM (gitignored)
│   ├── processed/           # Pipeline outputs (gitignored)
│   └── reference/           # Genome, annotation files
├── envs/
│   ├── main.yaml            # Main conda environment
│   ├── r_analysis.yaml      # R-based tools environment
│   └── deeptools.yaml       # DeepTools / ChIP-seq env
├── notebooks/
│   ├── 01_qc_overview.ipynb
│   ├── 02_genomic_drift.ipynb
│   ├── 03_transcriptomic_drift.ipynb
│   ├── 04_epigenomic_drift.ipynb
│   ├── 05_functional_drift.ipynb
│   ├── 06_mods_integration.ipynb
│   └── 07_freeze_thaw_analysis.ipynb
├── pipelines/
│   ├── wes/                 # WES Snakemake pipeline
│   ├── rnaseq/              # Bulk RNA-seq pipeline
│   ├── scrna/               # scRNA-seq pipeline
│   └── rrbs/                # RRBS methylation pipeline
├── src/
│   ├── qc/                  # QC modules
│   ├── genomics/            # WES analysis
│   ├── transcriptomics/     # RNA-seq analysis
│   ├── epigenomics/         # RRBS analysis
│   ├── functional/          # Drug sensitivity, growth kinetics
│   ├── integration/         # MODS score computation
│   └── utils/               # Shared utilities
├── tests/                   # Unit + integration tests
├── results/                 # Output figures, tables, reports
├── docs/                    # Extended documentation
├── setup.py
├── requirements.txt
└── README.md
```

---

## Installation & Environment Setup

### Prerequisites
- Linux / macOS (WSL2 on Windows)
- [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) or Mamba
- Git ≥ 2.30
- 50 GB disk space (raw data not included)

### Step 1 — Clone the Repository

```bash
git clone https://github.com/mtariqi/pdo-drift.git
cd pdo-drift
```

### Step 2 — Install Mamba (faster conda solver)

```bash
conda install -n base -c conda-forge mamba -y
```

### Step 3 — Create Environments

```bash
# Main Python environment
mamba env create -f envs/main.yaml
conda activate pdo-drift

# R analysis environment (DESeq2, Seurat)
mamba env create -f envs/r_analysis.yaml

# Install the local package
pip install -e .
```

### Step 4 — Configure Reference Genome

```bash
# Download hg38 reference (example using iGenomes)
bash scripts/download_references.sh hg38

# Or specify a custom path in config/config.yaml:
# reference:
#   genome: /path/to/hg38.fa
#   gtf: /path/to/gencode.v44.annotation.gtf
```

---

## Data Organization

### Sample Manifest (`config/samples.tsv`)

| sample_id | patient_id | tumor_type | passage | condition | batch | fastq_r1 | fastq_r2 |
|-----------|-----------|------------|---------|-----------|-------|----------|----------|
| PDO001_P2 | PT001 | CRC | 2 | fresh | B1 | data/raw/PDO001_P2_R1.fastq.gz | ... |
| PDO001_P5 | PT001 | CRC | 5 | fresh | B1 | ... | ... |

### Condition Codes
- `fresh` — never frozen, directly passaged
- `postthaw_24h` — 24 hours after thaw
- `postthaw_48h`, `postthaw_72h`, `postthaw_7d`

---

## Running the Pipelines

### WES Pipeline (Somatic Variant Calling)

```bash
cd pipelines/wes
snakemake --configfile ../../config/config.yaml \
          --use-conda \
          --cores 32 \
          --jobs 10 \
          -p all
```

### RNA-seq Pipeline

```bash
cd pipelines/rnaseq
snakemake --configfile ../../config/config.yaml \
          --use-conda --cores 32 -p all
```

### RRBS Methylation Pipeline

```bash
cd pipelines/rrbs
snakemake --configfile ../../config/config.yaml \
          --use-conda --cores 16 -p all
```

### scRNA-seq Pipeline

```bash
cd pipelines/scrna
snakemake --configfile ../../config/config.yaml \
          --use-conda --cores 32 -p all
```

### Run All Pipelines (HPC / SLURM)

```bash
snakemake --profile config/slurm_profile/ --jobs 200
```

---

## Python Analysis Modules

All modules are importable after `pip install -e .`:

```python
from pdo_drift.integration import MODSCalculator
from pdo_drift.genomics import ClonalEvolutionAnalyzer
from pdo_drift.transcriptomics import DriftTracker
from pdo_drift.functional import DrugSensitivityProfiler
```

See [docs/API.md](docs/API.md) for full API reference.

---

## Jupyter Notebooks

Run notebooks in order:

```bash
conda activate pdo-drift
jupyter lab notebooks/
```

Each notebook is self-contained with embedded documentation.

---

## Testing

```bash
pytest tests/ -v --cov=src/ --cov-report=html
```

---

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Commit with conventional commits: `git commit -m "feat: add MODS visualization"`
4. Push and open a Pull Request

See [CONTRIBUTING.md](docs/CONTRIBUTING.md) for full guidelines.

---

## Citation

If you use this framework, please cite:

```bibtex
@software{pdo_drift_2025,
  author = {Islam, Md Tariqi},
  title  = {PDO Temporal Drift Characterization Framework},
  year   = {2025},
  url    = {https://github.com/mtariqi/pdo-drift}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.
