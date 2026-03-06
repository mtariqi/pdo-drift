<div align="center">

<!-- ═══════════════════════════════════════════════════════════ -->
<!--                     ANIMATED BANNER                        -->
<!-- ═══════════════════════════════════════════════════════════ -->

<img width="100%" src="https://capsule-render.vercel.app/api?type=waving&color=0:1B3A6B,50:2563A8,100:0F6B62&height=220&section=header&text=PDO%20Drift&fontSize=72&fontColor=FFFFFF&fontAlignY=40&desc=Multi-Omic%20Temporal%20Drift%20Characterization%20of%20Patient-Derived%20Organoids&descAlignY=62&descSize=17&animation=fadeIn" />

<!-- ═══════════════════════════════════════════════════════════ -->
<!--                    PROJECT TITLE                           -->
<!-- ═══════════════════════════════════════════════════════════ -->

<h1>
  <img src="https://readme-typing-svg.demolab.com?font=Fira+Code&weight=700&size=28&pause=1000&color=2563A8&center=true&vCenter=true&width=800&lines=PDO+Temporal+Drift+Characterization;Genomic+%E2%80%A2+Transcriptomic+%E2%80%A2+Functional;Multi-Omic+Longitudinal+Framework;github.com%2Fmtariqi%2Fpdo-drift" alt="Typing SVG" />
</h1>

<br/>

<!-- ═══════════════════════════════════════════════════════════ -->
<!--                  STATUS BADGES — ROW 1                     -->
<!-- ═══════════════════════════════════════════════════════════ -->

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.32-00AA88?style=for-the-badge&logo=snakemake&logoColor=white)](https://snakemake.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-F59E0B?style=for-the-badge&logo=opensourceinitiative&logoColor=white)](LICENSE)
[![CI Status](https://img.shields.io/github/actions/workflow/status/mtariqi/pdo-drift/ci.yml?style=for-the-badge&logo=githubactions&logoColor=white&label=CI)](https://github.com/mtariqi/pdo-drift/actions)

<br/>

<!-- ═══════════════════════════════════════════════════════════ -->
<!--                  STATUS BADGES — ROW 2                     -->
<!-- ═══════════════════════════════════════════════════════════ -->

![Code style: black](https://img.shields.io/badge/code%20style-black-000000?style=flat-square&logo=python)
![Coverage](https://img.shields.io/badge/coverage-≥85%25-2EA043?style=flat-square&logo=pytest)
![GATK4](https://img.shields.io/badge/GATK-4.x-E03C31?style=flat-square&logo=broadinstitute&logoColor=white)
![STAR](https://img.shields.io/badge/STAR-2.7.11a-1B3A6B?style=flat-square)
![DESeq2](https://img.shields.io/badge/DESeq2-PyDESeq2-8B5CF6?style=flat-square)
![Scanpy](https://img.shields.io/badge/scRNA--seq-Scanpy-F59E0B?style=flat-square)
![Bismark](https://img.shields.io/badge/RRBS-Bismark-0F6B62?style=flat-square)
![Conda](https://img.shields.io/badge/Conda-environment-44A833?style=flat-square&logo=anaconda&logoColor=white)

<br/>

<!-- ═══════════════════════════════════════════════════════════ -->
<!--                  PROJECT META BADGES                       -->
<!-- ═══════════════════════════════════════════════════════════ -->

![GitHub repo size](https://img.shields.io/github/repo-size/mtariqi/pdo-drift?style=flat-square&logo=github&color=1B3A6B)
![GitHub last commit](https://img.shields.io/github/last-commit/mtariqi/pdo-drift?style=flat-square&logo=git&color=2563A8)
![GitHub issues](https://img.shields.io/github/issues/mtariqi/pdo-drift?style=flat-square&logo=github&color=E03C31)
![GitHub stars](https://img.shields.io/github/stars/mtariqi/pdo-drift?style=flat-square&logo=github&color=F59E0B)
![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen?style=flat-square&logo=github)
![Maintained](https://img.shields.io/badge/maintained-yes-2EA043?style=flat-square)

<br/><br/>

<!-- ═══════════════════════════════════════════════════════════ -->
<!--                   ONE-LINE DESCRIPTION                     -->
<!-- ═══════════════════════════════════════════════════════════ -->

> **A production-grade, Snakemake-powered multi-omic framework for quantifying and characterizing**  
> **the genomic, transcriptomic, epigenomic, and functional drift of PDOs across serial passages**  
> **and freeze-thaw cycles — centred on the novel Multi-Omic Drift Score (MODS).**

<br/>

---

</div>

<!-- ═══════════════════════════════════════════════════════════ -->
<!--                   TABLE OF CONTENTS                        -->
<!-- ═══════════════════════════════════════════════════════════ -->

## 📋 Table of Contents

| # | Section |
|---|---------|
| 1 | [🧬 Scientific Background](#-scientific-background) |
| 2 | [🎯 Project Aims](#-project-aims) |
| 3 | [🏗️ Repository Structure](#️-repository-structure) |
| 4 | [⚙️ Installation](#️-installation) |
| 5 | [🔬 Pipelines Overview](#-pipelines-overview) |
| 6 | [🐍 Python Modules](#-python-modules) |
| 7 | [📊 The MODS Framework](#-the-mods-framework) |
| 8 | [🚀 Quick Start](#-quick-start) |
| 9 | [📈 Expected Outputs](#-expected-outputs) |
| 10 | [🧪 Testing](#-testing) |
| 11 | [📁 Data Management](#-data-management) |
| 12 | [🤝 Contributing](#-contributing) |
| 13 | [📜 Citation](#-citation) |
| 14 | [📄 License](#-license) |

---

## 🧬 Scientific Background

**Patient-Derived Organoids (PDOs)** are 3D self-organizing tissue cultures that faithfully recapitulate the architecture, molecular identity, and drug-response profiles of their tumour-of-origin. They have become the gold standard pre-clinical model in precision oncology — yet their very biological richness makes them susceptible to **temporal drift**: progressive, culture-induced molecular divergence from the source tissue over serial passages.

<div align="center">

```
Primary Tumour  ──►  PDO Establishment (P0)  ──►  Serial Passaging
                                                        │
                          Genomic Drift ◄───────────────┤
                      Transcriptomic Drift ◄────────────┤
                        Epigenomic Drift ◄──────────────┤
                         Functional Drift ◄─────────────┘
                                │
                                ▼
                     Multi-Omic Drift Score (MODS)
                    MODS = 0.30·GD + 0.35·TD + 0.15·ED + 0.20·FD
```

</div>

This project is the **first systematic, longitudinal, multi-omic characterization of PDO temporal drift** across a diverse panel of 10–15 clinically annotated lines, spanning colorectal, pancreatic, lung, breast, and gastric carcinomas.

---

## 🎯 Project Aims

<div align="center">

| Aim | Focus | Key Question |
|:---:|-------|-------------|
| 🧩 **Aim 1** | Genomic & Epigenomic Drift | How do PDOs evolve at the DNA and methylome level across P2→P20? |
| 📡 **Aim 2** | Transcriptomic Drift & Growth Kinetics | Why do some lines grow faster, and what transcriptional programs drive it? |
| ❄️ **Aim 3** | Freeze-Thaw Impact | Does cryopreservation introduce persistent molecular or functional perturbations? |

</div>

---

## 🏗️ Repository Structure

```
pdo-drift/
│
├── 📂 .github/
│   └── workflows/
│       ├── ci.yml                  # Lint + tests + Snakemake dry-run
│       └── snakemake.yml           # Pipeline validation on push
│
├── 📂 config/
│   ├── config.yaml                 # 🔧 Master project configuration
│   ├── samples.tsv                 # 📋 Sample manifest (edit this first!)
│   └── resources.yaml              # 💻 HPC/SLURM resource settings
│
├── 📂 envs/
│   ├── main.yaml                   # 🐍 Main conda environment (Python + bioinformatics)
│   └── r_analysis.yaml             # 📊 R environment (DESeq2, Seurat, MethylKit)
│
├── 📂 pipelines/
│   ├── wes/         Snakefile      # 🔬 WES: fastp→bwa-mem2→GATK4→CNVKit→PyClone-VI
│   ├── rnaseq/      Snakefile      # 📡 RNA-seq: fastp→STAR→featureCounts→DESeq2→GSEA
│   ├── rrbs/        Snakefile      # 🧬 RRBS: TrimGalore→Bismark→methylation extraction
│   └── scrna/       Snakefile      # 🔭 scRNA-seq: CellRanger→Scanpy→Monocle3
│
├── 📂 src/
│   ├── integration/
│   │   └── mods_calculator.py      # ⭐ MODS framework (GD + TD + ED + FD)
│   ├── transcriptomics/
│   │   └── drift_tracker.py        # PyDESeq2 · GSEA · PCA · LASSO growth signature
│   ├── functional/
│   │   └── drug_sensitivity.py     # 4PL fitting · IC50/AUC · freeze-thaw recovery
│   ├── genomics/                   # Clonal evolution · CNV helpers
│   ├── epigenomics/                # Beta-matrix builders · DMR analysis
│   ├── qc/
│   │   └── qc_reporter.py          # Aggregate QC · pass/fail flagging
│   └── utils/                      # Shared helpers
│
├── 📂 notebooks/
│   ├── 01_qc_overview.ipynb
│   ├── 02_genomic_drift.ipynb
│   ├── 03_transcriptomic_drift.ipynb
│   ├── 04_epigenomic_drift.ipynb
│   ├── 05_functional_drift.ipynb
│   ├── 06_mods_integration.ipynb   # ⭐ End-to-end MODS walkthrough
│   └── 07_freeze_thaw_analysis.ipynb
│
├── 📂 tests/
│   └── test_mods.py                # pytest · 25+ unit tests for all modules
│
├── 📂 results/                     # Auto-generated outputs (gitignored)
├── 📂 docs/                        # Extended documentation
│
├── setup.py
├── requirements.txt
├── .gitignore
└── README.md
```

---

## ⚙️ Installation

### Prerequisites

| Requirement | Version | Purpose |
|-------------|---------|---------|
| Linux / macOS (WSL2) | — | Supported OS |
| [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) or Mamba | ≥23.x | Environment management |
| Git | ≥2.30 | Version control |
| Disk space | ≥50 GB | Raw data + indices |

### Step 1 — Clone

```bash
git clone git@github.com:mtariqi/pdo-drift.git
cd pdo-drift
```

### Step 2 — Install Mamba *(faster solver)*

```bash
conda install -n base -c conda-forge mamba -y
```

### Step 3 — Create Environments

```bash
# Main Python + bioinformatics environment
mamba env create -f envs/main.yaml
conda activate pdo-drift

# Install the local package in editable mode
pip install -e .
```

### Step 4 — Verify Installation

```bash
python -c "from pdo_drift.integration import MODSCalculator; print('✅ MODS ready')"
snakemake --version   # should print ≥7.32
pytest tests/ -q      # all tests should pass
```

---

## 🔬 Pipelines Overview

<div align="center">

### Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        RAW INPUT DATA                           │
│              (FASTQ pairs per sample × passage)                 │
└───────────┬──────────────┬──────────────┬───────────────────────┘
            │              │              │              │
            ▼              ▼              ▼              ▼
     ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────────┐
     │  WES     │  │ RNA-seq  │  │  RRBS    │  │  scRNA-seq   │
     │ Pipeline │  │ Pipeline │  │ Pipeline │  │  Pipeline    │
     └────┬─────┘  └────┬─────┘  └────┬─────┘  └──────┬───────┘
          │             │              │               │
          ▼             ▼              ▼               ▼
       VCF+CNV      Counts+DGE    Beta Matrix    AnnData h5ad
          │             │              │               │
          └─────────────┴──────────────┴───────────────┘
                                │
                                ▼
                    ┌───────────────────────┐
                    │   MODS Calculator     │
                    │  GD · TD · ED · FD    │
                    └───────────┬───────────┘
                                │
                                ▼
                    ┌───────────────────────┐
                    │   Drift Report &      │
                    │   Visualizations      │
                    └───────────────────────┘
```

</div>

### 🔬 WES Pipeline — Somatic Variant + CNV Calling

```
fastp  →  bwa-mem2  →  MarkDuplicates  →  BQSR  →  MuTect2  →  FilterMutectCalls
       →  CNVKit  →  PyClone-VI  →  MultiQC
```

```bash
cd pipelines/wes
snakemake --configfile ../../config/config.yaml \
          --use-conda --cores 32 -p all
```

### 📡 RNA-seq Pipeline — Transcriptomic Drift

```
fastp  →  STAR  →  featureCounts  →  Salmon  →  DESeq2  →  GSEA  →  MultiQC
```

```bash
cd pipelines/rnaseq
snakemake --configfile ../../config/config.yaml \
          --use-conda --cores 32 -p all
```

### 🧬 RRBS Pipeline — Epigenomic Drift

```
Trim Galore (RRBS-aware)  →  Bismark  →  methylation extractor  →  CpG beta matrix
```

```bash
cd pipelines/rrbs
snakemake --configfile ../../config/config.yaml \
          --use-conda --cores 16 -p all
```

### 🔭 scRNA-seq Pipeline — Single-Cell Heterogeneity

```
CellRanger count  →  Scanpy QC  →  Harmony integration  →  Monocle3 trajectory
```

```bash
cd pipelines/scrna
snakemake --configfile ../../config/config.yaml \
          --use-conda --cores 32 -p all
```

### 🖥️ HPC / SLURM — Run All Pipelines

```bash
snakemake --profile config/slurm_profile/ \
          --jobs 200 \
          --use-conda \
          --latency-wait 120
```

---

## 🐍 Python Modules

All modules are importable after `pip install -e .`:

### `MODSCalculator` — Multi-Omic Drift Score

```python
from pdo_drift.integration import MODSCalculator

calc = MODSCalculator(config_path="config/config.yaml")

# Compute MODS for all samples vs. P2 baseline
scores = calc.compute_all(
    vcf_dir="results/variants/",
    rnaseq_dir="results/deseq2/",
    rrbs_dir="results/methylation/",
    drug_dir="results/functional/"
)

# Visualize
calc.plot_heatmap(scores, outfile="results/figures/mods_heatmap.pdf")
calc.plot_passage_trajectory(scores, line="PDO001")

print(scores[["pdo_line","passage","MODS","drift_flag"]])
```

### `DriftTracker` — Transcriptomic Analysis

```python
from pdo_drift.transcriptomics import DriftTracker

tracker = DriftTracker(
    counts_csv="results/counts/merged_counts.csv",
    metadata_tsv="config/samples.tsv"
)

# Differential expression (PyDESeq2)
dge = tracker.run_deseq2(design_factor="passage", baseline_level="P2")

# GSEA across all contrasts
gsea = tracker.run_gsea(gene_sets=["MSigDB_Hallmark_2023", "KEGG_2021_Human"])

# PCA drift plot
tracker.plot_pca_drift(outfile="results/figures/pca_drift.pdf")

# Growth rate predictive signature (LASSO)
sig = tracker.identify_growth_signature(growth_rates=doubling_times)
```

### `DrugSensitivityProfiler` — Functional Drift

```python
from pdo_drift.functional import DrugSensitivityProfiler

profiler = DrugSensitivityProfiler(config=config["drug_panel"])
profiler.load_plate_data("results/functional/raw_viability.csv")
profiler.normalize()
profiler.fit_all_curves()              # 4PL model per drug × sample

# Compute functional drift vs P2
drift = profiler.compute_drift(baseline_passage=2)

# Freeze-thaw recovery kinetics
recovery = profiler.freeze_thaw_recovery(timepoints=["24h","48h","72h","7d"])

# Visualize
profiler.plot_heatmap()
profiler.plot_dose_response_grid("PDO001", passages=[2, 5, 10, 15, 20])
```

### `QCReporter` — Quality Control

```python
from pdo_drift.qc import QCReporter

qc = QCReporter("config/samples.tsv")
qc.load_fastp_reports("results/qc/fastp/")
qc.load_star_logs("data/star/")
qc.load_picard_metrics("results/qc/picard/")

report = qc.generate_report(outfile="results/qc/qc_report.csv")
qc.plot_qc_summary()

# Check how many passed
print(f"QC pass rate: {report['qc_pass'].mean():.1%}")
```

---

## 📊 The MODS Framework

The **Multi-Omic Drift Score (MODS)** is the central contribution of this project — a composite, quantitative metric integrating four molecular distance dimensions:

<div align="center">

```
╔══════════════════════════════════════════════════════════════╗
║                                                              ║
║   MODS  =  0.30 × GD  +  0.35 × TD  +  0.15 × ED  +  0.20 × FD    ║
║                                                              ║
╠══════════════════════════════════════════════════════════════╣
║  GD  Genomic Distance       Cosine(mutation) + CNV delta     ║
║  TD  Transcriptomic Dist.   1 − Pearson(rlog, top-5k genes)  ║
║  ED  Epigenomic Distance    Mean |Δβ| at variable CpGs       ║
║  FD  Functional Distance    Normalised ΔAUC (48-drug panel)  ║
╠══════════════════════════════════════════════════════════════╣
║  MODS ∈ [0, 1]   │   Drift flag threshold: MODS > 0.25       ║
╚══════════════════════════════════════════════════════════════╝
```

</div>

### MODS Score Interpretation

| MODS Range | Status | Recommended Action |
|:----------:|--------|-------------------|
| `0.00 – 0.15` | 🟢 **High Fidelity** | Safe for all assays; certified for pharmacological studies |
| `0.15 – 0.25` | 🟡 **Conditional** | Acceptable with documentation; re-bank recommended |
| `0.25 – 0.40` | 🟠 **Drifted** | Flag for review; restrict to non-translational studies |
| `> 0.40` | 🔴 **Compromised** | Restore from cryobank or re-derive from primary tissue |

---

## 🚀 Quick Start

### 1 — Configure your samples

Edit `config/samples.tsv` to add your PDO lines and FASTQ paths:

```tsv
sample_id     pdo_line  tumor_type  passage  condition  fastq_r1                    fastq_r2
PDO001_P2     PDO001    CRC         2        fresh       data/raw/PDO001_P2_R1.fq.gz  data/raw/PDO001_P2_R2.fq.gz
PDO001_P5     PDO001    CRC         5        fresh       data/raw/PDO001_P5_R1.fq.gz  data/raw/PDO001_P5_R2.fq.gz
```

### 2 — Update `config/config.yaml`

```yaml
reference:
  genome: "/path/to/hg38/GRCh38.primary_assembly.genome.fa"
  gtf: "/path/to/gencode.v44.annotation.gtf"
```

### 3 — Run

```bash
conda activate pdo-drift

# Single pipeline
cd pipelines/rnaseq && snakemake --use-conda --cores 32 -p all

# Full project (SLURM)
snakemake --profile config/slurm_profile/ --jobs 200
```

### 4 — Compute MODS

```bash
python -c "
from pdo_drift.integration import MODSCalculator
calc = MODSCalculator('config/config.yaml')
scores = calc.compute_all('results/variants/', 'results/deseq2/',
                          'results/methylation/', 'results/functional/')
calc.plot_heatmap(scores)
scores.to_csv('results/mods_scores.csv', index=False)
print(scores[['pdo_line','passage','MODS','drift_flag']].to_string())
"
```

### 5 — Explore notebooks

```bash
jupyter lab notebooks/
# Start with: 06_mods_integration.ipynb
```

---

## 📈 Expected Outputs

<div align="center">

| Output | Location | Description |
|--------|----------|-------------|
| `mods_scores.csv` | `results/` | MODS scores for all samples |
| `mods_heatmap.pdf` | `results/figures/` | PDO lines × passages heatmap |
| `pca_drift.pdf` | `results/figures/` | Transcriptomic drift PCA |
| `drug_heatmap.pdf` | `results/figures/` | IC50 fold-change heatmap |
| `drift_dge_results.csv` | `results/deseq2/` | Differential gene expression |
| `hallmark_gsea_results.csv` | `results/gsea/` | Pathway enrichment per contrast |
| `curve_parameters.csv` | `results/functional/` | IC50 / AUC / Emax per drug |
| `freeze_thaw_recovery.csv` | `results/functional/` | Post-thaw recovery kinetics |
| `cpg_beta_matrix.csv` | `results/methylation/` | CpG beta-value matrix |
| `qc_report.csv` | `results/qc/` | Per-sample QC pass/fail |
| `integrated_anndata.h5ad` | `results/scrna/` | Scanpy AnnData object |
| `multiqc_*.html` | `results/qc/` | Aggregated QC reports |

</div>

---

## 🧪 Testing

```bash
# Run all unit tests
pytest tests/ -v --cov=src/ --cov-report=html

# Run specific test class
pytest tests/test_mods.py::TestTranscriptomicDistanceCalculator -v

# View HTML coverage report
open htmlcov/index.html
```

### Test Coverage

| Module | Tests | Coverage Target |
|--------|-------|----------------|
| `mods_calculator.py` | 12 tests | ≥90% |
| `drift_tracker.py` | 8 tests | ≥85% |
| `drug_sensitivity.py` | 7 tests | ≥85% |
| `qc_reporter.py` | 5 tests | ≥80% |

---

## 📁 Data Management

> ⚠️ **Raw sequencing data (FASTQ, BAM, VCF) are gitignored and must NEVER be committed to Git.**

| Data Type | Storage Location | Notes |
|-----------|-----------------|-------|
| Raw FASTQ | HPC scratch / AWS S3 | Update paths in `config/samples.tsv` |
| Reference genome | `data/reference/` (gitignored) | Download via `scripts/download_references.sh` |
| Processed results | `results/` (gitignored) | Re-generated from pipelines |
| Final figures | `results/figures/` | Include in publications |
| Public deposit | dbGaP + GEO | After peer review |

### Sample Manifest Format

```tsv
sample_id    patient_id  pdo_line  tumor_type  passage  condition       assay   batch
PDO001_P2    PT001       PDO001    CRC         2        fresh           rnaseq  B1
PDO001_FT24h PT001       PDO001    CRC         5        postthaw_24h    rnaseq  B1
```

**Condition codes:** `fresh` · `postthaw_24h` · `postthaw_48h` · `postthaw_72h` · `postthaw_7d`

---

## 🤝 Contributing

Contributions are warmly welcomed! Please follow these steps:

1. **Fork** the repository
2. **Create** a feature branch: `git checkout -b feature/your-feature`
3. **Commit** using [Conventional Commits](https://www.conventionalcommits.org/):
   ```
   feat(genomics): add FACETS allele-specific CNV support
   fix(mods): handle missing epigenomic component gracefully
   docs: update freeze-thaw recovery section
   ```
4. **Push** and open a **Pull Request** against `develop`
5. Ensure **CI passes** (lint + tests + dry-run) before requesting review

### Branch Strategy

```
main        ← stable, peer-reviewed code only
develop     ← integration branch (PRs target here)
feature/*   ← individual features / analyses
hotfix/*    ← urgent bug fixes
```

---

## 📜 Citation

If you use this framework or the PDO Temporal Drift Atlas in your research, please cite:

```bibtex
@software{pdo_drift_2026,
  author    = {Islam,Md Tariqul.},
  title     = {{PDO Temporal Drift Characterization Framework}:
               A Multi-Omic Longitudinal Analysis of Patient-Derived Organoid Drift},
  year      = {2026},
  version   = {1.0.0},
  publisher = {GitHub},
  url       = {https://github.com/mtariqi/pdo-drift},
  note      = {Project Code: PDO-DRIFT-001}
}
```

---

## 📄 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

<div align="center">

<!-- ═══════════════════════════════════════════════════════════ -->
<!--                    FOOTER WAVE                             -->
<!-- ═══════════════════════════════════════════════════════════ -->

<img width="100%" src="https://capsule-render.vercel.app/api?type=waving&color=0:0F6B62,50:2563A8,100:1B3A6B&height=120&section=footer&animation=fadeIn"/>

<br/>

**Built with ❤️ for the organoid research community**

[![GitHub](https://img.shields.io/badge/GitHub-mtariqi-181717?style=for-the-badge&logo=github)](https://github.com/mtariqi)
[![Repository](https://img.shields.io/badge/Repo-pdo--drift-2563A8?style=for-the-badge&logo=github)](https://github.com/mtariqi/pdo-drift)

<sub>PDO-DRIFT-001 · Version 1.0 · March 2026</sub>

</div>
