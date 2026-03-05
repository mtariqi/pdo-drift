# 🚀 Step-by-Step GitHub Launch Guide
## PDO Drift Characterization Project → github.com/mtariqi/pdo-drift

---

## PHASE 1 — One-Time Setup (Do Once)

### Step 1: Install Git (if not already installed)

```bash
# Ubuntu / Debian
sudo apt update && sudo apt install git -y

# macOS (via Homebrew)
brew install git

# Verify
git --version   # should show git 2.30+
```

---

### Step 2: Configure Git with Your Identity

```bash
git config --global user.name "mtariqi"
git config --global user.email "your.email@institution.edu"
git config --global init.defaultBranch main
git config --global core.editor "nano"   # or vim, code, etc.

# Verify
git config --list
```

---

### Step 3: Create SSH Key for GitHub Authentication

```bash
# Generate SSH key (press Enter for all prompts)
ssh-keygen -t ed25519 -C "your.email@institution.edu"

# Start SSH agent
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519

# Copy public key to clipboard
cat ~/.ssh/id_ed25519.pub
# → Copy the entire output
```

**Then on GitHub:**
1. Go to → https://github.com/settings/keys
2. Click **New SSH key**
3. Title: `my-workstation` (or HPC cluster name)
4. Paste the public key → **Add SSH key**

```bash
# Test the connection
ssh -T git@github.com
# Expected: "Hi mtariqi! You've successfully authenticated..."
```

---

### Step 4: Create the GitHub Repository

**Option A — GitHub website:**
1. Go to https://github.com/new
2. Repository name: `pdo-drift`
3. Description: `Multi-omic characterization of PDO temporal drift`
4. Set to **Private** (recommended for unpublished data)
5. ✅ Do NOT initialize with README (we'll push our own)
6. Click **Create repository**

**Option B — GitHub CLI (faster):**
```bash
# Install GitHub CLI
sudo apt install gh -y   # or: brew install gh

# Login
gh auth login   # follow prompts

# Create repo directly
gh repo create mtariqi/pdo-drift \
    --private \
    --description "Multi-omic characterization of PDO temporal drift" \
    --clone=false
```

---

## PHASE 2 — Initialize and Push the Project

### Step 5: Navigate to the Project Directory

```bash
cd /path/to/pdo-drift    # wherever you saved the project files
# OR if starting fresh:
mkdir pdo-drift && cd pdo-drift
```

---

### Step 6: Initialize Git Repository

```bash
# Initialize local repo
git init

# Confirm the project files are present
ls -la
# Should show: README.md, src/, pipelines/, config/, envs/, tests/, etc.
```

---

### Step 7: Stage All Files

```bash
# Stage everything
git add .

# Verify what will be committed
git status
# Confirm .gitignore is working — you should NOT see:
#   data/raw/, data/aligned/, *.bam, *.fastq.gz, *.vcf
```

---

### Step 8: Create the Initial Commit

```bash
git commit -m "feat: initial project scaffold — PDO temporal drift characterization

- Multi-omic pipeline: WES (GATK4), RNA-seq (STAR+DESeq2), RRBS (Bismark), scRNA-seq (10x)
- MODS framework: genomic, transcriptomic, epigenomic, functional distance scoring
- Snakemake pipelines for all 4 assay types
- Python modules: MODSCalculator, DriftTracker, DrugSensitivityProfiler, QCReporter
- Unit tests for all core modules
- GitHub Actions CI: lint, pytest, Snakemake dry-run
- Conda environment definitions
- Full config system (config.yaml, samples.tsv)"
```

---

### Step 9: Connect to GitHub and Push

```bash
# Add GitHub remote
git remote add origin git@github.com:mtariqi/pdo-drift.git

# Rename branch to main (if needed)
git branch -M main

# Push to GitHub
git push -u origin main

# Verify — go to: https://github.com/mtariqi/pdo-drift
```

---

## PHASE 3 — Set Up Branching Strategy

### Step 10: Create Development Branch Structure

```bash
# Create and push a develop branch
git checkout -b develop
git push -u origin develop

# Create feature branches for each analysis aim
git checkout -b feature/aim1-genomic-drift
git push -u origin feature/aim1-genomic-drift

git checkout develop
git checkout -b feature/aim2-transcriptomic-drift
git push -u origin feature/aim2-transcriptomic-drift

git checkout develop
git checkout -b feature/aim3-freeze-thaw
git push -u origin feature/aim3-freeze-thaw

# Return to main
git checkout main
```

**Branch naming convention:**
```
main           → stable, peer-reviewed code only
develop        → integration branch
feature/aim1-* → Aim 1 work
feature/aim2-* → Aim 2 work
feature/aim3-* → Aim 3 work
hotfix/*       → urgent bug fixes
release/*      → release preparation
```

---

## PHASE 4 — Set Up Conda Environment Locally

### Step 11: Install Mamba and Create Environments

```bash
# Install Mamba (faster conda)
conda install -n base -c conda-forge mamba -y

# Create main analysis environment
mamba env create -f envs/main.yaml
conda activate pdo-drift

# Install local package in editable mode
pip install -e .

# Verify key tools
python -c "import scanpy, pydeseq2, gseapy; print('All imports OK')"
snakemake --version
```

---

## PHASE 5 — Configure Reference Data

### Step 12: Download Reference Genome

```bash
# Create reference directory
mkdir -p data/reference/hg38

# Download hg38 (GENCODE primary assembly)
wget -P data/reference/hg38 \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz

# Decompress
gunzip data/reference/hg38/GRCh38.primary_assembly.genome.fa.gz

# Download GENCODE v44 GTF
wget -P data/reference \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip data/reference/gencode.v44.annotation.gtf.gz

# Build BWA-MEM2 index (for WES)
bwa-mem2 index data/reference/hg38/GRCh38.primary_assembly.genome.fa

# Build STAR index (for RNA-seq, requires ~30GB RAM)
mkdir -p data/reference/STAR_index
STAR --runMode genomeGenerate \
    --genomeDir data/reference/STAR_index/ \
    --genomeFastaFiles data/reference/hg38/GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile data/reference/gencode.v44.annotation.gtf \
    --runThreadN 16

# Build Salmon index (for RNA-seq quantification)
salmon index \
    -t data/reference/gencode.v44.transcripts.fa.gz \
    -i data/reference/salmon_index \
    --gencode -p 16
```

---

## PHASE 6 — Run Tests and CI

### Step 13: Run Unit Tests Locally

```bash
conda activate pdo-drift

# Run all tests
pytest tests/ -v --cov=src/ --cov-report=html

# Open coverage report
open htmlcov/index.html   # macOS
xdg-open htmlcov/index.html   # Linux
```

### Step 14: Verify Snakemake Pipelines (Dry-Run)

```bash
# WES pipeline dry-run (no data needed)
cd pipelines/wes
snakemake --configfile ../../config/config.yaml --dryrun -p

# RNA-seq pipeline dry-run
cd ../rnaseq
snakemake --configfile ../../config/config.yaml --dryrun -p

cd ../..
```

---

## PHASE 7 — Configure GitHub Actions

### Step 15: Enable GitHub Actions

GitHub Actions runs automatically on every push. Check status at:
```
https://github.com/mtariqi/pdo-drift/actions
```

**First push will trigger:**
- ✅ Lint check (black + flake8)
- ✅ Unit tests (pytest)
- ✅ Snakemake dry-runs

If tests fail, fix locally and push again:
```bash
# Auto-format code with black
black src/ tests/

# Re-run tests
pytest tests/ -v

# Commit and push fix
git add .
git commit -m "fix: resolve test failures"
git push
```

---

## PHASE 8 — Configure Large File Storage (for BAMs / FASTQs)

### Step 16: Set Up Git LFS (for large reference files only)

```bash
# Install Git LFS
git lfs install

# Track large reference files (NOT raw data — use .gitignore for those)
git lfs track "data/reference/*.vcf.gz"
git lfs track "data/reference/*.fa"

# The .gitattributes file is auto-created — commit it
git add .gitattributes
git commit -m "chore: configure Git LFS for reference files"
git push
```

> ⚠️ **Raw FASTQ/BAM files should NEVER be in Git.** Store them on:
> - HPC scratch/project storage
> - AWS S3 / Google Cloud Storage
> - NCBI SRA (after publication)

---

## PHASE 9 — Daily Development Workflow

### Step 17: Standard Git Workflow for Daily Work

```bash
# 1. Start from develop, pull latest
git checkout develop
git pull origin develop

# 2. Work on a feature branch
git checkout feature/aim1-genomic-drift

# 3. Make changes, stage and commit frequently
git add src/genomics/clonal_evolution.py
git commit -m "feat(genomics): add PyClone-VI runner script

- Parses filtered MuTect2 VCF to PyClone input format
- Integrates CNVKit copy number calls
- Outputs clonal evolution TSV per sample"

# 4. Push feature branch
git push origin feature/aim1-genomic-drift

# 5. Open Pull Request on GitHub to merge into develop
# → https://github.com/mtariqi/pdo-drift/pulls → New pull request

# 6. After PR review + CI passes, merge into develop
# 7. When ready for release, merge develop → main
```

---

## PHASE 10 — Collaborate and Share

### Step 18: Add Collaborators

On GitHub:
1. Go to your repo → **Settings** → **Collaborators**
2. Add collaborators by GitHub username
3. Assign roles: Read / Write / Maintain / Admin

### Step 19: Create a Release

When Phase I is complete:
```bash
git checkout main
git pull

# Tag the release
git tag -a v1.0.0 -m "Phase I complete: baseline multi-omic dataset"
git push origin v1.0.0

# On GitHub: Releases → Draft a new release → Select tag v1.0.0
```

---

## Quick Reference — Most Used Commands

```bash
# Daily sync
git pull origin develop

# Save your work
git add -A && git commit -m "feat: describe what you did"

# Push
git push

# See what changed
git status
git log --oneline --graph --all

# Undo last commit (keep changes)
git reset --soft HEAD~1

# Switch branch
git checkout feature/aim2-transcriptomic-drift

# Pull request review (GitHub CLI)
gh pr create --title "Add MODS heatmap visualization" --body "..." --base develop
gh pr list
gh pr merge 3 --squash
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `Permission denied (publickey)` | Re-run `ssh -T git@github.com`; check SSH key is added to GitHub |
| `Large file detected, push rejected` | Add to `.gitignore` or use `git lfs track` |
| `Merge conflict` | Open conflicting file, resolve `<<<<<<` markers, `git add`, `git commit` |
| `conda env create fails` | Try `mamba env create` instead; check channel order |
| `Snakemake rule not found` | Run `snakemake --list` to see all available rules |
| GitHub Actions CI fails | Check the **Actions** tab logs for the specific error |

---

*PDO Drift Characterization Project — mtariqi/pdo-drift*
