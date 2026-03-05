"""
Multi-Omic Drift Score (MODS) Calculator
=========================================
Computes a composite drift score for each PDO line at each passage/condition
by integrating genomic, transcriptomic, epigenomic, and functional distances
from the passage-2 (P2) baseline.

    MODS = 0.30 × GD + 0.35 × TD + 0.15 × ED + 0.20 × FD

Usage
-----
    from pdo_drift.integration import MODSCalculator

    calc = MODSCalculator(config_path="config/config.yaml")
    scores = calc.compute_all(
        vcf_dir="results/variants/",
        rnaseq_dir="results/deseq2/",
        rrbs_dir="results/methylation/",
        drug_dir="results/functional/"
    )
    calc.plot_heatmap(scores, outfile="results/figures/mods_heatmap.pdf")
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import yaml
from scipy.spatial.distance import cosine, euclidean
from scipy.stats import pearsonr
from sklearn.preprocessing import MinMaxScaler

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=FutureWarning)


# ─────────────────────────────────────────────────────────────────────────────
#  Helper: load config
# ─────────────────────────────────────────────────────────────────────────────
def _load_config(path: str | Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


# ─────────────────────────────────────────────────────────────────────────────
#  Component calculators
# ─────────────────────────────────────────────────────────────────────────────
class GenomicDistanceCalculator:
    """
    Compute genomic distance (GD) between a test sample and the P2 baseline.

    Strategy
    --------
    1. Build a binary somatic mutation matrix (genes × samples) from filtered VCFs.
    2. Compute cosine distance between the test and baseline mutation vectors.
    3. Add a CNV burden component: scaled sum of |CNV segment log2 ratio| differences.
    4. GD = 0.6 × cosine_mut_distance + 0.4 × normalized_cnv_distance  (range 0–1)
    """

    def __init__(self, min_vaf: float = 0.05):
        self.min_vaf = min_vaf
        self._scaler = MinMaxScaler()

    def load_mutations(self, vcf_paths: list[Path]) -> pd.DataFrame:
        """
        Parse filtered VCF files and return a gene × sample binary matrix.

        Parameters
        ----------
        vcf_paths : list of Path
            One VCF per sample (filtered MuTect2 output).

        Returns
        -------
        DataFrame of shape (n_genes, n_samples) with 1/0 values.
        """
        import cyvcf2

        records: dict[str, dict[str, int]] = {}  # gene -> {sample_id: 1}
        sample_ids = []

        for vcf_path in vcf_paths:
            sample_id = vcf_path.stem.replace(".filtered", "")
            sample_ids.append(sample_id)
            vcf = cyvcf2.VCF(str(vcf_path))

            for variant in vcf:
                if variant.FILTER and variant.FILTER != "PASS":
                    continue
                af = variant.INFO.get("AF", [0])[0]
                if af < self.min_vaf:
                    continue
                gene = variant.INFO.get("GENE", variant.CHROM)
                records.setdefault(gene, {})[sample_id] = 1

        df = pd.DataFrame(records).T.reindex(columns=sample_ids).fillna(0).astype(int)
        logger.info(f"Loaded {df.shape[0]} mutated genes across {df.shape[1]} samples")
        return df

    def compute_distance(
        self,
        mut_matrix: pd.DataFrame,
        cnv_df: Optional[pd.DataFrame],
        baseline_col: str,
        test_col: str,
    ) -> float:
        """
        Compute GD between baseline and test sample.

        Parameters
        ----------
        mut_matrix  : gene × sample binary mutation matrix
        cnv_df      : DataFrame with columns [sample_id, gene, log2_ratio]
        baseline_col: column name for P2 reference sample
        test_col    : column name for test sample

        Returns
        -------
        GD score in [0, 1]
        """
        # Mutation cosine distance
        b_vec = mut_matrix[baseline_col].values.astype(float)
        t_vec = mut_matrix[test_col].values.astype(float)

        # Guard against zero vectors (no mutations at baseline — early passage)
        if b_vec.sum() == 0 and t_vec.sum() == 0:
            cos_dist = 0.0
        elif b_vec.sum() == 0 or t_vec.sum() == 0:
            cos_dist = 1.0
        else:
            cos_dist = cosine(b_vec, t_vec)

        # CNV distance component
        if cnv_df is not None:
            b_cnv = cnv_df[cnv_df["sample_id"] == baseline_col].set_index("gene")["log2_ratio"]
            t_cnv = cnv_df[cnv_df["sample_id"] == test_col].set_index("gene")["log2_ratio"]
            shared = b_cnv.index.intersection(t_cnv.index)
            if len(shared) > 0:
                cnv_dist_raw = np.mean(np.abs(t_cnv[shared].values - b_cnv[shared].values))
                cnv_dist = np.clip(cnv_dist_raw / 2.0, 0, 1)  # normalise by max log2=2
            else:
                cnv_dist = 0.0
        else:
            cnv_dist = cos_dist  # fallback

        gd = 0.6 * cos_dist + 0.4 * cnv_dist
        return float(np.clip(gd, 0, 1))


class TranscriptomicDistanceCalculator:
    """
    Compute transcriptomic distance (TD) using rlog-normalized counts.

    Strategy
    --------
    1. Load rlog-normalized count matrix (genes × samples) from DESeq2 output.
    2. Select top 5000 most variable genes across all samples.
    3. TD = 1 - Pearson correlation between test and baseline expression vectors.
       (Bounded in [0, 1] since correlation can be negative for drifted samples.)
    """

    def __init__(self, n_top_genes: int = 5000):
        self.n_top_genes = n_top_genes

    def load_expression(self, rlog_csv: Path) -> pd.DataFrame:
        """Load DESeq2 rlog-normalized count matrix."""
        df = pd.read_csv(rlog_csv, index_col=0)
        logger.info(f"Loaded expression matrix: {df.shape[0]} genes × {df.shape[1]} samples")
        return df

    def select_hvg(self, expr_df: pd.DataFrame) -> pd.DataFrame:
        """Retain top highly variable genes by row variance."""
        variances = expr_df.var(axis=1)
        top_genes = variances.nlargest(self.n_top_genes).index
        return expr_df.loc[top_genes]

    def compute_distance(
        self, expr_df: pd.DataFrame, baseline_col: str, test_col: str
    ) -> float:
        """
        Returns TD = 1 - r (Pearson), clipped to [0, 1].
        """
        hvg_df = self.select_hvg(expr_df)
        b_vec = hvg_df[baseline_col].values
        t_vec = hvg_df[test_col].values
        r, _ = pearsonr(b_vec, t_vec)
        td = 1.0 - r
        return float(np.clip(td, 0, 1))


class EpigenomicDistanceCalculator:
    """
    Compute epigenomic distance (ED) from RRBS beta-value matrices.

    Strategy
    --------
    1. Load per-sample CpG beta-value matrices (site × sample).
    2. Retain sites with ≥10× coverage in all samples.
    3. Focus on differentially variable CpGs (variance > 0.05).
    4. ED = mean absolute beta-value difference between test and baseline.
       Normalised to [0, 1] (maximum theoretical delta = 1.0).
    """

    def __init__(self, min_coverage: int = 10, min_variance: float = 0.05):
        self.min_coverage = min_coverage
        self.min_variance = min_variance

    def load_methylation(self, beta_csv: Path) -> pd.DataFrame:
        """
        Load a CpG × sample beta-value matrix.
        Missing values (< min_coverage) should be NaN.
        """
        df = pd.read_csv(beta_csv, index_col=0)
        logger.info(f"Loaded methylation matrix: {df.shape[0]} CpGs × {df.shape[1]} samples")
        return df

    def filter_sites(self, beta_df: pd.DataFrame) -> pd.DataFrame:
        """Drop sites with any missing coverage; retain variable sites."""
        complete = beta_df.dropna()
        variable = complete[complete.var(axis=1) >= self.min_variance]
        logger.info(
            f"After filtering: {variable.shape[0]} variable CpG sites retained "
            f"(from {beta_df.shape[0]} total)"
        )
        return variable

    def compute_distance(
        self, beta_df: pd.DataFrame, baseline_col: str, test_col: str
    ) -> float:
        filtered = self.filter_sites(beta_df)
        if filtered.empty:
            logger.warning("No variable CpG sites found — returning ED=0")
            return 0.0
        delta = np.abs(filtered[test_col].values - filtered[baseline_col].values)
        ed = float(np.mean(delta))
        return np.clip(ed, 0, 1)


class FunctionalDistanceCalculator:
    """
    Compute functional distance (FD) from drug sensitivity AUC matrices.

    Strategy
    --------
    1. Load AUC matrix (drugs × samples) from drug screen results.
    2. FD = mean normalized |ΔAUC| across the 48-compound panel.
       Normalised by compound-level range across all samples.
    3. Compound IC50 shift flags: report compounds where IC50 shifts >3-fold
       as an auxiliary output.
    """

    def __init__(self, ic50_shift_threshold: float = 3.0):
        self.ic50_shift_threshold = ic50_shift_threshold

    def load_auc(self, auc_csv: Path) -> pd.DataFrame:
        """Load drug × sample AUC matrix (range 0–1 per drug)."""
        df = pd.read_csv(auc_csv, index_col=0)
        logger.info(f"Loaded AUC matrix: {df.shape[0]} compounds × {df.shape[1]} samples")
        return df

    def compute_distance(
        self, auc_df: pd.DataFrame, baseline_col: str, test_col: str
    ) -> float:
        delta = np.abs(auc_df[test_col].values - auc_df[baseline_col].values)
        fd = float(np.mean(delta))
        return np.clip(fd, 0, 1)

    def flag_shifted_drugs(
        self,
        ic50_df: pd.DataFrame,
        baseline_col: str,
        test_col: str,
    ) -> pd.DataFrame:
        """
        Return a DataFrame of drugs where IC50 has shifted beyond threshold.
        """
        fold_change = ic50_df[test_col] / (ic50_df[baseline_col] + 1e-9)
        flagged = ic50_df[fold_change.abs() >= self.ic50_shift_threshold].copy()
        flagged["fold_change"] = fold_change[flagged.index]
        return flagged


# ─────────────────────────────────────────────────────────────────────────────
#  Master MODS Calculator
# ─────────────────────────────────────────────────────────────────────────────
class MODSCalculator:
    """
    Orchestrator that computes the Multi-Omic Drift Score for all samples.

    Parameters
    ----------
    config_path : str or Path
        Path to config/config.yaml

    Examples
    --------
    >>> calc = MODSCalculator("config/config.yaml")
    >>> scores = calc.compute_all(...)
    >>> calc.plot_heatmap(scores)
    >>> calc.plot_passage_trajectory(scores, line="PDO001")
    """

    def __init__(self, config_path: str | Path = "config/config.yaml"):
        self.cfg = _load_config(config_path)
        w = self.cfg["mods"]["weights"]
        self.weights = {
            "genomic": w["genomic"],
            "transcriptomic": w["transcriptomic"],
            "epigenomic": w["epigenomic"],
            "functional": w["functional"],
        }
        self.threshold = self.cfg["mods"]["drift_threshold"]
        self._gd_calc = GenomicDistanceCalculator(
            min_vaf=self.cfg["wes"]["min_allele_freq"]
        )
        self._td_calc = TranscriptomicDistanceCalculator()
        self._ed_calc = EpigenomicDistanceCalculator(
            min_coverage=self.cfg["rrbs"]["min_coverage"]
        )
        self._fd_calc = FunctionalDistanceCalculator(
            ic50_shift_threshold=self.cfg["drug_panel"]["ic50_shift_threshold"]
        )

    def compute_all(
        self,
        vcf_dir: str | Path,
        rnaseq_dir: str | Path,
        rrbs_dir: str | Path,
        drug_dir: str | Path,
    ) -> pd.DataFrame:
        """
        Compute MODS for all samples relative to their P2 baseline.

        Returns
        -------
        DataFrame with columns:
            sample_id, pdo_line, passage, condition,
            GD, TD, ED, FD, MODS, drift_flag
        """
        vcf_dir = Path(vcf_dir)
        rnaseq_dir = Path(rnaseq_dir)
        rrbs_dir = Path(rrbs_dir)
        drug_dir = Path(drug_dir)

        samples = pd.read_csv(self.cfg["samples"], sep="\t")
        results = []

        # ── Load shared matrices ──────────────────────────────
        vcf_paths = sorted(vcf_dir.glob("*.filtered.vcf.gz"))
        if vcf_paths:
            mut_matrix = self._gd_calc.load_mutations(vcf_paths)
            cnv_df = self._load_cnv(vcf_dir.parent / "cnv")
        else:
            mut_matrix, cnv_df = None, None

        rlog_csv = rnaseq_dir / "rlog_normalized_counts.csv"
        expr_df = self._td_calc.load_expression(rlog_csv) if rlog_csv.exists() else None

        beta_csv = rrbs_dir / "cpg_beta_matrix.csv"
        beta_df = self._ed_calc.load_methylation(beta_csv) if beta_csv.exists() else None

        auc_csv = drug_dir / "auc_matrix.csv"
        auc_df = self._fd_calc.load_auc(auc_csv) if auc_csv.exists() else None

        # ── Iterate samples ───────────────────────────────────
        for _, row in samples.iterrows():
            pdo_line = row["pdo_line"]
            passage = row["passage"]
            sample_id = row["sample_id"]
            condition = row.get("condition", "fresh")

            # Find the P2 baseline for this line
            baseline_mask = (
                (samples["pdo_line"] == pdo_line)
                & (samples["passage"] == 2)
                & (samples["condition"] == "fresh")
            )
            baseline_ids = samples.loc[baseline_mask, "sample_id"].tolist()
            if not baseline_ids or passage == 2 and condition == "fresh":
                continue
            baseline_id = baseline_ids[0]

            # ── Compute component distances ───────────────────
            gd = (
                self._gd_calc.compute_distance(
                    mut_matrix, cnv_df, baseline_id, sample_id
                )
                if mut_matrix is not None and sample_id in mut_matrix.columns
                else np.nan
            )
            td = (
                self._td_calc.compute_distance(expr_df, baseline_id, sample_id)
                if expr_df is not None and sample_id in expr_df.columns
                else np.nan
            )
            ed = (
                self._ed_calc.compute_distance(beta_df, baseline_id, sample_id)
                if beta_df is not None and sample_id in beta_df.columns
                else np.nan
            )
            fd = (
                self._fd_calc.compute_distance(auc_df, baseline_id, sample_id)
                if auc_df is not None and sample_id in auc_df.columns
                else np.nan
            )

            # ── Weighted MODS ─────────────────────────────────
            components = {"GD": gd, "TD": td, "ED": ed, "FD": fd}
            weight_keys = {"GD": "genomic", "TD": "transcriptomic", "ED": "epigenomic", "FD": "functional"}
            valid = {k: v for k, v in components.items() if not np.isnan(v)}
            if valid:
                total_w = sum(self.weights[weight_keys[k]] for k in valid)
                mods = sum(v * self.weights[weight_keys[k]] for k, v in valid.items()) / total_w
            else:
                mods = np.nan

            results.append({
                "sample_id": sample_id,
                "pdo_line": pdo_line,
                "passage": passage,
                "condition": condition,
                **components,
                "MODS": round(float(mods), 4) if not np.isnan(mods) else np.nan,
                "drift_flag": mods > self.threshold if not np.isnan(mods) else False,
            })

        df = pd.DataFrame(results).sort_values(["pdo_line", "passage"])
        logger.info(
            f"MODS computed for {len(df)} samples. "
            f"Flagged: {df['drift_flag'].sum()} ({df['drift_flag'].mean():.1%})"
        )
        return df

    # ── Plotting ──────────────────────────────────────────────
    def plot_heatmap(
        self,
        scores: pd.DataFrame,
        outfile: str | Path = "results/figures/mods_heatmap.pdf",
        figsize: tuple = (14, 8),
    ) -> None:
        """
        Plot a heatmap of MODS scores — PDO lines × passages.
        Annotates cells with the score value; flags cells exceeding threshold in red.
        """
        pivot = scores.pivot_table(
            index="pdo_line", columns="passage", values="MODS", aggfunc="mean"
        )
        fig, ax = plt.subplots(figsize=figsize)
        cmap = sns.diverging_palette(220, 20, as_cmap=True)

        sns.heatmap(
            pivot,
            ax=ax,
            cmap="YlOrRd",
            vmin=0,
            vmax=1,
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            linecolor="white",
            cbar_kws={"label": "MODS (0–1)", "shrink": 0.8},
        )
        ax.set_title(
            "Multi-Omic Drift Score (MODS) — PDO Lines × Passage",
            fontsize=14,
            fontweight="bold",
            pad=12,
        )
        ax.set_xlabel("Passage Number", fontsize=11)
        ax.set_ylabel("PDO Line", fontsize=11)

        # Highlight cells above threshold
        for i, line in enumerate(pivot.index):
            for j, passage in enumerate(pivot.columns):
                val = pivot.at[line, passage]
                if not np.isnan(val) and val > self.threshold:
                    ax.add_patch(
                        plt.Rectangle(
                            (j, i), 1, 1,
                            fill=False,
                            edgecolor="red",
                            lw=2.5,
                        )
                    )

        plt.tight_layout()
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"MODS heatmap saved to {outfile}")

    def plot_passage_trajectory(
        self,
        scores: pd.DataFrame,
        line: Optional[str] = None,
        outfile: str | Path = "results/figures/mods_trajectory.pdf",
        figsize: tuple = (12, 6),
    ) -> None:
        """
        Line plot of MODS and component scores across passages for one or all lines.
        """
        df = scores if line is None else scores[scores["pdo_line"] == line]
        if df.empty:
            logger.warning(f"No data for line: {line}")
            return

        fig, axes = plt.subplots(1, 2, figsize=figsize)
        palette = sns.color_palette("tab10")

        # Left: MODS trajectory per line
        ax = axes[0]
        for i, (lname, grp) in enumerate(df.groupby("pdo_line")):
            fresh = grp[grp["condition"] == "fresh"].sort_values("passage")
            ax.plot(
                fresh["passage"], fresh["MODS"],
                marker="o", linewidth=2,
                color=palette[i % len(palette)],
                label=lname,
            )
        ax.axhline(self.threshold, color="red", linestyle="--", alpha=0.7,
                   label=f"Drift threshold ({self.threshold})")
        ax.set_xlabel("Passage Number", fontsize=11)
        ax.set_ylabel("MODS", fontsize=11)
        ax.set_title("MODS Trajectory Across Passages", fontsize=12, fontweight="bold")
        ax.legend(fontsize=8, ncol=2)
        ax.set_ylim(0, 1)

        # Right: Component breakdown (stacked area for single line)
        ax2 = axes[1]
        if line and not df.empty:
            fresh = df[df["condition"] == "fresh"].sort_values("passage")
            for comp, color in zip(
                ["GD", "TD", "ED", "FD"],
                [palette[0], palette[1], palette[2], palette[3]],
            ):
                ax2.plot(
                    fresh["passage"], fresh[comp],
                    marker="s", linewidth=1.8,
                    color=color, label=comp,
                )
            ax2.set_title(f"Component Scores — {line}", fontsize=12, fontweight="bold")
            ax2.set_xlabel("Passage Number", fontsize=11)
            ax2.set_ylabel("Distance Score (0–1)", fontsize=11)
            ax2.legend(fontsize=10)
            ax2.set_ylim(0, 1)
        else:
            ax2.set_visible(False)

        plt.tight_layout()
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Trajectory plot saved to {outfile}")

    # ── Private helpers ───────────────────────────────────────
    @staticmethod
    def _load_cnv(cnv_dir: Path) -> Optional[pd.DataFrame]:
        """Load and concatenate CNVKit .cnr files into a tidy DataFrame."""
        frames = []
        for f in cnv_dir.glob("*.cnr"):
            df = pd.read_csv(f, sep="\t")
            df["sample_id"] = f.stem
            frames.append(df[["sample_id", "gene", "log2"]])
        if not frames:
            return None
        combined = pd.concat(frames, ignore_index=True)
        combined.rename(columns={"log2": "log2_ratio"}, inplace=True)
        return combined
