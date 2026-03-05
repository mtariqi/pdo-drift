"""
Transcriptomic Drift Tracker
==============================
Computes passage-resolved transcriptomic drift metrics from bulk RNA-seq data.

Features
--------
- PyDESeq2 differential expression (Python-native DESeq2)
- GSEA with gseapy (Hallmark, KEGG, WikiPathways)
- Drift trajectory visualization (PCA, UMAP)
- Growth kinetics gene signature identification (LASSO)
- EMT and cell-cycle pathway scoring
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    PYDESEQ2_AVAILABLE = True
except ImportError:
    PYDESEQ2_AVAILABLE = False
    logging.warning("pydeseq2 not installed. DESeq2 functions disabled.")

try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False
    logging.warning("gseapy not installed. GSEA functions disabled.")

logger = logging.getLogger(__name__)


class DriftTracker:
    """
    Tracks transcriptomic drift across serial PDO passages.

    Parameters
    ----------
    counts_csv : path to raw count matrix (genes × samples)
    metadata_tsv : path to sample metadata TSV

    Examples
    --------
    >>> tracker = DriftTracker("results/counts/merged_counts.csv", "config/samples.tsv")
    >>> tracker.run_deseq2(contrast=["passage", "P20", "P2"])
    >>> tracker.run_gsea()
    >>> tracker.plot_pca_drift()
    >>> sig = tracker.identify_growth_signature()
    """

    def __init__(self, counts_csv: str | Path, metadata_tsv: str | Path):
        self.counts = self._load_counts(counts_csv)
        self.metadata = pd.read_csv(metadata_tsv, sep="\t", index_col="sample_id")
        self._align_samples()
        self.dge_results: Optional[pd.DataFrame] = None
        self.gsea_results: dict[str, pd.DataFrame] = {}
        logger.info(
            f"DriftTracker initialized: {self.counts.shape[0]} genes, "
            f"{self.counts.shape[1]} samples"
        )

    # ── Data loading ──────────────────────────────────────────
    @staticmethod
    def _load_counts(path: str | Path) -> pd.DataFrame:
        """Load merged count matrix and strip ENSEMBL version suffixes."""
        df = pd.read_csv(path, index_col=0)
        df.index = df.index.str.replace(r"\.\d+$", "", regex=True)
        # Remove low-count genes (at least 10 reads in ≥2 samples)
        keep = (df >= 10).sum(axis=1) >= 2
        df = df.loc[keep]
        return df

    def _align_samples(self):
        """Ensure count matrix and metadata share the same samples."""
        shared = self.counts.columns.intersection(self.metadata.index)
        self.counts = self.counts[shared]
        self.metadata = self.metadata.loc[shared]
        if len(shared) < len(self.counts.columns):
            logger.warning(
                f"Dropped {len(self.counts.columns) - len(shared)} samples "
                f"not present in metadata."
            )

    # ── DESeq2 differential expression ───────────────────────
    def run_deseq2(
        self,
        design_factor: str = "passage",
        baseline_level: str = "P2",
        outfile: str | Path = "results/deseq2/drift_dge_results.csv",
    ) -> pd.DataFrame:
        """
        Run PyDESeq2 to identify genes differentially expressed vs. P2 baseline.

        Parameters
        ----------
        design_factor   : metadata column to use as the design variable
        baseline_level  : reference level (e.g., "P2")
        outfile         : CSV output path

        Returns
        -------
        DataFrame with columns: gene, baseMean, log2FC, lfcSE, stat, pvalue, padj
        """
        if not PYDESEQ2_AVAILABLE:
            raise ImportError("Install pydeseq2: pip install pydeseq2")

        meta = self.metadata[[design_factor]].copy()
        meta[design_factor] = meta[design_factor].astype("category")
        meta[design_factor] = meta[design_factor].cat.reorder_categories(
            [baseline_level]
            + [c for c in meta[design_factor].cat.categories if c != baseline_level]
        )

        dds = DeseqDataSet(
            counts=self.counts.T,
            clinical=meta,
            design_factors=design_factor,
            refit_cooks=True,
            n_cpus=8,
        )
        dds.deseq2()

        all_results = []
        for level in meta[design_factor].cat.categories:
            if level == baseline_level:
                continue
            stat_res = DeseqStats(
                dds,
                contrast=[design_factor, level, baseline_level],
                alpha=0.05,
                cooks_filter=True,
                independent_filter=True,
            )
            stat_res.summary()
            df = stat_res.results_df.copy()
            df["contrast"] = f"{level}_vs_{baseline_level}"
            all_results.append(df)

        combined = pd.concat(all_results)
        combined.index.name = "gene"
        combined.reset_index(inplace=True)

        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        combined.to_csv(outfile, index=False)
        self.dge_results = combined
        logger.info(f"DESeq2 complete. Saved to {outfile}")
        return combined

    # ── GSEA ──────────────────────────────────────────────────
    def run_gsea(
        self,
        gene_sets: list[str] = None,
        outdir: str | Path = "results/gsea/",
    ) -> dict[str, pd.DataFrame]:
        """
        Run pre-ranked GSEA for each contrast using gseapy.

        Parameters
        ----------
        gene_sets : list of MSigDB gene set names
        outdir    : directory to save GSEA output

        Returns
        -------
        Dictionary of {contrast: results_DataFrame}
        """
        if not GSEAPY_AVAILABLE:
            raise ImportError("Install gseapy: pip install gseapy")
        if self.dge_results is None:
            raise RuntimeError("Run run_deseq2() first.")

        if gene_sets is None:
            gene_sets = ["MSigDB_Hallmark_2023", "KEGG_2021_Human", "WikiPathway_2023_Human"]

        Path(outdir).mkdir(parents=True, exist_ok=True)
        results = {}

        for contrast in self.dge_results["contrast"].unique():
            subset = self.dge_results[self.dge_results["contrast"] == contrast].copy()
            subset = subset.dropna(subset=["stat"])

            # Rank metric: signed -log10(pvalue) × sign(LFC)
            subset["rank_metric"] = (
                np.sign(subset["log2FoldChange"])
                * -np.log10(subset["pvalue"].clip(lower=1e-300))
            )
            ranking = subset.set_index("gene")["rank_metric"].sort_values(ascending=False)

            for gs in gene_sets:
                try:
                    res = gp.prerank(
                        rnk=ranking,
                        gene_sets=gs,
                        threads=4,
                        permutation_num=1000,
                        outdir=str(Path(outdir) / contrast / gs.replace(" ", "_")),
                        seed=42,
                        verbose=False,
                    )
                    key = f"{contrast}_{gs}"
                    results[key] = res.res2d
                    logger.info(f"GSEA complete: {key} — {len(res.res2d)} gene sets tested")
                except Exception as exc:
                    logger.warning(f"GSEA failed for {contrast}/{gs}: {exc}")

        self.gsea_results = results
        return results

    # ── PCA drift visualization ───────────────────────────────
    def plot_pca_drift(
        self,
        outfile: str | Path = "results/figures/pca_drift.pdf",
        color_by: str = "passage",
        shape_by: str = "pdo_line",
        figsize: tuple = (10, 8),
    ) -> None:
        """
        PCA of rlog-normalized counts, colored by passage.
        Draws arrows showing passage trajectory per line.
        """
        # Log-normalize for PCA
        log_counts = np.log1p(self.counts)
        scaler = StandardScaler()
        scaled = scaler.fit_transform(log_counts.T)

        pca = PCA(n_components=10, random_state=42)
        coords = pca.fit_transform(scaled)
        var_exp = pca.explained_variance_ratio_

        pca_df = pd.DataFrame(
            coords[:, :2], columns=["PC1", "PC2"], index=self.counts.columns
        )
        pca_df = pca_df.join(self.metadata[[color_by, shape_by]])

        fig, ax = plt.subplots(figsize=figsize)
        palette = dict(zip(
            sorted(pca_df[color_by].unique()),
            sns.color_palette("viridis", n_colors=pca_df[color_by].nunique())
        ))

        for line, grp in pca_df.groupby(shape_by):
            sorted_grp = grp.sort_values(color_by)
            ax.plot(
                sorted_grp["PC1"], sorted_grp["PC2"],
                "k-", alpha=0.2, linewidth=0.8, zorder=1
            )

        scatter = ax.scatter(
            pca_df["PC1"], pca_df["PC2"],
            c=pca_df[color_by].map(palette),
            s=80, zorder=2, edgecolors="white", linewidths=0.5
        )

        ax.set_xlabel(f"PC1 ({var_exp[0]:.1%} variance)", fontsize=11)
        ax.set_ylabel(f"PC2 ({var_exp[1]:.1%} variance)", fontsize=11)
        ax.set_title("PCA of PDO Transcriptomes Across Passages", fontsize=13, fontweight="bold")

        from matplotlib.patches import Patch
        legend_handles = [Patch(color=c, label=lv) for lv, c in palette.items()]
        ax.legend(handles=legend_handles, title=color_by.capitalize(), fontsize=9, ncol=2)

        plt.tight_layout()
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"PCA drift plot saved to {outfile}")

    # ── Growth kinetics gene signature ────────────────────────
    def identify_growth_signature(
        self,
        growth_rates: pd.Series,
        outfile: str | Path = "results/transcriptomics/growth_signature.csv",
        max_iter: int = 10000,
    ) -> pd.DataFrame:
        """
        Use LASSO regression to identify genes predictive of PDO growth rate.

        Parameters
        ----------
        growth_rates : Series indexed by sample_id with doubling time (hours) values
        outfile      : output CSV for the growth rate gene signature

        Returns
        -------
        DataFrame of signature genes with LASSO coefficients
        """
        # Align expression with growth rate samples
        shared = self.counts.columns.intersection(growth_rates.index)
        X = np.log1p(self.counts[shared].T).values
        y = growth_rates[shared].values

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        lasso = LassoCV(cv=5, max_iter=max_iter, random_state=42, n_jobs=4)
        lasso.fit(X_scaled, y)

        coef_df = pd.DataFrame({
            "gene": self.counts.index,
            "lasso_coefficient": lasso.coef_,
        }).query("lasso_coefficient != 0").sort_values(
            "lasso_coefficient", key=abs, ascending=False
        )

        coef_df["direction"] = np.where(
            coef_df["lasso_coefficient"] > 0, "fast_growth", "slow_growth"
        )
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        coef_df.to_csv(outfile, index=False)
        logger.info(
            f"Growth signature: {len(coef_df)} genes selected "
            f"(α={lasso.alpha_:.4f}, CV R²={lasso.score(X_scaled, y):.3f})"
        )
        return coef_df

    # ── Pathway activity scoring ──────────────────────────────
    def score_pathway_activity(
        self,
        pathways: dict[str, list[str]] = None,
        outfile: str | Path = "results/transcriptomics/pathway_scores.csv",
    ) -> pd.DataFrame:
        """
        Score pathway activity per sample using mean z-scored gene expression.

        Parameters
        ----------
        pathways : dict mapping pathway_name -> list of gene symbols

        Returns
        -------
        DataFrame of shape (samples × pathways)
        """
        if pathways is None:
            pathways = {
                "WNT_signaling": ["CTNNB1", "APC", "AXIN1", "LRP5", "LRP6", "TCF7L2", "MYC"],
                "MAPK_ERK": ["KRAS", "BRAF", "MAP2K1", "MAPK1", "MAPK3", "ELK1"],
                "PI3K_AKT_mTOR": ["PIK3CA", "PTEN", "AKT1", "MTOR", "RPS6KB1", "EIF4EBP1"],
                "EMT": ["CDH1", "CDH2", "VIM", "FN1", "TWIST1", "ZEB1", "SNAI1"],
                "Cell_Cycle": ["CCND1", "CDK4", "CDK6", "RB1", "E2F1", "CDKN2A"],
                "DNA_Damage_Response": ["TP53", "BRCA1", "BRCA2", "ATM", "CHEK1", "CHEK2"],
                "Apoptosis": ["BCL2", "BAX", "CASP3", "CASP9", "APAF1", "CYCS"],
                "ER_Stress_UPR": ["ATF4", "DDIT3", "XBP1", "ERN1", "EIF2AK3", "ATF6"],
            }

        log_counts = np.log1p(self.counts)
        scores = {}

        for pathway, genes in pathways.items():
            present_genes = [g for g in genes if g in log_counts.index]
            if not present_genes:
                logger.warning(f"No genes found for pathway: {pathway}")
                continue
            subset = log_counts.loc[present_genes]
            # Z-score each gene across samples, then take mean per sample
            z_scored = ((subset.T - subset.mean(axis=1)) / (subset.std(axis=1) + 1e-9)).T
            scores[pathway] = z_scored.mean(axis=0)

        score_df = pd.DataFrame(scores)
        score_df.index.name = "sample_id"
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        score_df.to_csv(outfile)
        return score_df
