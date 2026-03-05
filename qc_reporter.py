"""
PDO QC Module
==============
Sample-level and run-level quality control checks for all assays.
Generates QC reports and flags samples that fail minimum thresholds.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)

QC_THRESHOLDS = {
    "viability_pct": 85.0,          # Minimum cell viability at harvest
    "rin": 7.0,                      # Minimum RNA integrity number
    "wes_mean_coverage": 100.0,      # Minimum WES coverage
    "wes_target_pct_100x": 90.0,     # % target bases at ≥100x
    "rnaseq_mapped_pct": 70.0,       # % reads mapped (STAR)
    "rnaseq_total_reads_M": 30.0,    # Minimum total reads (millions)
    "rrbs_cpg_covered_pct": 80.0,    # % CpG sites with ≥10x coverage
    "mycoplasma": "negative",        # Required mycoplasma status
}


class QCReporter:
    """
    Aggregate QC metrics from all assays and generate a unified QC report.

    Parameters
    ----------
    samples_tsv : path to config/samples.tsv

    Usage
    -----
    >>> qc = QCReporter("config/samples.tsv")
    >>> qc.load_fastp_reports("results/qc/fastp/")
    >>> qc.load_star_logs("data/star/")
    >>> qc.load_picard_metrics("results/qc/picard/")
    >>> report = qc.generate_report()
    >>> qc.plot_qc_summary()
    """

    def __init__(self, samples_tsv: str | Path, thresholds: dict = None):
        self.samples = pd.read_csv(samples_tsv, sep="\t", index_col="sample_id")
        self.thresholds = thresholds or QC_THRESHOLDS
        self.metrics: dict[str, dict] = {sid: {} for sid in self.samples.index}

    def load_fastp_reports(self, fastp_dir: str | Path) -> None:
        """Parse fastp JSON reports for read QC metrics."""
        for json_file in Path(fastp_dir).glob("*.json"):
            sample_id = json_file.stem
            if sample_id not in self.metrics:
                continue
            with open(json_file) as fh:
                data = json.load(fh)
            summary = data.get("summary", {})
            after = summary.get("after_filtering", {})
            self.metrics[sample_id].update({
                "total_reads_M": after.get("total_reads", 0) / 1e6,
                "q30_rate": after.get("q30_rate", 0),
                "duplication_rate": data.get("duplication", {}).get("rate", 0),
            })

    def load_star_logs(self, star_dir: str | Path) -> None:
        """Parse STAR Log.final.out for alignment statistics."""
        for log_file in Path(star_dir).rglob("Log.final.out"):
            sample_id = log_file.parent.name
            if sample_id not in self.metrics:
                continue
            with open(log_file) as fh:
                lines = fh.readlines()
            for line in lines:
                if "Uniquely mapped reads %" in line:
                    pct = float(line.split("|")[-1].strip().rstrip("%"))
                    self.metrics[sample_id]["rnaseq_mapped_pct"] = pct
                if "Number of input reads" in line:
                    n = int(line.split("|")[-1].strip())
                    self.metrics[sample_id].setdefault("total_reads_M", n / 1e6)

    def load_picard_metrics(self, picard_dir: str | Path) -> None:
        """Parse Picard duplication metrics files."""
        for met_file in Path(picard_dir).glob("*.dup_metrics.txt"):
            sample_id = met_file.stem.replace(".dup_metrics", "")
            if sample_id not in self.metrics:
                continue
            # Skip comment lines, parse METRICS CLASS header + data
            try:
                df = pd.read_csv(met_file, sep="\t", comment="#", nrows=1,
                                 skiprows=lambda x: x > 7 and x % 2 == 0)
                if "PERCENT_DUPLICATION" in df.columns:
                    self.metrics[sample_id]["duplication_rate"] = float(
                        df["PERCENT_DUPLICATION"].iloc[0]
                    )
            except Exception:
                pass

    def generate_report(
        self,
        outfile: str | Path = "results/qc/qc_report.csv",
    ) -> pd.DataFrame:
        """
        Compile all QC metrics and apply pass/fail flags.

        Returns
        -------
        DataFrame with QC metrics and per-sample pass/fail status.
        """
        rows = []
        for sample_id, m in self.metrics.items():
            row = {"sample_id": sample_id, **m}
            flags = []

            # Check numeric thresholds
            for metric, threshold in self.thresholds.items():
                if isinstance(threshold, (int, float)) and metric in m:
                    if m[metric] < threshold:
                        flags.append(f"{metric}<{threshold}")
                elif isinstance(threshold, str) and metric in m:
                    if str(m[metric]).lower() != threshold.lower():
                        flags.append(f"{metric}!={threshold}")

            row["qc_flags"] = "; ".join(flags) if flags else "PASS"
            row["qc_pass"] = len(flags) == 0
            rows.append(row)

        report = pd.DataFrame(rows).set_index("sample_id")
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        report.to_csv(outfile)

        n_pass = report["qc_pass"].sum()
        n_total = len(report)
        logger.info(f"QC report: {n_pass}/{n_total} samples passed all thresholds")
        return report

    def plot_qc_summary(
        self,
        report: Optional[pd.DataFrame] = None,
        outfile: str | Path = "results/figures/qc_summary.pdf",
    ) -> None:
        """
        Multi-panel QC overview plot.
        """
        if report is None:
            report = self.generate_report()

        numeric_cols = [c for c in report.columns
                        if report[c].dtype in [np.float64, np.int64]
                        and c not in ["qc_pass"]]
        if not numeric_cols:
            logger.warning("No numeric QC metrics to plot.")
            return

        n_plots = len(numeric_cols)
        n_cols = 3
        n_rows = int(np.ceil(n_plots / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, n_rows * 4))
        axes = axes.flatten()

        for i, col in enumerate(numeric_cols):
            ax = axes[i]
            vals = report[col].dropna()
            colors = ["#2E6DB4" if v else "#C95C00"
                      for v in report.loc[vals.index, "qc_pass"]]
            ax.bar(range(len(vals)), vals.values, color=colors, alpha=0.8)

            # Draw threshold line if available
            threshold = self.thresholds.get(col)
            if isinstance(threshold, (int, float)):
                ax.axhline(threshold, color="red", linestyle="--",
                           linewidth=1.2, label=f"Threshold: {threshold}")
                ax.legend(fontsize=8)

            ax.set_title(col.replace("_", " ").title(), fontsize=10, fontweight="bold")
            ax.set_ylabel(col, fontsize=8)
            ax.set_xticks([])
            ax.tick_params(labelsize=8)

        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        fig.suptitle("PDO Sample QC Summary", fontsize=14, fontweight="bold", y=1.02)
        plt.tight_layout()
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outfile, dpi=200, bbox_inches="tight")
        plt.close()
        logger.info(f"QC summary plot saved to {outfile}")
