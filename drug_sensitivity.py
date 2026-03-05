"""
Drug Sensitivity Profiler
==========================
Fits dose-response curves to organoid viability data, computes IC50 / AUC,
and detects functional drift across passages and freeze-thaw conditions.

Workflow
--------
1. Load raw viability plate data (CSV or Excel)
2. Normalize to DMSO and positive controls
3. Fit 4-parameter logistic (4PL) model per drug × sample
4. Compute IC50, AUC, Emax
5. Compare across passages → functional drift index
6. Plot dose-response grids and heatmaps

Dependencies: numpy, pandas, scipy, matplotlib, seaborn
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore", category=OptimizeWarning)
logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
#  4-Parameter Logistic Model
# ─────────────────────────────────────────────────────────────────────────────
def four_param_logistic(x: np.ndarray, top: float, bottom: float,
                        ic50: float, slope: float) -> np.ndarray:
    """
    4PL dose-response model.

    Parameters
    ----------
    x      : concentration (nM)
    top    : upper asymptote (max viability, ~100)
    bottom : lower asymptote (Emax)
    ic50   : inflection point (nM)
    slope  : Hill coefficient

    Returns
    -------
    Predicted viability (%)
    """
    return bottom + (top - bottom) / (1.0 + (ic50 / (x + 1e-9)) ** slope)


def fit_4pl(
    concentrations: np.ndarray,
    viabilities: np.ndarray,
    bounds: tuple = ((0, -10, 0.001, 0.1), (200, 110, 1e6, 10)),
) -> dict:
    """
    Fit 4PL model and return curve parameters + goodness-of-fit.

    Returns
    -------
    dict with keys: top, bottom, ic50_nM, slope, r_squared, auc, emax
    """
    try:
        p0 = [100, 0, concentrations.mean(), 1.0]
        popt, pcov = curve_fit(
            four_param_logistic, concentrations, viabilities,
            p0=p0, bounds=bounds, maxfev=10000
        )
        top, bottom, ic50, slope = popt

        # R²
        y_pred = four_param_logistic(concentrations, *popt)
        ss_res = np.sum((viabilities - y_pred) ** 2)
        ss_tot = np.sum((viabilities - viabilities.mean()) ** 2)
        r2 = 1 - ss_res / (ss_tot + 1e-9)

        # AUC (trapezoidal on log10 concentration axis, normalized 0–1)
        log_conc = np.log10(concentrations + 1e-9)
        norm_viab = np.clip(viabilities / 100.0, 0, 1)
        auc = float(np.trapz(norm_viab, log_conc) / (log_conc[-1] - log_conc[0]))

        return {
            "top": round(float(top), 3),
            "bottom": round(float(bottom), 3),
            "ic50_nM": round(float(ic50), 4),
            "slope": round(float(slope), 3),
            "r_squared": round(float(r2), 4),
            "auc": round(float(auc), 4),
            "emax": round(float(bottom), 3),
            "fit_success": True,
        }
    except Exception as exc:
        logger.debug(f"4PL fit failed: {exc}")
        return {
            "top": np.nan, "bottom": np.nan, "ic50_nM": np.nan,
            "slope": np.nan, "r_squared": np.nan, "auc": np.nan,
            "emax": np.nan, "fit_success": False,
        }


# ─────────────────────────────────────────────────────────────────────────────
#  DrugSensitivityProfiler
# ─────────────────────────────────────────────────────────────────────────────
class DrugSensitivityProfiler:
    """
    End-to-end drug sensitivity analysis for PDO passages.

    Parameters
    ----------
    config : dict from config/config.yaml (drug_panel section)

    Examples
    --------
    >>> profiler = DrugSensitivityProfiler(config["drug_panel"])
    >>> profiler.load_plate_data("results/functional/raw_viability.csv")
    >>> profiler.normalize()
    >>> profiler.fit_all_curves()
    >>> profiler.compute_drift(baseline_passage=2)
    >>> profiler.plot_heatmap()
    >>> profiler.plot_dose_response_grid("PDO001", passages=[2, 5, 10, 15, 20])
    """

    def __init__(self, config: dict):
        self.config = config
        self.concentrations = np.array(config["concentrations"], dtype=float)  # nM
        self.ic50_threshold = config["ic50_shift_threshold"]
        self.raw_data: Optional[pd.DataFrame] = None
        self.normalized: Optional[pd.DataFrame] = None
        self.curve_params: Optional[pd.DataFrame] = None
        self.drift_summary: Optional[pd.DataFrame] = None

    # ── Data loading ──────────────────────────────────────────
    def load_plate_data(self, path: str | Path) -> pd.DataFrame:
        """
        Load raw viability data.

        Expected columns: sample_id, pdo_line, passage, condition,
                          drug, replicate, conc_1, conc_2, ... conc_N
        (concentrations as separate columns, values = % viability raw)
        """
        df = pd.read_csv(path)
        required = {"sample_id", "pdo_line", "passage", "drug"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        self.raw_data = df
        logger.info(f"Loaded viability data: {len(df)} rows, {df['drug'].nunique()} drugs")
        return df

    def normalize(
        self,
        dmso_label: str = "DMSO",
        pos_ctrl_label: str = "Staurosporine_10uM",
    ) -> pd.DataFrame:
        """
        Normalize raw viability to percentage relative to DMSO (100%) and
        positive control (0%).

        Returns
        -------
        DataFrame with normalized viability columns
        """
        if self.raw_data is None:
            raise RuntimeError("Load data first with load_plate_data()")

        df = self.raw_data.copy()
        conc_cols = [c for c in df.columns if c.startswith("conc_")]

        # Plate-level normalization per sample_id × plate_id
        normalized_frames = []
        for (sample_id, plate_id), group in df.groupby(["sample_id", "passage"]):
            dmso_rows = group[group["drug"] == dmso_label][conc_cols].values
            pos_rows = group[group["drug"] == pos_ctrl_label][conc_cols].values

            dmso_mean = dmso_rows.mean() if len(dmso_rows) > 0 else 100.0
            pos_mean = pos_rows.mean() if len(pos_rows) > 0 else 0.0

            norm_group = group.copy()
            for col in conc_cols:
                norm_group[col] = (
                    (group[col] - pos_mean) / (dmso_mean - pos_mean + 1e-9) * 100
                ).clip(-10, 120)
            normalized_frames.append(norm_group)

        self.normalized = pd.concat(normalized_frames, ignore_index=True)
        logger.info("Normalization complete.")
        return self.normalized

    # ── Curve fitting ─────────────────────────────────────────
    def fit_all_curves(
        self,
        outfile: str | Path = "results/functional/curve_parameters.csv",
    ) -> pd.DataFrame:
        """
        Fit 4PL dose-response curves for each drug × sample combination.

        Returns
        -------
        DataFrame with IC50, AUC, Emax, slope, R² for all combinations.
        """
        if self.normalized is None:
            self.normalize()

        df = self.normalized.copy()
        conc_cols = [c for c in df.columns if c.startswith("conc_")]
        records = []

        # Average across replicates if present
        id_cols = ["sample_id", "pdo_line", "passage", "condition", "drug"]
        id_cols = [c for c in id_cols if c in df.columns]
        df_agg = df[id_cols + conc_cols].groupby(id_cols).mean().reset_index()

        for _, row in df_agg.iterrows():
            viab = row[conc_cols].values.astype(float)
            params = fit_4pl(self.concentrations, viab)
            records.append({**row[id_cols].to_dict(), **params})

        self.curve_params = pd.DataFrame(records)
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        self.curve_params.to_csv(outfile, index=False)
        success_rate = self.curve_params["fit_success"].mean()
        logger.info(
            f"Curve fitting complete: {len(self.curve_params)} fits, "
            f"{success_rate:.1%} successful."
        )
        return self.curve_params

    # ── Functional drift ──────────────────────────────────────
    def compute_drift(
        self,
        baseline_passage: int = 2,
        outfile: str | Path = "results/functional/functional_drift.csv",
    ) -> pd.DataFrame:
        """
        Compute functional drift relative to the baseline passage.

        For each drug × PDO line × passage, compute:
        - IC50 fold-change vs. baseline
        - AUC delta vs. baseline
        - Flag if |IC50 FC| >= threshold

        Returns
        -------
        Tidy DataFrame of drift metrics.
        """
        if self.curve_params is None:
            self.fit_all_curves()

        df = self.curve_params[self.curve_params["fit_success"]].copy()
        baseline = df[df["passage"] == baseline_passage][
            ["pdo_line", "drug", "ic50_nM", "auc"]
        ].rename(columns={"ic50_nM": "ic50_baseline", "auc": "auc_baseline"})

        merged = df.merge(baseline, on=["pdo_line", "drug"], how="left")
        merged["ic50_fold_change"] = merged["ic50_nM"] / (merged["ic50_baseline"] + 1e-9)
        merged["auc_delta"] = merged["auc"] - merged["auc_baseline"]
        merged["ic50_flag"] = merged["ic50_fold_change"].abs() >= self.ic50_threshold
        merged["direction"] = np.where(
            merged["ic50_fold_change"] > 1, "resistant", "sensitized"
        )

        self.drift_summary = merged
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(outfile, index=False)
        n_flagged = merged["ic50_flag"].sum()
        logger.info(
            f"Functional drift computed: {n_flagged} drug × sample pairs "
            f"exceed IC50 shift threshold (≥{self.ic50_threshold}x)"
        )
        return merged

    # ── Freeze-thaw recovery analysis ────────────────────────
    def freeze_thaw_recovery(
        self,
        timepoints: list[str] = None,
        outfile: str | Path = "results/functional/freeze_thaw_recovery.csv",
    ) -> pd.DataFrame:
        """
        Assess recovery of drug sensitivity after freeze-thaw cycling.
        Compares fresh vs. post-thaw at multiple timepoints.

        Returns
        -------
        DataFrame with per-drug, per-timepoint IC50 recovery metrics.
        """
        if self.curve_params is None:
            self.fit_all_curves()
        if timepoints is None:
            timepoints = ["24h", "48h", "72h", "7d"]

        df = self.curve_params.copy()
        fresh = df[df["condition"] == "fresh"][["pdo_line", "drug", "ic50_nM", "auc"]]\
            .rename(columns={"ic50_nM": "ic50_fresh", "auc": "auc_fresh"})

        records = []
        for tp in timepoints:
            thaw = df[df["condition"] == f"postthaw_{tp}"][[
                "pdo_line", "drug", "ic50_nM", "auc"
            ]].rename(columns={"ic50_nM": f"ic50_{tp}", "auc": f"auc_{tp}"})
            records.append(fresh.merge(thaw, on=["pdo_line", "drug"], how="outer"))

        result = records[0]
        for r in records[1:]:
            result = result.merge(r, on=["pdo_line", "drug"], how="outer")

        # Recovery ratio: 1.0 = fully recovered to fresh IC50
        for tp in timepoints:
            col = f"ic50_{tp}"
            if col in result.columns:
                result[f"recovery_ratio_{tp}"] = 1 - np.abs(
                    result[col] - result["ic50_fresh"]
                ) / (result["ic50_fresh"].abs() + 1e-9)

        result.to_csv(outfile, index=False)
        logger.info(f"Freeze-thaw recovery analysis saved to {outfile}")
        return result

    # ── Visualization ─────────────────────────────────────────
    def plot_heatmap(
        self,
        metric: str = "ic50_fold_change",
        outfile: str | Path = "results/figures/drug_heatmap.pdf",
        figsize: tuple = (20, 12),
    ) -> None:
        """
        Heatmap of IC50 fold-changes across drugs × passage × PDO line.
        """
        if self.drift_summary is None:
            self.compute_drift()

        df = self.drift_summary[self.drift_summary["fit_success"]].copy()
        df["log2_fc"] = np.log2(df["ic50_fold_change"].clip(0.01, 100))
        pivot = df.pivot_table(
            index="drug", columns=["pdo_line", "passage"], values="log2_fc", aggfunc="mean"
        )

        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(
            pivot, ax=ax,
            cmap=sns.diverging_palette(220, 20, as_cmap=True),
            center=0, vmin=-3, vmax=3,
            linewidths=0.3, linecolor="white",
            cbar_kws={"label": "log₂(IC50 fold-change vs P2)", "shrink": 0.6},
            xticklabels=True, yticklabels=True,
        )
        ax.set_title(
            "Drug Sensitivity Drift — log₂(IC50 FC vs P2 Baseline)",
            fontsize=13, fontweight="bold", pad=10
        )
        ax.set_xlabel("PDO Line / Passage", fontsize=10)
        ax.set_ylabel("Drug", fontsize=10)
        plt.xticks(rotation=45, ha="right", fontsize=7)
        plt.yticks(fontsize=8)
        plt.tight_layout()
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Drug sensitivity heatmap saved to {outfile}")

    def plot_dose_response_grid(
        self,
        pdo_line: str,
        passages: list[int] = None,
        n_drugs: int = 12,
        outfile: str | Path = None,
    ) -> None:
        """
        Grid of dose-response curves for a single PDO line across passages.
        """
        if self.normalized is None:
            raise RuntimeError("Run normalize() first.")
        if passages is None:
            passages = sorted(self.normalized["passage"].unique())

        df = self.normalized[self.normalized["pdo_line"] == pdo_line].copy()
        drugs = df["drug"].value_counts().head(n_drugs).index.tolist()
        conc_cols = [c for c in df.columns if c.startswith("conc_")]

        n_cols = 4
        n_rows = int(np.ceil(n_drugs / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 3))
        axes = axes.flatten()
        palette = dict(zip(passages, sns.color_palette("viridis", n_colors=len(passages))))

        for i, drug in enumerate(drugs):
            ax = axes[i]
            drug_df = df[df["drug"] == drug]

            for passage in passages:
                p_df = drug_df[drug_df["passage"] == passage]
                if p_df.empty:
                    continue
                viab = p_df[conc_cols].mean().values
                ax.semilogx(
                    self.concentrations, viab,
                    "o-", color=palette[passage],
                    markersize=4, linewidth=1.5,
                    label=f"P{passage}"
                )

            ax.axhline(50, color="gray", linestyle=":", alpha=0.5, linewidth=0.8)
            ax.set_title(drug, fontsize=8, fontweight="bold")
            ax.set_xlabel("Conc (nM)", fontsize=7)
            ax.set_ylabel("Viability (%)", fontsize=7)
            ax.set_ylim(-10, 120)
            ax.tick_params(labelsize=6)
            if i == 0:
                ax.legend(fontsize=6, ncol=2)

        # Hide unused axes
        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        fig.suptitle(
            f"Dose-Response Curves — {pdo_line} Across Passages",
            fontsize=13, fontweight="bold", y=1.01
        )
        plt.tight_layout()

        if outfile is None:
            outfile = f"results/figures/drc_{pdo_line}.pdf"
        Path(outfile).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outfile, dpi=200, bbox_inches="tight")
        plt.close()
        logger.info(f"Dose-response grid saved to {outfile}")
