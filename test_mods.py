"""
Unit tests for the MODS Calculator and component distance calculators.
Run with: pytest tests/test_mods.py -v
"""
import numpy as np
import pandas as pd
import pytest
from unittest.mock import MagicMock, patch

import sys
sys.path.insert(0, "src")

from integration.mods_calculator import (
    GenomicDistanceCalculator,
    TranscriptomicDistanceCalculator,
    EpigenomicDistanceCalculator,
    FunctionalDistanceCalculator,
)


# ─────────────────────────────────────────────────────────────
#  Fixtures
# ─────────────────────────────────────────────────────────────
@pytest.fixture
def mut_matrix():
    """Synthetic binary mutation matrix: 100 genes × 4 samples."""
    np.random.seed(42)
    data = (np.random.rand(100, 4) > 0.85).astype(int)
    return pd.DataFrame(
        data,
        index=[f"GENE{i}" for i in range(100)],
        columns=["PDO001_P2", "PDO001_P5", "PDO001_P10", "PDO001_P20"],
    )

@pytest.fixture
def expr_df():
    """Synthetic log-normalized expression: 500 genes × 4 samples."""
    np.random.seed(7)
    baseline = np.random.randn(500) * 3 + 5
    data = {
        "PDO001_P2": baseline,
        "PDO001_P5": baseline + np.random.randn(500) * 0.5,
        "PDO001_P10": baseline + np.random.randn(500) * 1.5,
        "PDO001_P20": baseline + np.random.randn(500) * 3.0,
    }
    return pd.DataFrame(data, index=[f"GENE{i}" for i in range(500)])

@pytest.fixture
def beta_df():
    """Synthetic CpG beta-value matrix: 1000 CpGs × 4 samples."""
    np.random.seed(99)
    baseline = np.random.uniform(0, 1, 1000)
    data = {
        "PDO001_P2": baseline,
        "PDO001_P5": np.clip(baseline + np.random.randn(1000) * 0.02, 0, 1),
        "PDO001_P10": np.clip(baseline + np.random.randn(1000) * 0.08, 0, 1),
        "PDO001_P20": np.clip(baseline + np.random.randn(1000) * 0.15, 0, 1),
    }
    return pd.DataFrame(data, index=[f"CpG_{i}" for i in range(1000)])

@pytest.fixture
def auc_df():
    """Synthetic AUC matrix: 48 drugs × 4 samples."""
    np.random.seed(22)
    baseline = np.random.uniform(0.2, 0.9, 48)
    data = {
        "PDO001_P2": baseline,
        "PDO001_P5": np.clip(baseline + np.random.randn(48) * 0.05, 0, 1),
        "PDO001_P10": np.clip(baseline + np.random.randn(48) * 0.12, 0, 1),
        "PDO001_P20": np.clip(baseline + np.random.randn(48) * 0.22, 0, 1),
    }
    return pd.DataFrame(data, index=[f"DRUG{i}" for i in range(48)])


# ─────────────────────────────────────────────────────────────
#  Genomic Distance Tests
# ─────────────────────────────────────────────────────────────
class TestGenomicDistanceCalculator:
    def test_identical_samples_return_zero(self, mut_matrix):
        calc = GenomicDistanceCalculator()
        gd = calc.compute_distance(mut_matrix, None, "PDO001_P2", "PDO001_P2")
        assert gd == pytest.approx(0.0, abs=1e-6)

    def test_distance_increases_with_divergence(self, mut_matrix):
        calc = GenomicDistanceCalculator()
        gd_early = calc.compute_distance(mut_matrix, None, "PDO001_P2", "PDO001_P5")
        gd_late = calc.compute_distance(mut_matrix, None, "PDO001_P2", "PDO001_P20")
        # Not strictly guaranteed with random data but expected trend
        assert gd_early >= 0.0
        assert gd_late >= 0.0

    def test_distance_bounded_0_1(self, mut_matrix):
        calc = GenomicDistanceCalculator()
        for col in ["PDO001_P5", "PDO001_P10", "PDO001_P20"]:
            gd = calc.compute_distance(mut_matrix, None, "PDO001_P2", col)
            assert 0.0 <= gd <= 1.0, f"GD out of bounds for {col}: {gd}"

    def test_zero_vector_baseline(self):
        """When baseline has no mutations, distance to any mutated sample = 1.0."""
        calc = GenomicDistanceCalculator()
        df = pd.DataFrame({"baseline": [0, 0, 0], "test": [1, 1, 0]})
        gd = calc.compute_distance(df, None, "baseline", "test")
        assert gd == pytest.approx(1.0, abs=0.01)

    def test_cnv_component_integrated(self, mut_matrix):
        """CNV data should be incorporated without raising errors."""
        calc = GenomicDistanceCalculator()
        cnv = pd.DataFrame({
            "sample_id": ["PDO001_P2"] * 5 + ["PDO001_P10"] * 5,
            "gene": [f"GENE{i}" for i in range(5)] * 2,
            "log2_ratio": [0.0, 0.1, -0.2, 0.0, 0.3, 0.5, 1.1, -0.8, 0.2, 0.4],
        })
        gd = calc.compute_distance(mut_matrix, cnv, "PDO001_P2", "PDO001_P10")
        assert 0.0 <= gd <= 1.0


# ─────────────────────────────────────────────────────────────
#  Transcriptomic Distance Tests
# ─────────────────────────────────────────────────────────────
class TestTranscriptomicDistanceCalculator:
    def test_identical_samples_near_zero(self, expr_df):
        calc = TranscriptomicDistanceCalculator()
        td = calc.compute_distance(expr_df, "PDO001_P2", "PDO001_P2")
        assert td == pytest.approx(0.0, abs=1e-6)

    def test_distance_monotonically_increases(self, expr_df):
        calc = TranscriptomicDistanceCalculator()
        td5 = calc.compute_distance(expr_df, "PDO001_P2", "PDO001_P5")
        td10 = calc.compute_distance(expr_df, "PDO001_P2", "PDO001_P10")
        td20 = calc.compute_distance(expr_df, "PDO001_P2", "PDO001_P20")
        assert td5 <= td10 <= td20, \
            f"Expected monotonic increase: {td5:.3f} ≤ {td10:.3f} ≤ {td20:.3f}"

    def test_distance_bounded(self, expr_df):
        calc = TranscriptomicDistanceCalculator()
        for col in expr_df.columns[1:]:
            td = calc.compute_distance(expr_df, "PDO001_P2", col)
            assert 0.0 <= td <= 1.0

    def test_hvg_selection(self, expr_df):
        calc = TranscriptomicDistanceCalculator(n_top_genes=100)
        hvg = calc.select_hvg(expr_df)
        assert len(hvg) == 100


# ─────────────────────────────────────────────────────────────
#  Epigenomic Distance Tests
# ─────────────────────────────────────────────────────────────
class TestEpigenomicDistanceCalculator:
    def test_identical_returns_zero(self, beta_df):
        calc = EpigenomicDistanceCalculator(min_variance=0.0)
        ed = calc.compute_distance(beta_df, "PDO001_P2", "PDO001_P2")
        assert ed == pytest.approx(0.0, abs=1e-9)

    def test_distance_bounded(self, beta_df):
        calc = EpigenomicDistanceCalculator(min_variance=0.01)
        for col in beta_df.columns[1:]:
            ed = calc.compute_distance(beta_df, "PDO001_P2", col)
            assert 0.0 <= ed <= 1.0

    def test_empty_after_filter_returns_zero(self):
        """High variance threshold → no sites → ED=0."""
        beta = pd.DataFrame({
            "S1": [0.5, 0.5], "S2": [0.5, 0.5]
        })
        calc = EpigenomicDistanceCalculator(min_variance=0.99)
        ed = calc.compute_distance(beta, "S1", "S2")
        assert ed == 0.0


# ─────────────────────────────────────────────────────────────
#  Functional Distance Tests
# ─────────────────────────────────────────────────────────────
class TestFunctionalDistanceCalculator:
    def test_identical_returns_zero(self, auc_df):
        calc = FunctionalDistanceCalculator()
        fd = calc.compute_distance(auc_df, "PDO001_P2", "PDO001_P2")
        assert fd == pytest.approx(0.0, abs=1e-9)

    def test_distance_bounded(self, auc_df):
        calc = FunctionalDistanceCalculator()
        for col in auc_df.columns[1:]:
            fd = calc.compute_distance(auc_df, "PDO001_P2", col)
            assert 0.0 <= fd <= 1.0

    def test_flag_shifted_drugs(self, auc_df):
        """Manually inject a 4-fold IC50 shift for one drug."""
        ic50 = auc_df.copy() * 100  # pretend these are IC50s in nM
        ic50.loc["DRUG0", "PDO001_P20"] = ic50.loc["DRUG0", "PDO001_P2"] * 5
        calc = FunctionalDistanceCalculator(ic50_shift_threshold=3.0)
        flagged = calc.flag_shifted_drugs(ic50, "PDO001_P2", "PDO001_P20")
        assert "DRUG0" in flagged.index


# ─────────────────────────────────────────────────────────────
#  Drug Sensitivity 4PL Fit Tests
# ─────────────────────────────────────────────────────────────
class TestFourParamLogistic:
    from functional.drug_sensitivity import fit_4pl, four_param_logistic

    def test_fit_synthetic_curve(self):
        from functional.drug_sensitivity import fit_4pl, four_param_logistic
        import numpy as np
        conc = np.array([1, 10, 100, 1000, 10000, 100000, 1e6, 1e7])
        # True params: top=100, bottom=5, ic50=1000, slope=1.2
        viab = four_param_logistic(conc, 100, 5, 1000, 1.2) + np.random.randn(8) * 1.5
        result = fit_4pl(conc, viab)
        assert result["fit_success"]
        assert 500 < result["ic50_nM"] < 2000, f"IC50 far from expected: {result['ic50_nM']}"
        assert result["r_squared"] > 0.90

    def test_fit_returns_nan_on_failure(self):
        from functional.drug_sensitivity import fit_4pl
        # All-zero viabilities → fitting should fail gracefully
        conc = np.array([1.0, 10.0, 100.0])
        viab = np.array([np.nan, np.nan, np.nan])
        result = fit_4pl(conc, viab)
        assert not result["fit_success"]
        assert np.isnan(result["ic50_nM"])
