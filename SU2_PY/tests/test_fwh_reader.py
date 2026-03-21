#!/usr/bin/env python3
"""Tests for SU2.io.fwh — FWH binary surface data reader.

Run with:
    cd /home/riddhi/SU2
    python3 -m pytest SU2_PY/tests/test_fwh_reader.py -v
"""

import struct
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from SU2.io.fwh import FWHData

REPO_ROOT = Path(__file__).parent.parent.parent
REFERENCE_FILE = REPO_ROOT / "TestCases/unsteady/square_cylinder/fwh_bin.dat"


def make_synthetic_fwh(path, n_t=10, n_pts=5, start_us=1000):
    rng = np.random.default_rng(42)
    pressure = (101325.0 + rng.standard_normal((n_t, n_pts)) * 100).astype("<f4")
    with open(path, "wb") as fh:
        fh.write(struct.pack("<3I", start_us, n_t, n_pts))
        fh.write(pressure.tobytes())
    return pressure


class TestFWHDataSynthetic:
    def test_roundtrip_shape(self, tmp_path):
        make_synthetic_fwh(tmp_path / "fwh_bin.dat", n_t=10, n_pts=5)
        fwh = FWHData.from_file(tmp_path / "fwh_bin.dat", time_step=0.0015)
        assert fwh.pressure.shape == (10, 5)
        assert fwh.n_timesteps == 10
        assert fwh.n_points == 5

    def test_roundtrip_values(self, tmp_path):
        p = make_synthetic_fwh(tmp_path / "fwh_bin.dat", n_t=10, n_pts=5)
        fwh = FWHData.from_file(tmp_path / "fwh_bin.dat", time_step=0.0015)
        np.testing.assert_allclose(fwh.pressure, p, rtol=1e-6)

    def test_time_axis_with_dt(self, tmp_path):
        make_synthetic_fwh(tmp_path / "fwh_bin.dat", n_t=10, n_pts=5, start_us=500000)
        fwh = FWHData.from_file(tmp_path / "fwh_bin.dat", time_step=0.0015)
        assert fwh.time[0] == pytest.approx(0.5, rel=1e-6)
        assert fwh.time[-1] == pytest.approx(0.5 + 9 * 0.0015, rel=1e-6)
        assert len(fwh.time) == 10

    def test_time_axis_no_dt(self, tmp_path):
        make_synthetic_fwh(tmp_path / "fwh_bin.dat", n_t=8, n_pts=3)
        fwh = FWHData.from_file(tmp_path / "fwh_bin.dat")
        np.testing.assert_array_equal(fwh.time, np.arange(8, dtype=float))

    def test_pressure_perturbation_zero_mean(self, tmp_path):
        make_synthetic_fwh(tmp_path / "fwh_bin.dat", n_t=10, n_pts=5)
        fwh = FWHData.from_file(tmp_path / "fwh_bin.dat")
        p_prime = fwh.pressure_perturbation()
        np.testing.assert_allclose(
            p_prime.mean(axis=0), 0.0, atol=0.01
        )  # float32 quantisation

    def test_pressure_perturbation_with_ref(self, tmp_path):
        make_synthetic_fwh(tmp_path / "fwh_bin.dat", n_t=10, n_pts=5)
        fwh = FWHData.from_file(tmp_path / "fwh_bin.dat")
        p_prime = fwh.pressure_perturbation(p_ref=101325.0)
        np.testing.assert_allclose(p_prime, fwh.pressure - 101325.0, rtol=1e-6)

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError, match="not found"):
            FWHData.from_file("/nonexistent/path/fwh_bin.dat")

    def test_truncated_file_raises(self, tmp_path):
        p = tmp_path / "fwh_bin.dat"
        p.write_bytes(struct.pack("<3I", 0, 10, 5) + b"\x00" * 100)
        with pytest.raises(ValueError, match="size mismatch"):
            FWHData.from_file(p)

    def test_repr(self, tmp_path):
        make_synthetic_fwh(tmp_path / "fwh_bin.dat")
        fwh = FWHData.from_file(tmp_path / "fwh_bin.dat", time_step=0.0015)
        r = repr(fwh)
        assert "FWHData" in r
        assert "n_timesteps=10" in r

    def test_to_dict_keys(self, tmp_path):
        make_synthetic_fwh(tmp_path / "fwh_bin.dat")
        fwh = FWHData.from_file(tmp_path / "fwh_bin.dat")
        d = fwh.to_dict()
        for key in ("n_timesteps", "n_points", "time", "pressure_mean", "pressure_rms"):
            assert key in d

    def test_time_mean_pressure_shape(self, tmp_path):
        make_synthetic_fwh(tmp_path / "fwh_bin.dat", n_t=10, n_pts=5)
        fwh = FWHData.from_file(tmp_path / "fwh_bin.dat")
        assert fwh.time_mean_pressure().shape == (5,)

    def test_pressure_rms_shape(self, tmp_path):
        make_synthetic_fwh(tmp_path / "fwh_bin.dat", n_t=10, n_pts=5)
        fwh = FWHData.from_file(tmp_path / "fwh_bin.dat")
        assert fwh.pressure_rms().shape == (5,)


@pytest.mark.skipif(
    not REFERENCE_FILE.exists(),
    reason="Reference fwh_bin.dat not found — run from repo root",
)
class TestFWHDataReference:
    @pytest.fixture(scope="class")
    def fwh(self):
        return FWHData.from_file(REFERENCE_FILE, time_step=0.0015)

    def test_shape(self, fwh):
        assert fwh.n_timesteps == 50
        assert fwh.n_points == 200
        assert fwh.pressure.shape == (50, 200)

    def test_pressure_range(self, fwh):
        assert fwh.pressure.min() > 800.0
        assert fwh.pressure.max() < 1200.0

    def test_unsteady_variation(self, fwh):
        assert fwh.pressure.std(axis=0).max() > 10.0

    def test_time_step_spacing(self, fwh):
        diffs = np.diff(fwh.time)
        np.testing.assert_allclose(diffs, 0.0015, rtol=1e-6)

    def test_pressure_rms_max(self, fwh):
        assert fwh.pressure_rms().max() == pytest.approx(74.22, abs=0.5)
