#!/usr/bin/env python3
# \file fwh.py
# \brief Reader for SU2 FWH binary surface data files (fwh_bin.dat).
#
# The binary format written by the SU2 FWH surface writer:
#   Header  : 3 × uint32 (little-endian)
#               [0] start_time_us  — simulation time at first stored timestep (microseconds)
#               [1] n_timesteps    — number of timesteps in this file
#               [2] n_points       — number of surface points
#   Data    : n_timesteps × n_points × float32 (little-endian), row-major
#               pressure[t, i]  — static pressure at surface point i, timestep t  (Pa)
#
# Usage:
#   from SU2.io.fwh import FWHData
#   fwh = FWHData.from_file("fwh_bin.dat", time_step=0.0015)
#   print(fwh.pressure.shape)   # (n_timesteps, n_points)
#   print(fwh.time)             # physical time array in seconds
#
# \author  Riddhi (GSoC 2026 candidate)
# \version 8.4.0 "Harrier"

import struct
import numpy as np
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class FWHData:
    """Container for FWH binary surface pressure data read from fwh_bin.dat.

    Attributes
    ----------
    pressure : np.ndarray, shape (n_timesteps, n_points), float32
        Static pressure at each FWH surface point for each timestep (Pa).
    time : np.ndarray, shape (n_timesteps,), float64
        Physical simulation time corresponding to each timestep (seconds).
    n_timesteps : int
        Number of timesteps stored in the file.
    n_points : int
        Number of surface quadrature points on the FWH surface.
    start_time_us : int
        Raw start-time tag from file header (microseconds, uint32).
    source_file : str
        Path to the file this data was read from.
    """

    pressure: np.ndarray
    time: np.ndarray
    n_timesteps: int
    n_points: int
    start_time_us: int
    source_file: str = ""

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_file(cls, path: str, time_step: float = 0.0) -> "FWHData":
        """Read a SU2 FWH binary file.

        Parameters
        ----------
        path : str or Path
            Path to fwh_bin.dat (or equivalent).
        time_step : float, optional
            CFD time step in seconds (TIME_STEP in the .cfg file).
            Used to build the physical time array.  If 0.0, the time
            array is integer timestep indices cast to float.

        Returns
        -------
        FWHData
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"FWH binary file not found: {path}")

        with open(path, "rb") as fh:
            raw = fh.read()

        # --- header ---------------------------------------------------
        if len(raw) < 12:
            raise ValueError(
                f"File too small to contain a valid FWH header: {len(raw)} bytes"
            )

        start_time_us, n_t, n_pts = struct.unpack_from("<3I", raw, 0)
        expected_data = n_t * n_pts * 4  # float32
        actual_data = len(raw) - 12

        if actual_data != expected_data:
            raise ValueError(
                f"FWH binary size mismatch: header says {n_t}×{n_pts} "
                f"({expected_data} bytes) but file has {actual_data} data bytes. "
                f"File may be truncated or use a different layout."
            )

        # --- data block -----------------------------------------------
        pressure = np.frombuffer(raw, dtype="<f4", offset=12).reshape(n_t, n_pts).copy()

        # --- time axis ------------------------------------------------
        if time_step > 0.0:
            t0 = start_time_us * 1e-6  # convert μs → s
            time = t0 + np.arange(n_t) * time_step
        else:
            time = np.arange(n_t, dtype=np.float64)

        return cls(
            pressure=pressure,
            time=time,
            n_timesteps=int(n_t),
            n_points=int(n_pts),
            start_time_us=int(start_time_us),
            source_file=str(path),
        )

    # ------------------------------------------------------------------
    # Derived quantities for FWH post-processing
    # ------------------------------------------------------------------

    def pressure_perturbation(self, p_ref: Optional[float] = None) -> np.ndarray:
        """Return p' = p - p_mean (or p - p_ref if supplied).

        Parameters
        ----------
        p_ref : float, optional
            Reference (free-stream) pressure in Pa.  If None, the
            time-mean pressure at each point is used.

        Returns
        -------
        np.ndarray, shape (n_timesteps, n_points)
        """
        if p_ref is not None:
            return self.pressure - p_ref
        return self.pressure - self.pressure.mean(axis=0, keepdims=True)

    def time_mean_pressure(self) -> np.ndarray:
        """Time-averaged pressure at each surface point (Pa).

        Returns
        -------
        np.ndarray, shape (n_points,)
        """
        return self.pressure.mean(axis=0)

    def pressure_rms(self) -> np.ndarray:
        """RMS of pressure fluctuation p' at each surface point (Pa).

        Returns
        -------
        np.ndarray, shape (n_points,)
        """
        return self.pressure_perturbation().std(axis=0)

    # ------------------------------------------------------------------
    # I/O helpers
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """Return a plain dict suitable for JSON serialisation or DataFrame construction."""
        return {
            "source_file": self.source_file,
            "n_timesteps": self.n_timesteps,
            "n_points": self.n_points,
            "start_time_us": self.start_time_us,
            "time": self.time.tolist(),
            "pressure_mean": self.time_mean_pressure().tolist(),
            "pressure_rms": self.pressure_rms().tolist(),
        }

    def to_dataframe(self):
        """Return a pandas DataFrame with columns [time, point_0, point_1, ...].

        Requires pandas.  Each column ``point_i`` contains the pressure
        time-series at surface point i.

        Returns
        -------
        pandas.DataFrame
        """
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                "pandas is required for FWHData.to_dataframe(). "
                "Install it with:  pip install pandas"
            ) from exc

        cols = {"time": self.time}
        cols.update({f"point_{i}": self.pressure[:, i] for i in range(self.n_points)})
        return pd.DataFrame(cols)

    # ------------------------------------------------------------------
    # Dunder
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"FWHData(n_timesteps={self.n_timesteps}, n_points={self.n_points}, "
            f"t=[{self.time[0]:.4f}, {self.time[-1]:.4f}] s, "
            f"p=[{self.pressure.min():.1f}, {self.pressure.max():.1f}] Pa)"
        )
