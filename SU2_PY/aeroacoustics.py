#!/usr/bin/env python3

## \file aeroacoustics.py
#  \brief Spectral analysis of unsteady pressure histories.
#  \version 8.5.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

"""Convert SU2 pressure histories into acoustic spectra and SPL metrics."""

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np


WINDOWS = {
    "hann": np.hanning,
    "hamming": np.hamming,
    "blackman": np.blackman,
    "rectangular": np.ones,
}


def read_history(filename):
    """Read a numeric SU2 CSV history file and return its columns."""
    with open(filename, newline="", encoding="utf-8-sig") as history_file:
        rows = (
            row
            for row in history_file
            if row.strip() and not row.lstrip().startswith(("%", "#"))
        )
        reader = csv.reader(rows)
        try:
            headers = [field.strip().strip('"') for field in next(reader)]
        except StopIteration as error:
            raise ValueError("history file is empty") from error

        if not headers or any(not field for field in headers):
            raise ValueError("history file contains an empty column name")
        if len(set(headers)) != len(headers):
            raise ValueError("history file contains duplicate column names")

        columns = {header: [] for header in headers}
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(headers):
                raise ValueError(
                    "line {} has {} values; expected {}".format(
                        line_number, len(row), len(headers)
                    )
                )
            try:
                for header, value in zip(headers, row):
                    columns[header].append(float(value))
            except ValueError as error:
                raise ValueError(
                    "line {} contains a non-numeric value".format(line_number)
                ) from error

    if not next(iter(columns.values())):
        raise ValueError("history file has a header but no data")
    return {name: np.asarray(values, dtype=float) for name, values in columns.items()}


def resolve_column(columns, requested):
    """Resolve a column exactly, or uniquely without case sensitivity."""
    if requested in columns:
        return requested
    matches = [name for name in columns if name.casefold() == requested.casefold()]
    if len(matches) == 1:
        return matches[0]
    raise ValueError(
        "column {!r} was not found; available columns: {}".format(
            requested, ", ".join(columns)
        )
    )


def sampling_rate(time, relative_tolerance=1.0e-5):
    """Return the sample rate after verifying monotonically uniform time data."""
    if time.size < 2:
        raise ValueError("at least two time samples are required")
    if not np.all(np.isfinite(time)):
        raise ValueError("time values must be finite")

    intervals = np.diff(time)
    if np.any(intervals <= 0.0):
        raise ValueError("time values must be strictly increasing")
    mean_interval = float(np.mean(intervals))
    maximum_error = float(np.max(np.abs(intervals - mean_interval)))
    if maximum_error > relative_tolerance * mean_interval:
        message = (
            "time samples are not uniform " "(maximum relative interval error {:.3g})"
        )
        raise ValueError(message.format(maximum_error / mean_interval))
    return 1.0 / mean_interval


def _detrend(signal, mode):
    if mode == "none":
        return signal.copy()
    if mode == "mean":
        return signal - np.mean(signal)
    positions = np.arange(signal.size, dtype=float)
    slope, intercept = np.polyfit(positions, signal, 1)
    return signal - (slope * positions + intercept)


def welch_psd(
    pressure,
    sample_rate,
    segment_length=None,
    overlap=0.5,
    window="hann",
    detrend="mean",
):
    """Estimate a one-sided pressure PSD using Welch's averaged periodogram."""
    pressure = np.asarray(pressure, dtype=float)
    if pressure.ndim != 1:
        raise ValueError("pressure data must be one-dimensional")
    if pressure.size < 2:
        raise ValueError("at least two pressure samples are required")
    if not np.all(np.isfinite(pressure)):
        raise ValueError("pressure values must be finite")
    if not math.isfinite(sample_rate) or sample_rate <= 0.0:
        raise ValueError("sample rate must be positive and finite")
    if not 0.0 <= overlap < 1.0:
        raise ValueError("overlap must be in the range [0, 1)")
    if window not in WINDOWS:
        raise ValueError("unsupported window {!r}".format(window))
    if detrend not in ("none", "mean", "linear"):
        raise ValueError("unsupported detrending mode {!r}".format(detrend))

    if segment_length is None:
        segment_length = min(1024, pressure.size)
    if segment_length < 2 or segment_length > pressure.size:
        raise ValueError("segment length must be between 2 and the number of samples")

    step = max(1, int(round(segment_length * (1.0 - overlap))))
    starts = range(0, pressure.size - segment_length + 1, step)
    weights = WINDOWS[window](segment_length)
    window_energy = float(np.sum(weights * weights))
    if window_energy <= 0.0:
        raise ValueError(
            "segment length {} is too short for the {} window".format(
                segment_length, window
            )
        )
    spectrum = np.zeros(segment_length // 2 + 1, dtype=float)
    segment_count = 0

    for start in starts:
        segment = _detrend(pressure[start : start + segment_length], detrend)
        transform = np.fft.rfft(segment * weights)
        periodogram = np.abs(transform) ** 2 / (sample_rate * window_energy)
        if segment_length % 2 == 0:
            periodogram[1:-1] *= 2.0
        else:
            periodogram[1:] *= 2.0
        spectrum += periodogram
        segment_count += 1

    spectrum /= segment_count
    frequencies = np.fft.rfftfreq(segment_length, d=1.0 / sample_rate)
    return frequencies, spectrum, segment_count


def acoustic_metrics(frequencies, psd, reference_pressure=20.0e-6):
    """Calculate PSD level, bin SPL, OASPL, and the dominant frequency."""
    if reference_pressure <= 0.0 or not math.isfinite(reference_pressure):
        raise ValueError("reference pressure must be positive and finite")
    if frequencies.size < 2:
        raise ValueError("an acoustic spectrum needs at least two frequency bins")

    frequency_step = float(frequencies[1] - frequencies[0])
    tiny = np.finfo(float).tiny
    psd_level = 10.0 * np.log10(np.maximum(psd, tiny) / reference_pressure**2)
    narrowband_spl = 10.0 * np.log10(
        np.maximum(psd * frequency_step, tiny) / reference_pressure**2
    )
    oaspl = 10.0 * math.log10(
        max(float(np.sum(psd) * frequency_step), tiny) / reference_pressure**2
    )
    dominant_index = 1 + int(np.argmax(psd[1:]))
    return {
        "psd_level_db_per_hz": psd_level,
        "narrowband_spl_db": narrowband_spl,
        "oaspl_db": oaspl,
        "dominant_frequency_hz": float(frequencies[dominant_index]),
    }


def _safe_column_name(name):
    normalized = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")
    return normalized or "pressure"


def write_spectrum(filename, frequencies, spectra):
    """Write frequency-domain results for one or more pressure signals."""
    fields = ["Frequency_Hz"]
    used_prefixes = set()
    for name in spectra:
        base = _safe_column_name(name)
        prefix = base
        suffix = 2
        while prefix in used_prefixes:
            prefix = "{}_{}".format(base, suffix)
            suffix += 1
        used_prefixes.add(prefix)
        fields.extend(
            [
                prefix + "_PSD_Pa2_per_Hz",
                prefix + "_PSD_level_dB_per_Hz",
                prefix + "_SPL_dB",
            ]
        )

    with open(filename, "w", newline="", encoding="utf-8") as spectrum_file:
        writer = csv.writer(spectrum_file)
        writer.writerow(fields)
        for index, frequency in enumerate(frequencies):
            row = [frequency]
            for result in spectra.values():
                row.extend(
                    [
                        result["psd"][index],
                        result["psd_level_db_per_hz"][index],
                        result["narrowband_spl_db"][index],
                    ]
                )
            writer.writerow(row)


def analyze(args):
    """Run the command-line analysis and return its summary."""
    columns = read_history(args.input)
    pressure_names = [resolve_column(columns, name) for name in args.pressure]
    if len(set(pressure_names)) != len(pressure_names):
        raise ValueError("each pressure column may only be requested once")
    if not math.isfinite(args.pressure_scale) or args.pressure_scale <= 0.0:
        raise ValueError("pressure scale must be positive and finite")
    if args.skip_samples < 0:
        raise ValueError("skip samples must not be negative")
    sample_count = len(next(iter(columns.values())))
    if args.skip_samples > sample_count - 2:
        raise ValueError("skip samples must leave at least two samples to analyze")

    if args.sample_rate is not None:
        sample_rate_hz = args.sample_rate
    else:
        time_name = resolve_column(columns, args.time)
        sample_rate_hz = sampling_rate(columns[time_name][args.skip_samples :])

    spectra = {}
    summary = {
        "input": str(args.input),
        "sample_rate_hz": sample_rate_hz,
        "reference_pressure_pa": args.reference_pressure,
        "samples": sample_count - args.skip_samples,
        "signals": {},
    }
    frequencies = None
    for name in pressure_names:
        dimensional_pressure = columns[name][args.skip_samples :] * args.pressure_scale
        frequencies, psd, segment_count = welch_psd(
            dimensional_pressure,
            sample_rate_hz,
            segment_length=args.segment_length,
            overlap=args.overlap,
            window=args.window,
            detrend=args.detrend,
        )
        metrics = acoustic_metrics(frequencies, psd, args.reference_pressure)
        spectra[name] = {"psd": psd, **metrics}
        summary["signals"][name] = {
            "oaspl_db": metrics["oaspl_db"],
            "dominant_frequency_hz": metrics["dominant_frequency_hz"],
            "segments": segment_count,
        }

    summary["frequency_resolution_hz"] = float(frequencies[1] - frequencies[0])
    summary["nyquist_frequency_hz"] = float(frequencies[-1])
    write_spectrum(args.output, frequencies, spectra)
    if args.summary:
        with open(args.summary, "w", encoding="utf-8") as summary_file:
            json.dump(summary, summary_file, indent=2)
            summary_file.write("\n")
    return summary


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Estimate acoustic spectra from pressure columns in an SU2 history CSV."
        )
    )
    parser.add_argument("input", type=Path, help="SU2 history CSV file")
    parser.add_argument(
        "-p",
        "--pressure",
        action="append",
        required=True,
        metavar="COLUMN",
        help="pressure column to analyze; repeat for multiple probes",
    )
    sampling = parser.add_mutually_exclusive_group()
    sampling.add_argument(
        "-t",
        "--time",
        default="Cur_Time",
        metavar="COLUMN",
        help="physical-time column (default: Cur_Time)",
    )
    sampling.add_argument(
        "--sample-rate",
        type=float,
        metavar="HZ",
        help="sampling rate when the file has no physical-time column",
    )
    parser.add_argument(
        "--pressure-scale",
        type=float,
        default=1.0,
        metavar="PA_PER_UNIT",
        help="multiply input pressures by this dimensionalization factor",
    )
    parser.add_argument(
        "--reference-pressure",
        type=float,
        default=20.0e-6,
        metavar="PA",
        help="SPL reference pressure (default: 20e-6 Pa)",
    )
    parser.add_argument(
        "-n",
        "--segment-length",
        type=int,
        metavar="SAMPLES",
        help="Welch segment length (default: min(1024, number of samples))",
    )
    parser.add_argument(
        "--overlap",
        type=float,
        default=0.5,
        metavar="FRACTION",
        help="Welch segment overlap as a fraction (default: 0.5)",
    )
    parser.add_argument(
        "--window",
        choices=tuple(WINDOWS),
        default="hann",
        help="spectral window (default: hann)",
    )
    parser.add_argument(
        "--detrend",
        choices=("none", "mean", "linear"),
        default="mean",
        help="detrending applied to each segment (default: mean)",
    )
    parser.add_argument(
        "--skip-samples",
        type=int,
        default=0,
        metavar="COUNT",
        help="discard initial transient samples before analysis (default: 0)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("acoustic_spectrum.csv"),
        help="output spectrum CSV (default: acoustic_spectrum.csv)",
    )
    parser.add_argument(
        "--summary", type=Path, help="optional JSON file for scalar acoustic metrics"
    )
    return parser


def main():
    args = build_parser().parse_args()
    try:
        summary = analyze(args)
    except (OSError, ValueError) as error:
        raise SystemExit("aeroacoustics: error: {}".format(error)) from error

    print("Sample rate: {:.8g} Hz".format(summary["sample_rate_hz"]))
    for name, result in summary["signals"].items():
        print(
            "{}: OASPL = {:.3f} dB, dominant frequency = {:.8g} Hz".format(
                name, result["oaspl_db"], result["dominant_frequency_hz"]
            )
        )
    print("Spectrum written to {}".format(args.output))


if __name__ == "__main__":
    main()
