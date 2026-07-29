#!/usr/bin/env python3

import csv
import importlib.util
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np


MODULE_PATH = Path(__file__).parents[1] / "aeroacoustics.py"
SPEC = importlib.util.spec_from_file_location("aeroacoustics", MODULE_PATH)
AEROACOUSTICS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AEROACOUSTICS)


class AeroacousticsTests(unittest.TestCase):
    def test_welch_finds_sine_frequency_and_level(self):
        sample_rate = 4096.0
        frequency = 256.0
        rms_pressure = 1.0
        time = np.arange(4096) / sample_rate
        pressure = (
            math.sqrt(2.0) * rms_pressure * np.sin(2.0 * np.pi * frequency * time)
        )

        frequencies, psd, segments = AEROACOUSTICS.welch_psd(
            pressure,
            sample_rate,
            segment_length=1024,
            overlap=0.5,
            window="hann",
        )
        metrics = AEROACOUSTICS.acoustic_metrics(frequencies, psd)

        self.assertEqual(segments, 7)
        self.assertEqual(metrics["dominant_frequency_hz"], frequency)
        expected_oaspl = 20.0 * math.log10(rms_pressure / 20.0e-6)
        self.assertAlmostEqual(metrics["oaspl_db"], expected_oaspl, places=3)

    def test_sampling_rate_rejects_nonuniform_time(self):
        with self.assertRaisesRegex(ValueError, "not uniform"):
            AEROACOUSTICS.sampling_rate(np.array([0.0, 0.1, 0.21, 0.3]))

    def test_read_history_accepts_su2_quoted_headers_and_comments(self):
        with tempfile.TemporaryDirectory() as directory:
            history = Path(directory) / "history.csv"
            history.write_text(
                '% SU2 history\n"Time","probe 1"\n0.0,101325.0\n0.1,101326.0\n',
                encoding="utf-8",
            )
            columns = AEROACOUSTICS.read_history(history)

        self.assertEqual(list(columns), ["Time", "probe 1"])
        np.testing.assert_allclose(columns["probe 1"], [101325.0, 101326.0])

    def test_write_spectrum_supports_multiple_signals(self):
        frequencies = np.array([0.0, 1.0])
        result = {
            "psd": np.array([1.0, 2.0]),
            "psd_level_db_per_hz": np.array([3.0, 4.0]),
            "narrowband_spl_db": np.array([5.0, 6.0]),
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "spectrum.csv"
            AEROACOUSTICS.write_spectrum(
                output, frequencies, {"probe 1": result, "probe/2": result}
            )
            with output.open(newline="", encoding="utf-8") as output_file:
                rows = list(csv.reader(output_file))

        self.assertEqual(len(rows), 3)
        self.assertIn("probe_1_PSD_Pa2_per_Hz", rows[0])
        self.assertIn("probe_2_SPL_dB", rows[0])

    def test_write_spectrum_disambiguates_normalized_names(self):
        frequencies = np.array([0.0, 1.0])
        result = {
            "psd": np.ones(2),
            "psd_level_db_per_hz": np.ones(2),
            "narrowband_spl_db": np.ones(2),
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "spectrum.csv"
            AEROACOUSTICS.write_spectrum(
                output, frequencies, {"probe 1": result, "probe/1": result}
            )
            header = output.read_text(encoding="utf-8").splitlines()[0]

        self.assertIn("probe_1_PSD_Pa2_per_Hz", header)
        self.assertIn("probe_1_2_PSD_Pa2_per_Hz", header)

    def test_analyze_writes_spectrum_and_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            history = directory / "history.csv"
            output = directory / "spectrum.csv"
            summary_file = directory / "summary.json"
            time = np.arange(256) / 128.0
            pressure = np.sin(2.0 * np.pi * 16.0 * time)
            with history.open("w", newline="", encoding="utf-8") as history_file:
                writer = csv.writer(history_file)
                writer.writerow(["Cur_Time", "mic"])
                writer.writerows(zip(time, pressure))

            summary = AEROACOUSTICS.analyze(
                SimpleNamespace(
                    input=history,
                    pressure=["mic"],
                    sample_rate=None,
                    time="Cur_Time",
                    pressure_scale=1.0,
                    reference_pressure=20.0e-6,
                    skip_samples=64,
                    segment_length=64,
                    overlap=0.5,
                    window="hann",
                    detrend="mean",
                    output=output,
                    summary=summary_file,
                )
            )

            self.assertTrue(output.is_file())
            self.assertTrue(summary_file.is_file())
            self.assertEqual(summary["samples"], 192)
            self.assertEqual(summary["frequency_resolution_hz"], 2.0)
            self.assertEqual(summary["nyquist_frequency_hz"], 64.0)
            self.assertEqual(summary["signals"]["mic"]["dominant_frequency_hz"], 16.0)


if __name__ == "__main__":
    unittest.main()
