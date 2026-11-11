# Aeroacoustic analysis of unsteady pressure

`aeroacoustics.py` converts pressure histories from an unsteady SU2
simulation into frequency-domain acoustic quantities. It is intended for
pressure probes from LES, DES, DDES, ZDES, or EDDES cases, but it can analyze
any uniformly sampled pressure history.

The tool estimates the one-sided power spectral density (PSD) with Welch's
method and reports:

- pressure PSD in Pa²/Hz;
- PSD level in dB/Hz;
- bin-integrated narrowband sound pressure level (SPL) in dB;
- overall sound pressure level (OASPL) in dB; and
- the dominant non-zero frequency.

Only NumPy is required.

## Recording pressure probes

SU2 custom outputs can record pressure at one or more points. For example:

```ini
CUSTOM_OUTPUTS= 'mic1 : Probe{PRESSURE}[1.0, 0.0, 0.0]; \
                 mic2 : Probe{PRESSURE}[0.0, 1.0, 0.0]'
HISTORY_OUTPUT= TIME_ITER, CUR_TIME, mic1, mic2
```

Use a fixed physical time step and write every desired acoustic sample. The
history column containing physical time is commonly named `Cur_Time`; inspect
the CSV header and pass its exact name with `--time`.

## Usage

Analyze two pressure columns that are already dimensional:

```bash
python SU2_PY/aeroacoustics.py history.csv \
  --pressure mic1 \
  --pressure mic2 \
  --skip-samples 10000 \
  --segment-length 2048 \
  --summary acoustic_summary.json \
  --output acoustic_spectrum.csv
```

If the file contains an iteration rather than physical time, provide the
sampling rate explicitly:

```bash
python SU2_PY/aeroacoustics.py history.csv \
  --sample-rate 50000 \
  --pressure mic1
```

Pressure must be dimensional in pascals for physical SPL. If an SU2 case uses
nondimensional pressure, apply its dimensionalization factor in pascals per
input unit:

```bash
python SU2_PY/aeroacoustics.py history.csv \
  --time Cur_Time \
  --pressure mic1 \
  --pressure-scale 101325
```

The default acoustic reference is 20 µPa, appropriate for airborne sound. Use
`--reference-pressure` for another medium or convention.

## Sampling guidance

- Resolve the highest frequency of interest with adequate samples per period;
  the absolute Nyquist limit is half the sample rate.
- Start collecting samples after initial transients have left the domain, or
  remove them with `--skip-samples`.
- Increase `--segment-length` for finer frequency resolution. Frequency-bin
  width is the sample rate divided by segment length.
- Longer signals give more Welch segments and a smoother PSD estimate.
- The default Hann window, 50% overlap, and mean removal are suitable general
  choices. `--detrend linear` can remove slow numerical drift.
- Avoid interpreting probe pressure as a far-field acoustic result when the
  probe is inside a vortical or hydrodynamic near field. Permeable-surface
  Ffowcs Williams–Hawkings propagation is a larger, separate capability.

The spectrum CSV contains one frequency column followed by PSD, PSD-level, and
narrowband-SPL columns for every requested pressure signal. The optional JSON
summary contains scalar metrics suitable for automated comparisons.
