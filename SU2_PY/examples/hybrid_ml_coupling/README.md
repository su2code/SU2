# Hybrid ML-SU2 Coupling Example

This example demonstrates how to couple SU2 with PyTorch for Physics-Informed Machine Learning (PIML).

## Features
- Real-time data extraction from SU2 using `GetOutputValue()`
- Online training of ML surrogate model
- Integration with `CSinglezoneDriver` and `mpi4py`

## Requirements
- SU2 with Python wrapper
- PyTorch
- mpi4py

## Usage
```bash
python hybrid_ml_example.py
```

## Description
Extracts flow variables (e.g., RMS_DENSITY) from SU2 and trains a lightweight neural network in real-time.
