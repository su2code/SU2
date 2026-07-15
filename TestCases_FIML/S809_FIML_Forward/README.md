# S809 Airfoil FIML Test Case

## Description
S809 airfoil with Field Inversion Machine Learning (FIML) turbulence model corrections.
Ported from SU2 v5.0 FIML implementation to v8.5.0 "Harrier".

## Test Case Details
- **Airfoil**: S809 wind turbine airfoil
- **Mesh**: 35,550 points (coarse structured grid)
- **Flow Conditions**:
  - Incompressible RANS
  - Reynolds number: Re = 2.0E6 (based on chord)
  - Angle of Attack: ~14 degrees
  - Freestream velocity: 47.8 m/s
  - Target CL: 1.0546 (from DNS/experiments)

## Configuration Files

### 1. config_fiml_forward.cfg
**Purpose**: Forward simulation with pre-trained neural network weights

**Settings**:
- `TRAIN_NN= NO` - Use existing weights (no training)
- `KIND_TRAIN_NN= WEIGHTS` - Weights as design variables mode
- `ITER_START_NN= 100` - Start applying beta corrections after 100 iterations

**Usage**:
```bash
SU2_CFD config_fiml_forward.cfg
```

### 2. config_fiml_training.cfg
**Purpose**: Train neural network to match beta target field

**Settings**:
- `TRAIN_NN= YES` - Enable training
- `KIND_TRAIN_NN= BACKPROP` - Use backpropagation
- `BETA_TARGET_FILE= beta_target.dat` - Load target beta values
- `NUM_EPOCH= 100` - Train for 100 epochs
- `LEARNING_RATE= 0.01` - Gradient descent learning rate
- `ITER_START_NN= 1000` - Start training after flow converges (1000 iters)

**Usage**:
```bash
SU2_CFD config_fiml_training.cfg
```

## Neural Network Architecture
- **Input features**: 4 (prod/dest ratio, chi, delta criterion, S/Omega ratio)
- **Hidden layers**: 3
- **Neurons per layer**: 20
- **Activation**: tanh
- **Scaling**: Box-Cox transformation + Z-score normalization
- **Output**: beta_fiml correction factor (applied to SA production term)

## Files

### Required Files
- `S809_struct_coarse.su2` - Mesh file (35,550 points)
- `S809_struct_coarse.geo` - Geometry file
- `config_fiml_forward.cfg` - Forward mode configuration
- `config_fiml_training.cfg` - Training mode configuration

### Optional Files
- `beta_target.dat` - Beta target values for training (35,550 lines)
  - Format: One beta value per line
  - Order: Same as mesh point ordering
  - Example: Values in range [0.5, 1.5], mean=1.0

## Workflow

### Option A: Forward Mode (Use Pre-trained Weights)
1. Initialize NN weights to zero (or load from file)
2. Run forward simulation:
   ```bash
   SU2_CFD config_fiml_forward.cfg
   ```
3. Monitor convergence and CL/CD values

### Option B: Training Mode (Learn from Beta Targets)
1. Generate or provide `beta_target.dat` file
2. Run training simulation:
   ```bash
   SU2_CFD config_fiml_training.cfg
   ```
3. Monitor training loss convergence
4. Extract trained weights for future use

### Option C: Complete Workflow (Not yet implemented)
1. Run external optimization to find optimal beta field
2. Save optimized beta as `beta_target.dat`
3. Train NN using config_fiml_training.cfg
4. Use trained NN for design studies with config_fiml_forward.cfg

## Expected Results

### Baseline SA Model (without FIML)
- CL ≈ 0.95 (under-predicts by ~10%)
- CD ≈ 0.025

### With FIML Correction
- CL → 1.0546 (matches DNS/experimental target)
- CD → improved prediction
- Better prediction of separation/stall behavior

## Notes
- The `beta_target.dat` file in this directory contains random example values for testing
- For actual applications, beta targets should come from:
  - DNS/LES data assimilation
  - Experimental data matching
  - External optimization (Python, IPOPT, etc.)
- Training converges in ~100 epochs typically
- Forward mode is fast (no training overhead)

## References
- Original FIML implementation: Parish & Duraisamy, "A paradigm for data-driven predictive modeling using field inversion and machine learning" (2016)
- SU2 FIML v5.0: https://github.com/su2code/SU2/tree/v5.0
