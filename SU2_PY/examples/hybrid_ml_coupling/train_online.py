"""
Hybrid Machine Learning Coupling Example for SU2
-------------------------------------------------------------------------
This script demonstrates how to couple the SU2 CFD solver with a 
PyTorch training loop in real-time using the Python Wrapper (pysu2).

It performs the following:
1. Initializes a lightweight PyTorch Surrogate Model.
2. Launches the SU2 solver in Single-Zone mode.
3. Performs a time-stepping loop where:
   - SU2 advances the physics simulation.
   - Flow data (RMS Density) is extracted from memory.
   - The PyTorch model is trained online to predict this state.

Usage:
    python3 train_online.py

Requirements:
    - SU2 compiled with -Denable-pywrapper=true -Dwith-mpi=enabled
    - mpi4py
    - torch
"""

import pysu2
from mpi4py import MPI
import torch
import torch.nn as nn
import sys

# --- 1. Define Surrogate Model ---
class SurrogateModel(nn.Module):
    def __init__(self):
        super(SurrogateModel, self).__init__()
        self.fc = nn.Linear(1, 1) 
    
    def forward(self, x):
        return self.fc(x)

def main():
    # Initialize MPI Communicator
    # SU2 requires a valid communicator even for single-zone runs
    comm = MPI.COMM_WORLD
    
    # Configuration
    config_file = "inv_NACA0012.cfg"
    
    # Initialize Driver
    try:
        driver = pysu2.CSinglezoneDriver(config_file, 1, comm)
    except TypeError:
        print("[Error] Failed to initialize driver. Ensure SU2 is compiled with MPI support.")
        sys.exit(1)

    driver.Preprocess(0)

    # Initialize ML Model
    model = SurrogateModel()
    optimizer = torch.optim.SGD(model.parameters(), lr=0.01)
    print("\n[Hybrid-ML] Models initialized. Starting coupled loop...\n")

    # --- Hybrid Simulation Loop ---
    n_iterations = 10
    
    for i in range(n_iterations):
        # A. Run Physics Step
        driver.Run()
        
        # B. Extract Physics Data (The Bridge)
        try:
            # We use the History Output value for this demonstration
            # Note: Ensure "RMS_DENSITY" is present in your config's HISTORY_OUTPUT
            physics_state = driver.GetOutputValue("RMS_DENSITY")
        except Exception:
            print(f"  Iter {i}: [Warning] Could not extract RMS_DENSITY. Using dummy value.")
            physics_state = 0.0
            
        print(f"  Iter {i}: Physics State (Log Rho) = {physics_state:.6f}")
        
        # C. Online Training Step
        # Train the model to predict the current physics state from the iteration index
        input_tensor = torch.tensor([[float(i)]])
        target_tensor = torch.tensor([[physics_state]])
        
        optimizer.zero_grad()
        prediction = model(input_tensor)
        loss = (prediction - target_tensor) ** 2
        loss.backward()
        optimizer.step()
        
        print(f"  Iter {i}: ML State (MSE Loss)     = {loss.item():.6f}")

    # Finalize
    driver.Postprocess()
    print("\n[Hybrid-ML] Simulation Complete.")

if __name__ == "__main__":
    main()
