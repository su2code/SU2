"""
Hybrid ML-SU2 Coupling Example
Demonstrates real-time coupling between SU2 solver and PyTorch ML model
"""
from mpi4py import MPI
import torch
import torch.nn as nn

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

class SimpleSurrogate(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(SimpleSurrogate, self).__init__()
        self.fc1 = nn.Linear(input_dim, 64)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(64, output_dim)

    def forward(self, x):
        return self.fc2(self.relu(self.fc1(x)))

if rank == 0:
    print("Initializing SU2-PyTorch Hybrid Coupling Example...")
    
    surrogate_model = SimpleSurrogate(input_dim=1, output_dim=1)
    optimizer = torch.optim.Adam(surrogate_model.parameters(), lr=0.001)
    criterion = nn.MSELoss()
    
    for i in range(10):
        dummy_input = torch.randn(1, 1)
        dummy_target = torch.randn(1, 1)
        optimizer.zero_grad()
        output = surrogate_model(dummy_input)
        loss = criterion(output, dummy_target)
        loss.backward()
        optimizer.step()
    
    print("Setup complete. Ready for hybrid simulation.")

MPI.Finalize()
