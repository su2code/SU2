"""Hybrid ML-SU2 Coupling Example"""
from mpi4py import MPI
import torch
import torch.nn as nn
try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    print("PyTorch not available - ML features disabled")

try:
    import mpi4py.MPI as MPI
    MPI_AVAILABLE = True
except ImportError:
    MPI_AVAILABLE = False
    print("mpi4py not available - MPI features disabled")

# Later in main execution:
if __name__ == "__main__":
    if not TORCH_AVAILABLE or not MPI_AVAILABLE:
        print("Skipping example - missing dependencies")
        sys.exit(0)
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
    model = SimpleSurrogate(1, 1)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.MSELoss()

    for _ in range(5):
        x = torch.randn(1, 1)
        y = torch.randn(1, 1)
        optimizer.zero_grad()
        loss = criterion(model(x), y)
        loss.backward()
        optimizer.step()

MPI.Finalize()
