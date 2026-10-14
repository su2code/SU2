import subprocess
import sys

# SU2_MPI::Error exits the test process, so run each hidden Catch case in isolation.
# Check the boundary diagnostic as well as the exit code: a later flux NaN is a different failure.
for tag, expected in [
    (
        "Turbo-invalid-pressure",
        "Unable to construct an admissible static-pressure outflow state",
    ),
    ("Turbo-backflow", "requires an inflow state for backflow"),
    (
        "Ordinary-invalid-pressure",
        "Unable to construct an admissible static-pressure outflow state",
    ),
    ("Ordinary-backflow", "requires an inflow state for backflow"),
    (
        "Boundary-invalid-pressure",
        "Unable to construct an admissible static-pressure outflow state",
    ),
    (
        "Boundary-backflow",
        "requires an inflow state for backflow",
    ),
]:
    result = subprocess.run(
        [sys.argv[1], f"[{tag}]", "-r", "compact"],
        capture_output=True,
        text=True,
        timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert expected in output, output
    assert "NaN detected" not in output, output
print("Invalid characteristic boundary states are rejected before flux assembly.")
