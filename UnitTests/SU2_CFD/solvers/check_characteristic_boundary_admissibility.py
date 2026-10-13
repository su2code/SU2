import subprocess
import sys

for tag, expected in [
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
