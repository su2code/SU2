import subprocess
import sys

for tag, expected in [
    ("Boundary-negative-density", "Nonphysical characteristic boundary density"),
    (
        "Boundary-invalid-thermo",
        "Nonphysical characteristic boundary thermodynamic state",
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
