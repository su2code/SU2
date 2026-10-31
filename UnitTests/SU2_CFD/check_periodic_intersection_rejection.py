"""Run the MPI-aborting geometry check in a separate process."""

import subprocess
import sys

result = subprocess.run(
    sys.argv[1:] + ["Invalid rotational intersection"],
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True,
    timeout=60,
)
print(result.stdout)
assert result.returncode != 0, "Rotational intersection was accepted"
assert (
    "Intersecting periodic pairs involving rotation are not supported yet."
    in result.stdout
)
