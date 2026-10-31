"""Structured quadrilateral mesh of a 2D channel in SU2 format."""
import sys


def write_channel(filename, length, height, nx, ny):
    def node(i, j):
        return j * (nx + 1) + i

    with open(filename, "w") as f:
        f.write("NDIME= 2\n")
        f.write(f"NELEM= {nx * ny}\n")
        e = 0
        for j in range(ny):
            for i in range(nx):
                f.write(f"9 {node(i, j)} {node(i + 1, j)} {node(i + 1, j + 1)} {node(i, j + 1)} {e}\n")
                e += 1
        f.write(f"NPOIN= {(nx + 1) * (ny + 1)}\n")
        for j in range(ny + 1):
            for i in range(nx + 1):
                f.write(f"{length * i / nx:.10e} {height * j / ny:.10e} {node(i, j)}\n")
        markers = {
            "inlet": [(node(0, j + 1), node(0, j)) for j in range(ny)],
            "outlet": [(node(nx, j), node(nx, j + 1)) for j in range(ny)],
            "bottom": [(node(i, 0), node(i + 1, 0)) for i in range(nx)],
            "top": [(node(i + 1, ny), node(i, ny)) for i in range(nx)],
        }
        f.write(f"NMARK= {len(markers)}\n")
        for tag, lines in markers.items():
            f.write(f"MARKER_TAG= {tag}\nMARKER_ELEMS= {len(lines)}\n")
            for a, b in lines:
                f.write(f"3 {a} {b}\n")


if __name__ == "__main__":
    write_channel(sys.argv[1] if len(sys.argv) > 1 else "channel.su2", 0.2, 0.01, 200, 4)
