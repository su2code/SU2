"""Structured meshes of a laminar half channel (symmetry plane at y = 0, wall at y = h), and of a solid wall layer."""
import argparse


def stretched(n, ratio):
    """Coordinates in [0, 1] clustered toward 1 with a geometric cell growth ratio."""
    if abs(ratio - 1.0) < 1e-12:
        return [i / n for i in range(n + 1)]
    sizes = [ratio ** (n - 1 - i) for i in range(n)]
    total = sum(sizes)
    y = [0.0]
    for s in sizes:
        y.append(y[-1] + s / total)
    return y


def write_block(f, xs, ys, markers):
    nx, ny = len(xs) - 1, len(ys) - 1

    def node(i, j):
        return j * (nx + 1) + i

    edges = {
        "south": [(node(i, 0), node(i + 1, 0)) for i in range(nx)],
        "east": [(node(nx, j), node(nx, j + 1)) for j in range(ny)],
        "north": [(node(i + 1, ny), node(i, ny)) for i in range(nx)],
        "west": [(node(0, j + 1), node(0, j)) for j in range(ny)],
    }
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
            f.write(f"{xs[i]:.12e} {ys[j]:.12e} {node(i, j)}\n")
    f.write(f"NMARK= {len(markers)}\n")
    for side, tag in markers.items():
        f.write(f"MARKER_TAG= {tag}\nMARKER_ELEMS= {len(edges[side])}\n")
        for a, b in edges[side]:
            f.write(f"3 {a} {b}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--length", type=float, default=0.02)
    p.add_argument("--half_height", type=float, default=1.25e-4)
    p.add_argument("--solid_thickness", type=float, default=1.25e-4)
    p.add_argument("--nx", type=int, default=160)
    p.add_argument("--ny", type=int, default=16)
    p.add_argument("--ny_solid", type=int, default=8)
    p.add_argument("--ratio", type=float, default=1.1, help="wall-normal growth ratio toward the wall")
    p.add_argument("--output", default="channel_fluid.su2")
    p.add_argument("--cht", action="store_true", help="two zones, the fluid and a solid layer above the wall")
    a = p.parse_args()

    xs = [a.length * i / a.nx for i in range(a.nx + 1)]
    ys = [a.half_height * y for y in stretched(a.ny, a.ratio)]
    fluid = {"south": "symmetry", "east": "outlet", "north": "wall", "west": "inlet"}
    with open(a.output, "w") as f:
        if not a.cht:
            write_block(f, xs, ys, fluid)
            return
        ys_solid = [a.half_height + a.solid_thickness * y for y in stretched(a.ny_solid, 1.0)]
        solid = {"south": "solid_wall", "east": "solid_outlet", "north": "solid_outer", "west": "solid_inlet"}
        f.write("NZONE= 2\nIZONE= 1\n")
        write_block(f, xs, ys, fluid)
        f.write("IZONE= 2\n")
        write_block(f, xs, ys_solid, solid)


if __name__ == "__main__":
    main()
