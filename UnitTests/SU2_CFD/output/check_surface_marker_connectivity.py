#!/usr/bin/env python3
"""Check that surface sorting preserves cells across distinct marker-local indices.

Run against SU2_CFD, optionally with --launcher mpirun -np 4. Each probe runs in
an isolated temporary directory. Select --case box9 for larger MPI rank counts.
"""
import argparse
from pathlib import Path
import subprocess
import tempfile


CUBE_MARKERS = ("x_minus", "x_plus", "y_minus", "y_plus", "z_minus", "z_plus")


def canonical(nodes):
    return tuple(sorted(tuple(round(v, 12) for v in p) for p in nodes))


def box(n, markers):
    values = [i / (n - 1) for i in range(n)]
    faces = {}
    for axis in range(3):
        other = [i for i in range(3) if i != axis]
        for side in range(2):
            name = CUBE_MARKERS[2 * axis + side]
            faces[name] = []
            for i in range(n - 1):
                for j in range(n - 1):
                    nodes = []
                    for a, b in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)):
                        p = [0.0] * 3
                        p[axis] = side
                        p[other[0]] = values[a]
                        p[other[1]] = values[b]
                        nodes.append(p)
                    faces[name].append(canonical(nodes))
    return set(face for name in markers for face in faces[name])


def mesh_fixture(kind):
    if kind == "square":
        points = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
        dim = 2
        volume = "9 0 1 2 3 0"
        boundaries = [(3, 0), (1, 2), (0, 1), (2, 3)]
        vtk = 3
    else:
        points = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
        dim = 3
        volume = "10 0 1 2 3 0"
        boundaries = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)]
        vtk = 5
    markers = tuple("face" + str(i) for i in range(4))
    lines = ["NDIME= " + str(dim), "NELEM= 1", volume, "NPOIN= 4"]
    lines += [" ".join(map(str, p[:dim])) + " " + str(i) for i, p in enumerate(points)]
    lines += ["NMARK= 4"]
    for name, nodes in zip(markers, boundaries):
        lines += [
            "MARKER_TAG= " + name,
            "MARKER_ELEMS= 1",
            str(vtk) + " " + " ".join(map(str, nodes)),
        ]
    return (
        "\n".join(lines) + "\n",
        markers,
        {canonical([points[i] for i in face]) for face in boundaries},
    )


def read_surface(path):
    tokens = path.read_text().split()
    start = tokens.index("POINTS")
    n = int(tokens[start + 1])
    offset = start + 3
    points = [
        tuple(map(float, tokens[offset + 3 * i : offset + 3 * i + 3])) for i in range(n)
    ]
    start = tokens.index("CELLS")
    count = int(tokens[start + 1])
    offset = start + 3
    faces = []
    for _ in range(count):
        width = int(tokens[offset])
        offset += 1
        ids = list(map(int, tokens[offset : offset + width]))
        offset += width
        faces.append(canonical([points[i] for i in ids]))
    return faces


def run(binary, launcher, name):
    markers = CUBE_MARKERS
    with tempfile.TemporaryDirectory(prefix="su2-surface-marker-") as temporary:
        work = Path(temporary)
        if name in ("square", "tetrahedron"):
            mesh, markers, expected = mesh_fixture(name)
            (work / "mesh.su2").write_text(mesh)
            mesh_config = "MESH_FORMAT= SU2\nMESH_FILENAME= mesh.su2\n"
        else:
            n = {
                "single_hex": 2,
                "box3": 3,
                "box9": 9,
                "subset": 2,
                "single_marker": 2,
            }[name]
            if name == "subset":
                markers = ("x_minus", "z_plus")
            elif name == "single_marker":
                markers = ("x_minus",)
            expected = box(n, markers)
            mesh_config = (
                "MESH_FORMAT= BOX\nMESH_FILENAME= unused.su2\nMESH_BOX_SIZE= (%s,%s,%s)\nMESH_BOX_LENGTH= (1,1,1)\n"
                % (n, n, n)
            )
        far = ",".join(
            CUBE_MARKERS if name not in ("square", "tetrahedron") else markers
        )
        cfg = (
            """SOLVER= EULER
MACH_NUMBER= 0.5
ITER= 0
MGLEVEL= 0
CONV_NUM_METHOD_FLOW= ROE
MUSCL_FLOW= NO
TIME_DISCRE_FLOW= RUNGE-KUTTA_EXPLICIT
OUTPUT_FILES= (SURFACE_PARAVIEW_ASCII)
SURFACE_FILENAME= surface
HISTORY_OUTPUT= (ITER)
SCREEN_OUTPUT= (INNER_ITER)
MARKER_FAR= ("""
            + far
            + """)
MARKER_PLOTTING= ("""
            + ",".join(markers)
            + ")\n"
            + mesh_config
        )
        (work / "case.cfg").write_text(cfg)
        result = subprocess.run(
            launcher + [str(binary), "case.cfg"],
            cwd=work,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=120,
        )
        if result.returncode:
            raise RuntimeError(name + " solver failed:\n" + result.stdout[-4000:])
        faces = read_surface(work / "surface.vtk")
        if (
            len(faces) != len(expected)
            or len(set(faces)) != len(faces)
            or set(faces) != expected
        ):
            raise AssertionError(
                "%s: wrote %s cells, expected %s; missing %s, unexpected %s"
                % (
                    name,
                    len(faces),
                    len(expected),
                    len(expected - set(faces)),
                    len(set(faces) - expected),
                )
            )
        print(
            "PASS: "
            + name
            + " preserves all "
            + str(len(expected))
            + " distinct surface cells",
            flush=True,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument(
        "--case",
        action="append",
        choices=(
            "single_hex",
            "box3",
            "box9",
            "subset",
            "single_marker",
            "square",
            "tetrahedron",
        ),
    )
    parser.add_argument("--launcher", nargs=argparse.REMAINDER, default=[])
    args = parser.parse_args()
    for name in args.case or (
        "single_hex",
        "box3",
        "subset",
        "single_marker",
        "square",
        "tetrahedron",
    ):
        run(args.binary.resolve(), args.launcher, name)


if __name__ == "__main__":
    main()
