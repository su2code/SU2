#!/usr/bin/env python3
"""Structured C-mesh around an airfoil in Selig format, built with gmsh and written in SU2 format.

All boundary and grid-line nodes are placed explicitly; gmsh fills each block between two grid
lines by transfinite interpolation. Lengths are in chords unless stated otherwise.
"""
import argparse

import gmsh
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq


def read_selig(fn):
    pts = []
    for line in open(fn):
        try:
            pts.append([float(v) for v in line.split()[:2]])
        except ValueError:
            continue
    return np.array([p for p in pts if len(p) == 2])


class Airfoil:
    """Cubic spline through the coordinates, parametrized by arc length, split at the leading edge."""

    def __init__(self, xy):
        t = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(xy, axis=0).T))])
        self.spline = CubicSpline(t, xy)
        fine = np.linspace(0.0, t[-1], 200001)
        p = self.spline(fine)
        self.s = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
        self.t = fine
        self.s_le = self.s[np.argmin(p[:, 0])]
        self.length = self.s[-1]

    def point(self, s):
        return self.spline(np.interp(s, self.s, self.t))

    def tangent(self, s):
        d = self.spline(np.interp(s, self.s, self.t), 1)
        return d / np.linalg.norm(d, axis=-1, keepdims=True)


def cluster(length, ds0, ds1, ds_max, w0, w1):
    """Node positions on [0, length] with spacing ds0 at 0, ds1 at length and at most ds_max between."""
    x = np.linspace(0.0, length, 20001)
    f = ds_max - (ds_max - ds0) * np.exp(-x / w0) - (ds_max - ds1) * np.exp(-(length - x) / w1)
    g = np.concatenate([[0.0], np.cumsum(0.5 * (1 / f[1:] + 1 / f[:-1]) * np.diff(x))])
    n = int(np.ceil(g[-1]))
    return np.interp(np.linspace(0.0, g[-1], n + 1), g, x)


def geometric(length, ds0, ratio):
    """Node positions on [0, length] with first spacing ds0 and constant growth ratio."""
    n = int(np.ceil(np.log(1 + length * (ratio - 1) / ds0) / np.log(ratio)))
    r = brentq(lambda r: ds0 * (r**n - 1) / (r - 1) - length, 1 + 1e-9, 2 * ratio)
    return np.concatenate([[0.0], np.cumsum(ds0 * r ** np.arange(n))])


def stretch(length, ds0, n, r1, height):
    """n node positions on [0, length]: growth r1 up to height, then a smooth rise to the ratio that fits."""
    kb = np.log(1 + height * (r1 - 1) / ds0) / np.log(r1)
    k = np.arange(n - 1)

    def spacing(r2):
        r = r1 + (r2 - r1) * 0.5 * (1 + np.tanh((k - kb) / 6.0))
        return ds0 * np.concatenate([[1.0], np.cumprod(r[:-1])])

    r2 = brentq(lambda r2: spacing(r2).sum() - length, 0.5, 3.0)
    return np.concatenate([[0.0], np.cumsum(spacing(r2))])


def resample(curve, s_nodes):
    s = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(curve, axis=0).T))])
    return np.column_stack([np.interp(s_nodes * s[-1], s, curve[:, 0]), np.interp(s_nodes * s[-1], s, curve[:, 1])])


def c_line(foil, args, sc):
    """Wall arc lengths from the lower trailing edge, wake distances, and the leading-edge index."""
    lower = cluster(foil.s_le, args.ds_te / sc, args.ds_le / sc, args.ds_lower / sc, 0.04, 0.08)
    upper = cluster(foil.length - foil.s_le, args.ds_le / sc, args.ds_te / sc, args.ds_upper / sc, 0.04, 0.08)
    d = geometric(args.wake, args.ds_te / sc, args.wake_ratio ** (1 / sc))
    return np.concatenate([lower, foil.s_le + upper[1:]]), d, len(lower) - 1


def c_coord(foil, s_wall, d):
    """Position of the C-line nodes, from the lower outlet to the upper outlet."""
    return np.concatenate([-d[:0:-1], s_wall, foil.length + d[1:]])


def build(args, ref):
    """Grid lines of the C-mesh; far-field nodes are spread as for the spacings in ref at scale 1."""
    foil = Airfoil(read_selig(args.airfoil)[::-1])
    sc = args.scale

    # Wall nodes, ordered from the lower trailing edge over the leading edge to the upper trailing edge.
    s_wall, d, i_le = c_line(foil, args, sc)
    wall = foil.point(s_wall)
    tang = foil.tangent(s_wall)
    wall[0] = wall[-1] = 0.5 * (wall[0] + wall[-1])
    n_wall = np.column_stack([-tang[:, 1], tang[:, 0]])

    # Wake cut leaves along the trailing-edge bisector and bends back to the trailing-edge height.
    te = wall[-1]
    bis = tang[-1] - tang[0]
    bis /= np.linalg.norm(bis)
    slope = bis[1] / bis[0]
    yw = te[1] + slope * d * np.exp(-d / 0.5)
    wake = np.column_stack([te[0] + d, yw])
    dy = slope * np.exp(-d / 0.5) * (1 - d / 0.5)
    t_wake = np.column_stack([np.ones_like(dy), dy]) / np.hypot(1, dy)[:, None]
    n_wake_up = np.column_stack([-t_wake[:, 1], t_wake[:, 0]])
    n_te = np.array([-bis[1], bis[0]])

    n_wall[0], n_wall[-1] = -n_te, n_te

    # C-line from the lower outlet along the wake, around the airfoil and back along the wake.
    nw, nk = len(wall), len(wake)
    c_pts = np.vstack([wake[:0:-1], wall, wake[1:]])
    c_dir = np.vstack([-n_wake_up[:0:-1], n_wall, n_wake_up[1:]])
    c_extra = np.concatenate([d[:0:-1], np.zeros(nw), d[1:]]) * args.wake_y1_rate
    c_key = [("k", j) for j in range(nk - 1, 0, -1)] + [("w", i) for i in range(nw)] + [("k", j) for j in range(1, nk)]
    c_key[nk - 1] = ("w", nw - 1)
    i_le = nk - 1 + i_le

    # Spread the turn of the start directions at the trailing-edge corners, where grid lines converge.
    s_c = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(c_pts, axis=0).T))])
    ang = np.unwrap(np.arctan2(c_dir[:, 1], c_dir[:, 0]))
    for k in (nk - 1, nk - 1 + nw - 1):
        win = np.abs(s_c - s_c[k]) < 6 * args.te_smooth
        idx = np.flatnonzero(win)
        kern = np.exp(-0.5 * ((s_c[idx, None] - s_c[None, idx]) / args.te_smooth) ** 2)
        smooth = kern @ ang[idx] / kern.sum(axis=1)
        w = np.exp(-0.5 * ((s_c[idx] - s_c[k]) / (2 * args.te_smooth)) ** 2)
        ang[idx] = (1 - w) * ang[idx] + w * smooth
    c_dir = np.column_stack([np.cos(ang), np.sin(ang)])

    # Far field: straight lines from the outlet to an arc of radius R centred at the trailing edge.
    R, x_out = args.farfield, wake[-1, 0]
    th = np.linspace(1.5 * np.pi, 0.5 * np.pi, 4001)
    far = np.vstack([np.column_stack([np.linspace(x_out, te[0], 4001), np.full(4001, -R)])[:-1],
                     np.column_stack([te[0] + R * np.cos(th), R * np.sin(th)]),
                     np.column_stack([np.linspace(te[0], x_out, 4001), np.full(4001, R)])[1:]])
    sigma = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(far, axis=0).T))])

    # Each grid line ends where its wall-normal ray meets the far field, made monotone and blended with uniform.
    q = c_pts + c_dir * args.tension
    hit = np.empty(len(q))
    for i in range(len(q)):
        r = far - q[i]
        hit[i] = sigma[np.argmax(r @ c_dir[i] / np.linalg.norm(r, axis=1))]
    hit = 0.5 * (np.maximum.accumulate(hit) + np.minimum.accumulate(hit[::-1])[::-1])
    s_ref, d_ref, _ = c_line(foil, ref, 1.0)
    xi_ref = c_coord(foil, s_ref, d_ref)
    uniform = np.interp(c_coord(foil, s_wall, d), xi_ref, np.linspace(0.0, sigma[-1], len(xi_ref)))
    hit = (1 - args.far_uniform) * hit + args.far_uniform * uniform
    c_far = np.column_stack([np.interp(hit, sigma, far[:, 0]), np.interp(hit, sigma, far[:, 1])])

    y1 = args.y1 / sc
    nn = int(round((args.n_normal - 1) * sc)) + 1
    r1 = args.growth ** (1 / sc)
    kb = np.log(1 + args.bl_height * (r1 - 1) / y1) / np.log(r1)
    ramp = np.minimum(np.arange(nn), kb)

    def grid_line(i):
        """Wall-normal ray for the tension length, blended over the blend length into the line to the far point."""
        p0, d0, p1 = c_pts[i], c_dir[i], c_far[i]
        dist = np.linalg.norm(p1 - p0)
        t = np.linspace(0.0, 1.0, 20001)[:, None]
        w = np.clip((t * dist - args.tension) / args.blend, 0.0, 1.0)
        w = w * w * (3 - 2 * w)
        curve = (1 - w) * (p0 + d0 * dist * t) + w * (p0 + (p1 - p0) * t)
        arc = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(curve, axis=0).T))])[-1]
        h = stretch(arc, y1, nn, r1, args.bl_height)
        return resample(curve, (h + c_extra[i] * ramp * (1 - h / arc)) / arc)

    lines = [grid_line(i) for i in range(len(c_pts))]

    # The two halves of the C-line, both running from the leading edge to the outlet.
    half_up = [(c_key[i], lines[i]) for i in range(i_le, len(c_pts))]
    half_lo = [(c_key[i], lines[i]) for i in range(i_le, -1, -1)]
    i_te_up, i_te_lo = nk - 1 + nw - 1 - i_le, i_le - (nk - 1)

    return dict(nn=nn, halves=[(half_up, i_te_up, +1), (half_lo, i_te_lo, -1)])


def mesh(g, args):
    """Blocks between every args.block-th grid line, filled by transfinite interpolation into quads."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("cmesh")
    geo = gmsh.model.geo
    points, lines = {}, {}
    nn = g["nn"]

    def point(key, xy):
        if key not in points:
            points[key] = geo.addPoint(xy[0] * args.chord, xy[1] * args.chord, 0.0)
        return points[key]

    def line(a, b):
        if (a, b) in lines:
            return lines[(a, b)]
        if (b, a) in lines:
            return -lines[(b, a)]
        lines[(a, b)] = geo.addLine(a, b)
        return lines[(a, b)]

    airfoil, farfield = [], []
    for half, i_te, side in g["halves"]:
        def node(i, j):
            """Point tag of node j on grid line i, created on first use."""
            key, gl = half[i]
            if j > 0:
                key = ("le", j) if i == 0 else (key, side, j)
            return point(key, gl[j])

        top_j = nn - 1
        cuts = sorted(set(list(range(0, len(half), args.block)) + [i_te, len(half) - 1]))
        for a, b in zip(cuts[:-1], cuts[1:]):
            bottom = [line(node(i, 0), node(i + 1, 0)) for i in range(a, b)]
            right = [line(node(b, j), node(b, j + 1)) for j in range(nn - 1)]
            top = [line(node(i + 1, top_j), node(i, top_j)) for i in range(b - 1, a - 1, -1)]
            left = [line(node(a, j + 1), node(a, j)) for j in range(nn - 2, -1, -1)]
            loop = bottom + right + top + left
            if side < 0:
                loop = [-c for c in reversed(loop)]
            s = geo.addPlaneSurface([geo.addCurveLoop(loop)])
            geo.mesh.setTransfiniteSurface(s, cornerTags=[node(a, 0), node(b, 0), node(b, top_j), node(a, top_j)])
            geo.mesh.setRecombine(2, s)
            if b <= i_te:
                airfoil += [abs(c) for c in bottom]
            farfield += [abs(c) for c in top]
            if b == len(half) - 1:
                farfield += [abs(c) for c in right]
    geo.synchronize()
    for c in gmsh.model.getEntities(1):
        gmsh.model.mesh.setTransfiniteCurve(c[1], 2)
    gmsh.model.addPhysicalGroup(1, sorted(set(airfoil)), name="airfoil")
    gmsh.model.addPhysicalGroup(1, sorted(set(farfield)), name="farfield")
    gmsh.model.addPhysicalGroup(2, [s[1] for s in gmsh.model.getEntities(2)], name="fluid")
    gmsh.model.mesh.generate(2)


def quality():
    """Print node and cell counts, corner angles, inverted cells and the largest aspect ratio."""
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xy = np.zeros((int(tags.max()) + 1, 2))
    xy[tags.astype(int)] = xyz.reshape(-1, 3)[:, :2]
    _, _, conn = gmsh.model.mesh.getElements(2)
    q = xy[np.concatenate(conn).astype(int).reshape(-1, 4)]
    e = np.roll(q, -1, axis=1) - q
    ln = np.linalg.norm(e, axis=2)
    cosang = -np.einsum("ijk,ijk->ij", e, np.roll(e, 1, axis=1)) / (ln * np.roll(ln, 1, axis=1))
    ang = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
    cross = e[:, :, 0] * np.roll(e, 1, axis=1)[:, :, 1] - e[:, :, 1] * np.roll(e, 1, axis=1)[:, :, 0]
    print(f"nodes {len(tags)}, quads {len(q)}")
    print(f"corner angle min {ang.min():.1f} deg, max {ang.max():.1f} deg, "
          f"cells with an angle below 45 deg: {(ang.min(axis=1) < 45).sum()}")
    print(f"inverted cells: {(np.sign(cross) != np.sign(cross[:, :1])).any(axis=1).sum()}")
    print(f"aspect ratio max {(ln.max(axis=1) / ln.min(axis=1)).max():.0f}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--airfoil", default="s809.dat")
    p.add_argument("--output", default="s809_cmesh.msh")
    p.add_argument("--chord", type=float, default=0.6, help="chord in mesh units")
    p.add_argument("--scale", type=float, default=1.0, help="refinement factor applied to all spacings")
    p.add_argument("--y1", type=float, default=3e-6, help="first wall spacing")
    p.add_argument("--growth", type=float, default=1.1, help="wall-normal growth ratio in the boundary layer")
    p.add_argument("--bl_height", type=float, default=0.04, help="height of the region with constant growth")
    p.add_argument("--n_normal", type=int, default=129, help="nodes on each wall-normal grid line")
    p.add_argument("--ds_le", type=float, default=1e-4, help="wall spacing at the leading edge")
    p.add_argument("--ds_te", type=float, default=6e-4, help="wall spacing at the trailing edge")
    p.add_argument("--ds_upper", type=float, default=3.5e-3, help="largest spacing on the upper surface")
    p.add_argument("--ds_lower", type=float, default=8e-3, help="largest spacing on the lower surface")
    p.add_argument("--wake", type=float, default=100.0, help="wake length behind the trailing edge")
    p.add_argument("--wake_ratio", type=float, default=1.06, help="growth ratio along the wake")
    p.add_argument("--wake_y1_rate", type=float, default=2e-3, help="increase of the wake-cut spacing per chord downstream")
    p.add_argument("--farfield", type=float, default=100.0, help="far-field radius")
    p.add_argument("--tension", type=float, default=0.1, help="length over which grid lines leave the wall normally")
    p.add_argument("--te_smooth", type=float, default=0.06, help="smoothing length of start directions at the trailing edge")
    p.add_argument("--blend", type=float, default=1.5, help="length over which grid lines turn to the far point")
    p.add_argument("--far_uniform", type=float, default=0.5, help="blend of far-field node spacing towards uniform")
    p.add_argument("--block", type=int, default=8, help="C-line nodes per transfinite block")
    p.add_argument("--gui", action="store_true")
    args = p.parse_args()
    keys = ("ds_le", "ds_te", "ds_upper", "ds_lower", "wake_ratio")
    ref = argparse.Namespace(**{**vars(args), **{k: p.get_default(k) for k in keys}})

    g = build(args, ref)
    mesh(g, args)
    quality()
    gmsh.write(args.output)
    if args.gui:
        gmsh.fltk.run()
    gmsh.finalize()


if __name__ == "__main__":
    main()
