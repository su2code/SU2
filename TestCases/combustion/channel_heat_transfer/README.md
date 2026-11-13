# Laminar channel heat transfer with Cantera (compressible solver)

Half of a planar channel (full height 0.25 mm, length 20 mm, symmetry plane at y = 0) with an O2/N2
mixture from `h2o2.yaml`. The inlet and outlet are Riemann boundaries (total conditions, static pressure).

| Case | Config | Check |
|------|--------|-------|
| Isothermal wall, M = 0.05 | `channel_isothermal.cfg` | `nusselt.py --wall_T 300`, fully developed Nu = 7.541 |
| Uniform heat flux, M = 0.05 | `channel_heatflux.cfg` | `nusselt.py`, fully developed Nu = 8.235 |
| Isothermal wall, M = 0.01, `LOW_MACH_CORR` | `channel_lowmach.cfg` | `compare_solvers.py` with `channel_lowmach_inc.cfg` |
| Conjugate heat transfer with a solid wall | `cht_channel.cfg` | `energy_balance.py` |

The Nusselt numbers are those of Shah & London, *Laminar Flow Forced Convection in Ducts* (1978), for
parallel plates with a fully developed velocity and temperature profile. The post-processing scripts read
the Tecplot ASCII volume output and need the Python interface of Cantera.

Results on these meshes:

- Isothermal wall: Nu = 7.40 to 7.61 for 0.014 < x* < 0.06, with x* = x / (D_h Re Pr).
- Uniform heat flux: Nu = 8.26 to 8.31 for 0.03 < x* < 0.08.
- Low Mach: mass flow 0.65 % and enthalpy flux 0.02 to 1.4 % from the incompressible solver.
- Conjugate heat transfer: enthalpy drop of the fluid 232.03 W/m, wall heat flux 231.93 W/m, heat through
  the outer face of the solid 232.72 W/m.

Near the outlet the local Nu rises above the fully developed value, where the bulk to wall temperature
difference is small and viscous heating is no longer negligible.

Meshes: `python3 make_mesh.py` (fluid) and `python3 make_mesh.py --cht --output channel_cht.su2`.
