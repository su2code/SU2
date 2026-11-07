# SST turbulent kinetic energy in the compressible solver

The compressible SST solver stores total energy per unit mass as
`E = e + |u|²/2 + k`. Thermodynamic pressure and temperature use `e`;
total enthalpy and its convective flux include `k`. Flow Jacobians hold `k`
fixed because it is advanced by the separate turbulence solver.

This energy convention does not imply that every term of the full
Favre-averaged energy equation is modeled. The current SSTm formulation:

- omits the isotropic `-2 rho k I/3` term from turbulent stress;
- uses approximate strain-based production `P = mu_t 2S:S` before applying
  the production limiter (other production options modify this expression);
- omits molecular diffusion and turbulent transport of `k` from the total
  energy flux. Temperature heat conduction is still present.

These are modeling approximations, distinct from interpreting `k` as
thermal energy inside a numerical flux or its Jacobian. The
[Turbulence Modeling Resource SST definitions](https://tmbwg.github.io/turbmodels/sst.html)
explicitly describe the SSTm production approximation, and its
[compressible RANS implementation notes](https://tmbwg.github.io/turbmodels/implementrans.html)
allow omission of the energy transport term involving the gradient of `k`.

For unmodified strain production, its difference from the contraction of
SU2's deviatoric turbulent stress with the velocity gradient is
`2 mu_t (div u)²/3`, before limiting. At fixed total energy, additional
production of `k` reduces thermal energy. This difference vanishes for
incompressible flow and can matter in shocks and strong expansions. The
production limiter, sustaining terms, transition coupling, and production
modifiers prevent interpreting this expression as the actual production
excess in every configuration. Changing these approximations requires
choosing and validating a different SST formulation.
