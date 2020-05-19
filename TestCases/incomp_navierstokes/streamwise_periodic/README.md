# Streamwise Periodicity testcases

All Testcases use the incompressible solver implemented by Thomas Economon.

## `pipe_slice_3D`

Overview: Hagen Poiseuille flow through a 1-primal-cell thick pipe slice in 3D.

Analytical solution of the velocity magnitude for steady laminar pipe flow in a round pipe `v_mag (r) = -1/(4*mu) * (Delta p / Delta x) * (R**2 - r**2)` therefore a pressure drop Delta p is prescribed, heated walls

`Re = rho * v * L / mu = 1.0 * ? * 5e-3 / 1.8e-5` bei v -> averaged (mass or area weighted?) velocity the critical Re ~= 2300 (v = 0.6 for now) makes Re=167

It would nice to have a Re ~= 1500 to have a better testcase (achieve that with v~5 or 6 i.e. scale Delta P by factor 10 from 0.001 to 0.01)

##  `half_cylinder_2D`  
half cylinder massflow prescribed heated cylinder - probably discontinued

## 2D_pinArray_dp_hf

## 2D_pinArray_mf

## 2D_pinArray_cht_dp_hf

### Discrete Adjoint

## 3D_pinArray_mf_hf

## 3D_pinArray_cht_dp_hf

### Discrete Adjoint