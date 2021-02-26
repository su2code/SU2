# Gradient validation from start to finish

This guide steps you through the steps necessary to perform a validation of the discrete adjoint sensitivities using finite differences.

All necessary config files are present and this guide steps through the different tasks to do.

If you are lucky enough too have some cores to spare, 14 is a suitable substitution for the `<#cores>` placeholder.

## FFD-box creation
This step is optional as the provided mesh already contains FFD box. This is for completeness if a new mesh e.g. with different resolution is created.
In `configMaster.cfg` the mentioned options have to be uncommented and others commented if they appear twice in the config.
Note that (only!) for the FFD-box creation a `MARKER_HEATFLUX= ( fluid_symmetry ) is artificially is set to avoid an error. This has to be done to make the config-Postprocessing aware that this marker exists as it is used in `DV_MARKER`.
Call `SU2_DEF configMaster.cfg` which creates the new mesh with the name given in 'MESH_OUT_FILENAME'.

## Primal run
Run `mpirun -n <#cores> SU2_CFD configMaster.cfg`

## Discrete-Adjoint run
Rename\copy\symlink `restart_*.dat` -> `solution_*.dat`
Run `mpirun -n <#cores> SU2_CFD_AD DA_configMaster.cfg` and afterwards `SU2_DOT_AD DA_configMaster.cfg`

## Finite-Differences run
The `OUTER_ITER` is set low in order to be suitable for the regression test. Set that back the number given in the config.
For the full gradient validation uncomment all design variables of the `DEFINITION_DV` config option.
Run `finite_differences.py -f FD_configMaster.cfg -z 2 -n <#cores>`.

## Comparing results
Just plot the `of_grad.csv` and `FINDIFF/of_grad_findiff.csv` with your tool of choice. Paraview's `Line Chart View` is one option.
