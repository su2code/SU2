# PR Details

## Branch Name
`fix/2d-actuator-disk-and-periodic-solver-bugs`

## Commit Message
```text
Fix bugs in 2D Actuator Disk and Periodic Solver

- CSU2ASCIIMeshReaderFVM.cpp: Fix uniqueness logic and array indexing in 2D Actuator Disk splitting.
- CSolver.cpp: Fix matrix corruption in Periodic Least Squares (replace duplicate r23_b with 0.0).

Signed-off-by: shbhmexe <shubhushukla586@gmail.com>
```

## PR Title
Fix logic bugs in 2D Actuator Disk preprocessing and Periodic Solver

## PR Message
## Proposed Changes
This PR fixes two core logical issues:
1. **2D Actuator Disk Splitting**: Corrects the perimeter detection logic (removing an incorrect `unique()` call) and fixes two array indexing errors (`iEdge` -> `iPoint`, `EdgeBegin` -> `ActDiskPoint_Front`) that could cause incorrect connectivity or segmentation faults.
2. **Periodic Least Squares**: Corrects the communication buffer construction by removing a redundant and misplaced term (`r23_b`), ensuring the Upper Triangular matrix structure is preserved (setting the ZY term to `0.0`).

## Related Work
Fixes logical bugs found during code review. No related PRs or issues.

## PR Checklist
- [x] I am submitting my contribution to the develop branch.
- [x] My contribution generates no new compiler warnings (try with --warnlevel=3 when using meson).
- [x] My contribution is commented and consistent with SU2 style (https://su2code.github.io/docs_v7/Style-Guide/).
- [x] I used the pre-commit hook to prevent dirty commits and used `pre-commit run --all` to format old commits.
- [ ] I have added a test case that demonstrates my contribution, if necessary.
- [ ] I have updated appropriate documentation (Tutorials, Docs Page, config_template.cpp), if necessary.
```
