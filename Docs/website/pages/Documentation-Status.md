---
title: Documentation Improvements and Status
permalink: /docs_v7/Documentation-Status/
---

This page tracks the documentation overhaul: what has been implemented, what still needs review, and the next work to undertake. Update it when a batch of work is completed, with links to its changes and checks.

**Status date: 22 October 2026.** This work is implemented on the `docs_site` branches of [rois1995/SU2](https://github.com/rois1995/SU2/tree/docs_site) and [rois1995/su2code.github.io](https://github.com/rois1995/su2code.github.io/tree/docs_site). It is a development preview; adoption by the official repositories remains a separate review.

For the current build, see [Docs Check](../Docs-Check/), the [source and tool manifest]({{ '/build-info.json' | relative_url }}) and the [automated check results]({{ '/docs-check.json' | relative_url }}). Those reports complement this maintained account of the work; a successful build does not establish the scientific accuracy of every explanation or recommendation.

- [What has improved](#what-has-improved)
- [Corrections from the content and code review](#corrections-from-the-content-and-code-review)
- [Evidence and current limits](#evidence-and-current-limits)
- [Guide and model review](#guide-and-model-review)
- [What is still missing](#what-is-still-missing)
- [How to go further](#how-to-go-further)
- [Change history](#change-history)
- [Keeping this page current](#keeping-this-page-current)

## What has improved

**Implemented** means that the feature or text is present in the fork preview. **Draft / review needed** means that content exists but the stated review is incomplete. **Planned** means that the work has not been completed. These labels do not imply numerical validation.

| Area | Status | Improvement and where to find it |
| --- | --- | --- |
| Documentation structure | Implemented | [Using the Documentation](../Using-the-Documentation/) explains where to find setup choices, worked examples, model definitions, numerical evidence and implementation details. |
| Configuration reference | Implemented | The [Option Reference](../Option-Reference/) is generated from option registrations, enum maps and source comments, with template and regression examples. Unknown defaults and additional applicability restrictions are stated explicitly. |
| Output reference | Implemented | The [Output Fields Reference](../Output-Fields-Reference/) collects registered screen, history and volume fields. Runtime-generated names and solver-specific availability still require attention to the stated limits. |
| Errors and migration | Implemented | [Error Messages](../Error-Messages/), [Renamed and Removed Options](../Renamed-Options/) and [Changes in v8](../Guide-to-v8/) help interpret diagnostics and update old configurations. |
| Developer explanations | Implemented; ongoing maintenance | [Code Structure](../Code-Structure/), [regression testing](../Running-Regression-Tests/), [Automatic Differentiation](../Automatic-Differentiation/), [Debugging](../Debugging/), [Glossary](../Glossary/) and eight **How to Add...** guides connect explanations to the implementation. |
| Code reference | Implemented | [Doxygen](../Code-Reference/) is built with the site. [Docs Check](../Docs-Check/) reports approximate source-documentation coverage and Doxygen warnings. Publishing it does not mean every class is adequately documented. |
| Physical model explanations | Source/theory follow-up completed; expert validation pending | Dedicated pages cover [turbulence](../Turbulence-Models/), [transition](../Transition-Models/), [hybrid RANS/LES](../Hybrid-RANS-LES/) and [wall treatment](../Wall-Treatment/), including variants, implementation differences and known limits. Scientific completeness still needs expert review. |
| Case setup guidance | Nine drafts; source review completed | Guides cover external RANS, incompressible flow, unsteady/DES, transition, turbomachinery, CHT, species/combustion, shape optimization and NICFD. Shared introductions explain evidence tags and review status. |
| Worked cases and numerical comparisons | Existing collections retained | [Tutorials]({{ '/tutorials/home/' | relative_url }}) remain the home for complete worked examples; [V&V]({{ '/vandv/home/' | relative_url }}) remains the home for numerical comparisons. Setup guides link to those collections rather than creating a second case collection. |
| Search and navigation | Implemented | [Search]({{ '/search/' | relative_url }}) indexes the current documentation, Tutorials and V&V without the former per-page text truncation. Old v6 pages remain available for existing links but are excluded from this search. |
| Website maintenance | Implemented | Feed and sitemap generation were restored by consolidating the plugin configuration. Added/modified menu sections, entries and page titles remain red in the fork preview. |
| Automated checks | Implemented | Generator regression tests, local files/anchors, documented option names, enum keywords in configuration blocks, descriptions for new options and search coverage are checked during builds. Existing findings are recorded separately from new failures. |
| Local and online builds | Implemented | One build command stages both sources, generates references, builds Doxygen/Jekyll and checks the assembled site. Ruby and gem versions are fixed; each artifact records both source revisions, build flags and actual tool versions. |
| Review and publishing | Implemented on the forks | Both repositories have push and pull-request checks. Pull requests produce previews without deploying. Successful SU2 branch builds publish to the fork's Pages site; website-only changes are checked separately and can be published by selecting their commit SHA in the SU2 workflow. |

## Corrections from the content and code review

The 22 October review corrected 16 documentation or generation findings. The main changes were:

- **Allowed keywords:** commented-out enum entries no longer appear as accepted settings. Migration examples now show accepted SA option combinations and explain the exceptional combinations.
- **Defaults:** array storage variables are no longer presented as default values. Template examples are distinguished from parser defaults; values initialized elsewhere are left explicitly unextracted.
- **Physical explanations:** turbulent heat transport is described separately for compressible and incompressible implementations; SA inlet values are distinguished from eddy-viscosity ratios; the SST stress explanation identifies the implemented variant and exceptions.
- **Support and diagnostics:** transition pages state the lack of LM discrete-adjoint support and the unreliable BCM intermittency output. These are documented limitations, not solver fixes delivered by this documentation work.
- **Examples and troubleshooting:** flamelet option/default names, CGNS format guidance, migration anchors and misleading transition-tutorial explanations were corrected.
- **Developer guidance:** AD comparison advice now includes convergence and matching settings; regression coverage and the existing warnings-as-errors workflow wiring are described accurately.
- **Search:** removal of the text limit restored searchable output names near the ends of long reference pages.

The changes are recorded in SU2 commit [7862504f9d](https://github.com/rois1995/SU2/commit/7862504f9d77ee708e002f125b52b3c1eaf6f641) and website commit [be90edc](https://github.com/rois1995/su2code.github.io/commit/be90edcd0c520f57c7381f19a2cc18433916c22c).

## Evidence and current limits

The latest completed build-system milestone before this status page used SU2 [a14b4df35d](https://github.com/rois1995/SU2/commit/a14b4df35d3863fc79b48c781fc83b8cb222d912) and website [602189b](https://github.com/rois1995/su2code.github.io/commit/602189b07ccee63ed1a99dae1ecc3b792d8ed782). Its dated results were:

- Nine regression checks passed; 25,607 local references and 178 search entries were checked, with zero new findings.
- Both the [SU2 build/deployment](https://github.com/rois1995/SU2/actions/runs/39788838345) and [website checks](https://github.com/rois1995/su2code.github.io/actions/runs/39788836318) passed. The published artifact was independently checked against the local source revisions, dependency lock and highlighting settings.
- Ruby 3.3.12, Bundler 2.4.22 and Jekyll 3.9.5 matched locally and in CI. Doxygen and Graphviz versions are recorded rather than pinned, so identical source revisions do not promise byte-identical output across environments.

These numbers describe that milestone, not a permanent total. Use the current generated reports linked above for this build's source revisions and counts.

The checks cover internal references, source extraction and selected configuration facts. They do not check external URL availability, prove convergence, validate physics or guarantee that every solver-specific option combination is supported. Source extraction uses text parsing and can miss facts established later at runtime. Draft recommendations and model explanations still need review against their stated assumptions, implementation and references.

## Guide and model review

**Review record, 22 October 2026:** Codex performed a documentation/source comparison of all nine setup guides and four physical-model pages against [SU2 `41ea6a219a`](https://github.com/rois1995/SU2/tree/41ea6a219ac038574f3d4a5b130c0c1d8af66457). Option spellings and enum values were checked automatically; support, units, equations and selected examples were followed into the implementation. This is an assisted source and reference review, not an independent CFD expert sign-off. No new CFD cases, gradient comparisons or mesh/time-refinement campaigns were run.

| Page | Reviewed and corrected where needed | Evidence and remaining work |
| --- | --- | --- |
| [External RANS](../Setup-Guide-RANS-External-Aero/) | SST ambient terms, inlet defaults, model sensitivity, log10 linear residuals, wall and far-field advice. | `CConfig`, `ParseSSTOptions`, `turb_sources.hpp`, `CFlowCompOutput`; compare the selected SA/SST variant with TMR and test mesh/domain sensitivity. |
| [Incompressible](../Setup-Guide-Incompressible/) | Density/pressure coupling, thermal feedback through viscosity, enthalpy versus temperature outputs, scaling, inlet direction, 2D mass flow and FDS reconstruction. | `CIncEulerSolver`, `CConfig`; independently run the chosen density-based or pressure-based case and check conservation. |
| [Unsteady/DES](../Setup-Guide-Unsteady-DES/) | Original versus delayed DES length scales, grid measures, shielding and statistical convergence. | `CTurbSASolver::SetDES_LengthScale` and the DDES paper; a 2D regression does not validate 3D DES. Time-step, span and sampling studies remain. |
| [Transition setup](../Setup-Guide-Transition/) | Base-model coupling, turbulence decay, supported correlations and interpretation of wall shear. | `CTransLMSolver`, transition source terms and TMR; calibration ranges, roughness scaling and moving-frame use still need specialist review. |
| [Turbomachinery](../Setup-Guide-Turbomachinery/) | Axis conventions for axial/radial and 2D examples, boundary kinds, averaging and ramps. | Turbomachinery geometry, `CEulerSolver` boundaries and linked stage configurations; validate radial/3D performance and interface conservation. |
| [Heat transfer/CHT](../Setup-Guide-Heat-Transfer-CHT/) | Heat-flux sign, steady versus transient material roles, coupling and heat balance. | `CHeatSolver`, `CIncNSSolver` and CHT transfer routines; quantify non-matching-interface heat conservation and transient convergence. |
| [Species/combustion](../Setup-Guide-Species-Combustion/) | Mixture properties, Lewis-model limits, raw wall-flux scaling, clipping and source-inclusive species balances. | `CSpeciesSolver::BC_Wall_Generic`, fluid/diffusion models, `CSpeciesFlameletSolver`; the modified fresh-start flame example and manifold accuracy still require runs. |
| [Shape optimization](../Setup-Guide-Shape-Optimization/) | Adjoint cost/derivative scope, legacy versus new mesh deformation, physical boundary preservation and disabled validation flag. | `CDeformationDriver`, `CLinearElasticity`, `CMeshSolver`, `discrete_adjoint.py`; perform converged finite-difference comparisons and mesh-validity checks. |
| [NICFD](../Setup-Guide-NICFD/) | Phase stability, critical/acentric definitions, EOS accuracy limits and modified-example status. | SU2 fluid models and CoolProp primary documentation; validate the actual fluid/state range, transport closure and surrogate derivatives. |
| [Turbulence models](../Turbulence-Models/) | SA/SST variants and restrictions, Kato–Launder production, ambient terms and inlet definitions. | `option_structure.hpp`, turbulence source/solver files and TMR; implementation correspondence is not a complete numerical verification. |
| [Transition models](../Transition-Models/) | LM/SA coupling, BCM coefficient convention and documented output/adjoint limitations. | Transition source/solver files and TMR; solver limitations remain open, and specialist assessment is pending. |
| [Hybrid RANS/LES](../Hybrid-RANS-LES/) | Shielded length scale, grid measures, backscatter caveats and averaging limits. | SA length-scale and backscatter implementations plus original papers; full 3D validation and the noted backscatter unit discrepancy remain. |
| [Wall treatment](../Wall-Treatment/) | Wall-function restrictions, roughness/heat-transfer scope and conditional y+ guidance. | Flow/turbulence wall routines, TMR SA roughness and cited papers; low-y+ fallback and mixed-wall limitations remain solver work. |

Primary equation references checked include the TMR pages for [SA](https://tmbwg.github.io/turbmodels/spalart.html), [SST](https://tmbwg.github.io/turbmodels/sst.html), [LM-SST](https://tmbwg.github.io/turbmodels/langtrymenter_4eqn.html), [LM-SA](https://tmbwg.github.io/turbmodels/sa_gamma_retheta_3eqn.html) and [SA-BCM](https://tmbwg.github.io/turbmodels/sa-bc_1eqn.html), the [Spalart et al. DDES paper](https://doi.org/10.1007/s00162-006-0015-0), and [CoolProp EOS documentation](https://coolprop.org/fluid_properties/PurePseudoPure.html). All 18 TMR links/anchors on the reviewed pages resolved at this review. These checks do not endorse every empirical range in the guides; **Draft** and **to be reviewed** labels remain.

The Python configuration helper still writes the unsupported buffet-monitoring key in `SU2_PY/SU2/io/config.py`; that caller needs a separate fix and regression check before using its buffet path. Removing the template assignment does not repair that helper.

The same follow-up repaired all 180 recorded link findings: image/media paths, two image enlargement links, old page names and anchors, version paths, two tutorial file links, and 59 archived downloads. The latter now point to their official archive; all 59 downloads and both replacement case files returned HTTP 200. It also fixed four template spellings (`DNETTHRUST_DBCTHRUST`, `FLAME_INIT_IGNITION`, `HOLD_GRID_FIXED_COORD`, `DEFORM_POISSONS_RATIO`), removed four unsupported assignments, explained solid diffusivity through material properties, and changed the structural example to `PARAVIEW_LEGACY` to retain its stated `.vtk` output. The enum checker now handles hyphenated keywords, with a reproducing test.

## What is still missing

The current [accepted-findings list](https://github.com/rois1995/SU2/blob/docs_site/Docs/website/tools/check-baseline.json) retains 24 missing source descriptions. The link, template-name and enum baselines are empty after the follow-up above. Historical build counts earlier on this page remain dated evidence; use the generated reports for current counts.

| Priority | Remaining work | What completion requires |
| --- | --- | --- |
| Next | Improve option descriptions: 24 registrations lacked source descriptions; approximately 296 descriptions were very short. Enum choices also need explanations. | Describe effect, units, applicability, restrictions, interactions and conditional defaults in the source, prioritizing options used by the setup guides. Check the generated reference afterwards. |
| Next | Obtain independent expert review and numerical evidence for the source-reviewed guides and model pages. | Use the review table above to select representative cases, check convergence/conservation and compare with appropriate theory/data. Retain draft labels until the stated review is complete. |
| Then | Refresh older tutorials, version claims and contribution instructions as they are checked. | Confirm linked case files, accepted settings, run instructions and stated results for the declared version. Keep worked examples in Tutorials and numerical comparisons in V&V. |
| Then | Fill missing model explanations and keep known support limitations current. | Link exact variants and equation references to the SU2 implementation, define symbols/units, and distinguish verified behavior from proposed or unsupported capabilities. Extend the existing model pages where appropriate. |
| Then | Strengthen extraction and build checks when a new mismatch is demonstrated. | Add a focused reproducing test and fix its cause. Solver-specific configuration compatibility still needs appropriate SU2 checks; a general parser or framework rewrite is not an agreed requirement. |
| Then | Decide whether the remaining build tools need tighter version control. | Demonstrate that an unpinned tool causes a relevant difference, then document or pin it. Continue recording actual versions and both source SHAs. |
| Separate work | Verify and fix underlying solver defects noted during documentation review. | Reproduce the behavior, add an appropriate failing check, fix the implementation and update the affected documentation. A documentation caveat alone is not a solver repair. |
| Separate review | Propose the improvements to the official repositories. | Submit focused, reviewable changes with evidence and coordinated source/website revisions. Fork deployment does not establish upstream approval. |

## How to go further

1. **Keep the new checks passing.** Keep the resolved link/template/enum baselines empty. Correct new findings at their source; do not add them to an accepted-debt list.
2. **Improve the source descriptions used by the guides.** Work through one option family at a time. This improves the reference and explanations together without maintaining a competing option catalog.
3. **Review the existing guidance and evidence.** Use the linked Tutorials and V&V cases. Distinguish source checking, a successful run, numerical verification and physical validation; record which has actually been performed.
4. **Resolve implementation issues separately.** Update known-limit sections only when the relevant code change and its checks support doing so.
5. **Prepare upstream changes in coherent batches.** Keep tooling, corrected explanations and reviewed guidance understandable on their own, and identify the matching changes in the other repository.

These are proposed next priorities. The page records completed work; it does not imply that every backlog item is already assigned or underway.

## Change history

| Date | Delivered | Evidence |
| --- | --- | --- |
| 21 October 2026 | Generated references and Doxygen integration; developer and extension guides; model pages; nine draft setup guides; search, navigation and feed/sitemap repairs. | The linked pages above and the history of the two `docs_site` branches. Draft guidance remained explicitly marked for review. |
| 22 October 2026 | Corrected 16 documentation/generation findings from comparison with source and model references. | SU2 `7862504f9d` and website `be90edc`, linked above. Solver limitations were documented, not silently changed. |
| 22 October 2026 | Clarified page roles, introduced permanent checks and made local/CI builds use a shared entry point and locked Ruby/gem dependencies. | SU2 `a14b4df35d`, website `602189b` and the successful milestone workflows linked above. |
| 22 October 2026 | Added this maintained status page under Developer Docs and contributor guidance for keeping the record current. Fixed report self-link ordering with a reproducing regression check. | This page's source history. The current build manifest identifies the revisions that include it. |
| 22 October 2026 | Repaired all 180 historical link findings, eight template names and the obsolete structural output keyword; reviewed nine setup guides and four model pages; added a hyphenated-enum regression check. | Review scope and primary references above; current build manifest and check results identify this batch. Eleven tooling tests; CFD validation remains pending. |
| 23 October 2026 | Added the V&V case "NACA 0021 Beyond Stall: Hybrid RANS/LES" (SA and SST DDES/EDDES against experiment and published DES). Described the SST DES variants and the new `SA_EDDES` cell size of SU2 PR #2150 on the hybrid RANS/LES, turbulence model and unsteady setup pages. | Website `7279f4d`. Page content checked against PR #2150 at `e9e85686ad`. The V&V data come from local runs; the flow picture is still to be added. |

## Keeping this page current

Edit [Documentation-Status.md in the SU2 repository](https://github.com/rois1995/SU2/blob/docs_site/Docs/website/pages/Documentation-Status.md). The build copies it into the site and adds it under **Developer Docs**. Do not create a second maintained copy in the website repository.

For each completed batch:

1. Update the status date and the relevant implemented/backlog row. State what remains unreviewed or unsupported.
2. Add a dated history entry with source and website commit or pull-request links, the checks performed, and any numerical evidence actually obtained.
3. Remove resolved historical findings individually. Keep current generated counts in the generated reports; label any numbers quoted here with their milestone and date.
4. Run the shared documentation build and the relevant checks. Use **Docs Website** with the reviewed website commit SHA when publishing website-only changes.

See [Contribute](../contribute/) for page ownership and editing steps. Build commands and contributor details are in [Docs/website/README.md](https://github.com/rois1995/SU2/blob/docs_site/Docs/website/README.md).

Maintained against the two `docs_site` branches. Current source revisions and check totals are recorded in this build’s manifest and check report.
