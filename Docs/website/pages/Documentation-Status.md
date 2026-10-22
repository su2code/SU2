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
| Physical model explanations | Implemented; review needed | Dedicated pages cover [turbulence](../Turbulence-Models/), [transition](../Transition-Models/), [hybrid RANS/LES](../Hybrid-RANS-LES/) and [wall treatment](../Wall-Treatment/), including variants, implementation differences and known limits. Scientific completeness still needs expert review. |
| Case setup guidance | Nine drafts; review needed | Guides cover external RANS, incompressible flow, unsteady/DES, transition, turbomachinery, CHT, species/combustion, shape optimization and NICFD. Shared introductions explain evidence tags and review status. |
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

## What is still missing

The quantities below are a snapshot from the same milestone. The source list of accepted historical findings is [check-baseline.json](https://github.com/rois1995/SU2/blob/a14b4df35d3863fc79b48c781fc83b8cb222d912/Docs/website/tools/check-baseline.json); current link findings also appear in the generated check results.

| Priority | Remaining work | What completion requires |
| --- | --- | --- |
| Next | Repair the 180 historical link findings, including missing assets, old anchors and paths escaping the preview base URL. | Correct destinations or restore the intended assets, pass the checks, and remove each resolved baseline entry. Do not hide a new failure by regenerating the baseline. |
| Next | Investigate eight template names unknown to the C++ option registrations and the historical enum keyword `PARAVIEW_BINARY`. | Establish which parser owns each setting, correct obsolete examples where appropriate, and check the resulting configuration against the applicable SU2 version. |
| Next | Improve option descriptions: 24 registrations lacked source descriptions; approximately 296 descriptions were very short. Enum choices also need explanations. | Describe effect, units, applicability, restrictions, interactions and conditional defaults in the source, prioritizing options used by the setup guides. Check the generated reference afterwards. |
| Next | Review the nine draft setup guides and physical model explanations. | Record the reviewer, checked version, scope and evidence. Retain draft/review labels until that review is complete; do not generalize a short regression configuration into an engineering recommendation. |
| Then | Refresh older tutorials, version claims and contribution instructions as they are checked. | Confirm linked case files, accepted settings, run instructions and stated results for the declared version. Keep worked examples in Tutorials and numerical comparisons in V&V. |
| Then | Fill missing model explanations and keep known support limitations current. | Link exact variants and equation references to the SU2 implementation, define symbols/units, and distinguish verified behavior from proposed or unsupported capabilities. Extend the existing model pages where appropriate. |
| Then | Strengthen extraction and build checks when a new mismatch is demonstrated. | Add a focused reproducing test and fix its cause. Solver-specific configuration compatibility still needs appropriate SU2 checks; a general parser or framework rewrite is not an agreed requirement. |
| Then | Decide whether the remaining build tools need tighter version control. | Demonstrate that an unpinned tool causes a relevant difference, then document or pin it. Continue recording actual versions and both source SHAs. |
| Separate work | Verify and fix underlying solver defects noted during documentation review. | Reproduce the behavior, add an appropriate failing check, fix the implementation and update the affected documentation. A documentation caveat alone is not a solver repair. |
| Separate review | Propose the improvements to the official repositories. | Submit focused, reviewable changes with evidence and coordinated source/website revisions. Fork deployment does not establish upstream approval. |

## How to go further

1. **Keep the new checks passing.** Repair historical links and configuration examples in small groups; remove only the entries actually resolved. This gives a measurable reduction in existing debt.
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

## Keeping this page current

Edit [Documentation-Status.md in the SU2 repository](https://github.com/rois1995/SU2/blob/docs_site/Docs/website/pages/Documentation-Status.md). The build copies it into the site and adds it under **Developer Docs**. Do not create a second maintained copy in the website repository.

For each completed batch:

1. Update the status date and the relevant implemented/backlog row. State what remains unreviewed or unsupported.
2. Add a dated history entry with source and website commit or pull-request links, the checks performed, and any numerical evidence actually obtained.
3. Remove resolved historical findings individually. Keep current generated counts in the generated reports; label any numbers quoted here with their milestone and date.
4. Run the shared documentation build and the relevant checks. Use **Docs Website** with the reviewed website commit SHA when publishing website-only changes.

See [Contribute](../contribute/) for page ownership and editing steps. Build commands and contributor details are in [Docs/website/README.md](https://github.com/rois1995/SU2/blob/docs_site/Docs/website/README.md).

Last checked against SU2 develop `a14b4df35d`.
