# SU2 documentation website

The combined site uses two repositories. Developer pages and the generation/build/check tools live here; user pages, model explanations, setup guides, tutorials, V&V, navigation and search live in the website repository. Builds stage copies of both sources; the source checkouts are not overwritten.

## Where an explanation belongs

- **Setup guides:** how to choose settings for a class of problems, why they matter, convergence checks and troubleshooting. Link to worked cases instead of repeating their full instructions.
- **Tutorials:** complete worked cases, with meshes, configurations, commands and post-processing.
- **V&V:** numerical evidence, comparisons and their scope; distinguish implementation verification from physical validation.
- **Model pages:** assumptions, exact model variant, definitions and units, governing equations or precise equation references, and differences between the literature and SU2.
- **Developer pages:** source locations, control flow, extension steps and appropriate regression or unit checks.
- **Generated references:** option/output/error facts extracted from the source. Explain extraction limits rather than guessing defaults or applicability.

The website's **Using the Documentation** page explains this structure to readers. Keep general explanations in one place and cross-link them. In a model or developer explanation, identify the source file and checked commit; for a theoretical claim, cite the relevant paper or equation. A regression config establishes usage, not a universal numerical recommendation. Retain draft labels until the stated review is complete.

## Editing

Developer pages live in `pages/`, with `title` and `permalink: /docs_v7/<File-Name>/` front matter. `Adding-...` pages appear under **How to Add...**; other new pages appear under **Developer Docs**. Images in `pages/img/` are copied to `docs_files/website/`.

After checking a developer page against the implementation, update its `Last checked against SU2 develop <commit>` line. The Docs Check page flags named files changed since that commit; this is a review prompt, not a proof that the page is wrong or complete.

Option descriptions belong above their registration in `Common/src/CConfig.cpp`, using the existing Doxygen `DESCRIPTION:`, `OPTIONS:` and `DEFAULT:` conventions. Keep examples in `config_template.cfg` consistent with the parser. Defaults set elsewhere may not be extractable; a template example is not necessarily a default. The reference also links to relevant pages and regression configurations.

The generator copies these pages and Markdown notes from code folders into a staged website. It produces the option, output, error and migration references and the Docs Check page, then updates navigation. Do not edit these generated copies in the website repository.

## Local build

Requirements: Python 3.8 or newer, Git, Ruby **3.3.12** (the website's `.ruby-version`), Bundler **2.4.22**, Doxygen and Graphviz. Use an existing Ruby version manager to select the required Ruby; the system Ruby need not change. The website's tracked `Gemfile.lock` fixes gem versions.

With the two `docs_site` checkouts, run from the SU2 source root:

```sh
# Set these paths to your two existing checkouts/output folder.
WEBSITE_SOURCE=/path/to/su2code.github.io
DOCS_BUILD=/path/to/docs-build

gem install bundler -v 2.4.22
BUNDLE_GEMFILE="$WEBSITE_SOURCE/Gemfile" bundle _2.4.22_ install
python3 Docs/website/tools/build_docs.py --site "$WEBSITE_SOURCE" --build-dir "$DOCS_BUILD" --baseurl /SU2 --mark-changes
mkdir -p "$DOCS_BUILD/serve"
ln -sfn ../site "$DOCS_BUILD/serve/SU2"
python3 -m http.server 8765 --directory "$DOCS_BUILD/serve"
```

Open `http://localhost:8765/SU2/`. For a direct server without the symlink/prefix, build with `--baseurl ''`, serve `$DOCS_BUILD/site` and open `http://localhost:8765/`.

`--mark-changes` preserves the red highlighting of added/modified menu items, menu sections and page titles compared with the official site. Both CI workflows pass it. This debug comparison fetches upstream website history; if unavailable, generated/developer pages remain highlighted and the build reports the limitation.

The build directory must be outside both source trees. It must be empty on first use; subsequent runs replace only its managed `source`, `site` and `doxygen` subdirectories. Keep source files elsewhere. Omitted download folders link to the official site.

The same command builds Doxygen, stages the site, generates references, builds Jekyll, adjusts subpath links, and checks the assembled result. It runs the Python tests first and fails on new documentation problems. There is no separate local assembly recipe.

## Checks and existing debt

Run the small checks directly when changing tooling:

```sh
python3 -m unittest discover -s Docs/website/tools
```

The full build checks local links and anchors in current documentation, Tutorials and V&V; configuration option names and enum keywords in configuration blocks; descriptions for newly added options; search entries and searchable output-field names. It does not request external URLs or establish that numerical examples are converged. Solver-specific option combinations still require appropriate SU2 tests.

`tools/check-baseline.json` records pre-existing broken links, eight template names unknown to the C++ registration, and 24 options missing source descriptions, and historical enum-keyword findings. Existing link findings remain visible in the machine-readable `docs-check.json`; new findings fail the build. Remove an entry when its problem is fixed. Do not regenerate the baseline to approve new regressions. The generated Docs Check page also lists source-description and stale-page work.

Each assembled site includes:

- `build-info.json`: both source commit SHAs, whether either working tree was dirty, the gem lock hash, base path, build flags and tool versions.
- `docs-check.json`: new problems, checked-reference/search counts and existing link findings.

For a reproducible rebuild, check out the two recorded SHAs, select the recorded Ruby and locked gems, and use the recorded base path and flags. Dirty working trees are useful for previews but cannot be reconstructed from commit SHAs alone. Doxygen/Graphviz versions are recorded, not pinned; changes in those tools can change generated output. The date on the option reference comes from the source commits through `SOURCE_DATE_EPOCH`.

## CI and publishing

`.github/workflows/docs-website.yml` checks SU2 `docs_site` pushes and pull requests. Pull requests build an artifact without deploying. Successful pushes and manual runs publish the checked artifact to this fork's Pages site. Checkout references are resolved once; the manifest records the exact revisions used.

The website repository's `docs-check.yml` checks its own pushes/pull requests against SU2 `docs_site` and uploads a checked preview. Website-only changes therefore get checks without cross-repository write credentials. After those checks pass, publish from the SU2 **Docs Website** workflow, choosing the website commit SHA as `site_ref`:

```sh
gh workflow run docs-website.yml --repo rois1995/SU2 --ref docs_site -f site_repo=rois1995/su2code.github.io -f site_ref=<website-commit-sha>
```

To reproduce a historical deployment, select its SU2 revision as well. Normal branch builds follow the current website branch; use a commit SHA when reviewing or reproducing a particular pair of changes. Changes to workflow files require GitHub credentials allowed to update workflows; ordinary repository-only OAuth scopes may reject them.

Both workflows target the `rois1995` forks. Do not push or publish these development changes to the official repositories.
