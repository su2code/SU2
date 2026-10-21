# SU2 website: pages kept with the code

This folder holds the parts of the SU2 website (https://su2code.github.io) that must change together with the code.
The workflow `.github/workflows/docs-website.yml` builds the website and publishes it on the GitHub Pages of this repository.

## What the build does

1. Gets the website (by default `rois1995/su2code.github.io`, branch `docs_site`; other ones can be chosen when the workflow is started by hand).
2. Copies every page in `pages/` into the user guide (`_docs_v7/`). A page with the same name as a website page replaces it.
3. Makes the **Option Reference** from `Common/src/CConfig.cpp`. Each option gets its description, type, allowed values, default value and a link to the code.
4. Makes the **Docs Check** page: options named in the docs, tutorials and V&V pages that do not exist in the code.
5. Builds the website with Jekyll and the **Code Reference** with Doxygen.

## How to change the docs

Pages that explain the code (developer pages) live here, so they change in the same PR as the code.
User guide pages and tutorials live in the website repository.

- **New option:** write the comment above the `add...Option` call in `CConfig.cpp`. Use the usual form:
  `/*!\brief NAME \n DESCRIPTION: what it does \n DEFAULT: value \ingroup Config*/`.
  The Option Reference updates by itself. Also add the option to `config_template.cfg`.
- **New or changed page:** add or edit a Markdown file in `pages/`. It needs a front matter with `title` and
  `permalink: /docs_v7/<File-Name>/`. New pages appear at the end of "Developer Docs" in the menu.
- **Images:** put them in `pages/img/`. Use them as `../../docs_files/website/<name>.png`.

## Run it on your computer

```
git clone --depth 1 https://github.com/su2code/su2code.github.io.git site
python3 Docs/website/tools/build_site.py --su2 . --site site
cd site && bundle install && bundle exec jekyll serve
```

Open http://localhost:4000. Run `doxygen doxyfile` in the SU2 folder for the Code Reference (`Docs/html/index.html`).
