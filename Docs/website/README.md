# SU2 website: pages kept with the code

This folder holds the parts of the SU2 website (https://su2code.github.io) that must change together with the code.
The workflow `.github/workflows/docs-website.yml` builds the website and publishes it on the GitHub Pages of this repository.

## What the build does

1. It gets the website. By default this is `rois1995/su2code.github.io`, branch `docs_site`. When you start the workflow by hand, you can choose another repository and branch.
2. It copies every page in `pages/` into the user guide (`_docs_v7/`). A page with the same name as a website page replaces it.
3. It copies the Markdown pages kept in the code folders (for example `Docs/Aeroacoustics.md`) into the user guide.
4. It makes the **Option Reference** from `Common/src/CConfig.cpp`. Each option shows:
   - its description, type, allowed values and default value;
   - the comment of `config_template.cfg`, when it adds something;
   - the website pages and tutorials that mention it;
   - up to three test case files that use it, with their values;
   - a link to the code.
5. It makes the **Docs Check** page. This page lists:
   - options that the docs, tutorials and V&V pages name, but that do not exist in the code;
   - options with no description, a very short one, or no example.
6. It builds the website with Jekyll, and the **Code Reference** with Doxygen.

## How to change the docs

Pages that explain the code (developer pages) live here, so they change in the same PR as the code.
User guide pages and tutorials live in the website repository.

- **New option:** write the comment above the `add...Option` call in `CConfig.cpp`, in the usual form:
  `/*!\brief NAME \n DESCRIPTION: what it does \n DEFAULT: value \ingroup Config*/`.
  The Option Reference updates by itself. Also add the option to `config_template.cfg`.
- **New or changed page:** add or edit a Markdown file in `pages/`. The file must start with a front matter: a block
  between two `---` lines that holds the page settings. It needs `title` and `permalink: /docs_v7/<File-Name>/`.
  Pages named `Adding-...` appear under "How to Add..." in the menu. Other new pages appear at the end of
  "Developer Docs".
- **Images:** put them in `pages/img/`. Link to them as `../../docs_files/website/<name>.png`.

## Run it on your computer

```
git clone --depth 1 https://github.com/su2code/su2code.github.io.git site
python3 Docs/website/tools/build_site.py --su2 . --site site
cd site && bundle install && bundle exec jekyll serve
```

Then open http://localhost:4000. For the Code Reference, run `doxygen doxyfile` in the SU2 folder and open `Docs/html/index.html`.
