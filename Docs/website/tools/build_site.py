#!/usr/bin/env python3

## \file build_site.py
#  \brief Prepare a checkout of the SU2 website (su2code.github.io) for a Jekyll build:
#         copy the pages kept in this repository, generate the configuration option
#         reference from Common/src/CConfig.cpp and check the options named in the docs.
#  \version 8.5.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import argparse
import datetime
import glob
import json
import os
import re
import shutil
import subprocess
import sys

DOCS = "_docs_v7"
URL = "docs_v7"
NAV = os.path.join("_data", "docs_v7.yml")
REF_PREFIX = "Option-Reference"

# Order of the "How to Add..." pages in the menu, from small to large changes; other "Adding-" pages go last.
HOWTO_ORDER = ["Adding-a-Config-Option", "Adding-Output-Fields", "Adding-a-Boundary-Condition",
               "Adding-a-Numerical-Scheme", "Adding-a-Fluid-Model", "Adding-a-Turbulence-Model",
               "Adding-a-Transition-Model", "Adding-a-Solver"]

# Folders of the website that are only downloads; the test site links to the real site instead.
HEAVY_FOLDERS = ["documents", "download", "emails"]

# Upper case words followed by "=" in the docs that are not configuration options
# (keywords of the SU2 mesh format, etc.).
NOT_OPTIONS = {"NDIME", "NELEM", "NPOIN", "NMARK", "NZONE", "IZONE", "MARKER_TAG", "MARKER_ELEMS", "FFD_NBOX",
               "FFD_NLEVEL", "FFD_TAG", "FFD_LEVEL", "FFD_DEGREE_I", "FFD_DEGREE_J", "FFD_DEGREE_K",
               "FFD_PARENTS", "FFD_CHILDREN", "FFD_CORNER_POINTS", "FFD_CONTROL_POINTS", "FFD_SURFACE_POINTS",
               "NMARK_FFD", "NPERIODIC", "PERIODIC_INDEX", "NELEM_BOUND", "VALUE", "OPTION",
               # keys of the Python FSI configuration file (Unsteady_FSI_Python tutorial)
               "IMPOSED_MODES", "IMPOSED_PARAMETERS"}

TYPE_NAMES = {
  "Double": "real number", "Bool": "YES / NO", "String": "string", "UnsignedLong": "integer",
  "UnsignedShort": "integer", "Long": "integer", "DoubleList": "list of real numbers",
  "DoubleArray": "fixed-size list of real numbers", "StringList": "list of strings",
  "UShortList": "list of integers", "ShortList": "list of integers", "ULongList": "list of integers",
  "ULongArray": "fixed-size list of integers", "UShortArray": "fixed-size list of integers",
  "StringDoubleList": "list of (string, real number) pairs", "Python": "used by the Python scripts only",
}


def slugify(text):
  return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def split_args(text):
  """Split the arguments of a C++ call at the commas that are not inside brackets."""
  args, depth, start = [], 0, 0
  for i, c in enumerate(text):
    if c in "([{<":
      depth += 1
    elif c in ")]}>":
      depth -= 1
    elif c == "," and depth == 0:
      args.append(text[start:i].strip())
      start = i + 1
  args.append(text[start:].strip())
  return args


def call_body(text, open_paren):
  """Return the text between the parenthesis at open_paren and its matching parenthesis."""
  depth = 0
  for i in range(open_paren, len(text)):
    if text[i] == "(":
      depth += 1
    elif text[i] == ")":
      depth -= 1
      if depth == 0:
        return text[open_paren + 1:i], i
  return text[open_paren + 1:], len(text)


def read_maps(su2):
  """Map name -> list of (key, C++ value) for every MapType<std::string, ...> in Common/include."""
  maps = {}
  pattern = re.compile(r"MapType\s*<\s*std::string\s*,[^>]*>\s+(\w+)\s*=\s*\{(.*?)\};", re.S)
  for path in glob.glob(os.path.join(su2, "Common", "include", "**", "*.*pp"), recursive=True):
    with open(path, encoding="utf-8", errors="replace") as f:
      text = f.read()
    for m in pattern.finditer(text):
      pairs = re.findall(r'MakePair\(\s*"([^"]+)"\s*,\s*([^)]*?)\s*\)', m.group(2))
      maps[m.group(1)] = pairs
  return maps


def clean_comment(raw):
  """Turn a /*! ... */ comment of CConfig.cpp into (description, allowed values text, default text)."""
  text = raw.replace("\\ingroup Config", "")
  text = re.sub(r"\\brief\s+[A-Z0-9_]+\s*", "", text)
  text = text.replace("\\brief", "")
  text = re.sub(r"\\link\s+(\w+)\s*\\endlink", r"\1", text)
  text = re.sub(r"\s*\*\s*", " ", text) if "\n" in text else text
  parts = {"DESCRIPTION": "", "OPTIONS": "", "DEFAULT": ""}
  key = "DESCRIPTION"
  for chunk in re.split(r"\\n", text):
    chunk = chunk.strip()
    m = re.match(r"(DESCRIPTION|Options|OPTIONS|DEFAULT|Default)\s*:\s*(.*)", chunk, re.S)
    if m:
      key = m.group(1).upper()
      chunk = m.group(2)
    if chunk:
      parts[key] = (parts[key] + " " + chunk).strip()
  return parts["DESCRIPTION"], parts["OPTIONS"], parts["DEFAULT"]


def read_options(su2):
  """Read every option registered in CConfig::SetConfig_Options, in order, with its category."""
  with open(os.path.join(su2, "Common", "src", "CConfig.cpp"), encoding="utf-8") as f:
    text = f.read()

  # Option comments start with "/*!" or "/*"; only the ones with a description or a \brief are kept,
  # so that plain code comments such as "/*--- ... ---*/" do not become descriptions.
  token = re.compile(r"/\*!?(.*?)\*/|\badd(\w*)Option\s*\(", re.S)
  options, category, comment = [], "General", None
  pos = 0
  while True:
    m = token.search(text, pos)
    if not m:
      break
    if m.group(1) is not None:
      body = m.group(1)
      cat = re.search(r"\\par\s+CONFIG_CATEGORY\s*:\s*(.*?)\s*\\ingroup", body, re.S)
      if cat:
        category, comment = " ".join(cat.group(1).split()), None
      elif "DESCRIPTION" in body or "\\brief" in body:
        comment = body
      pos = m.end()
      continue
    args_text, end = call_body(text, m.end() - 1)
    pos = end + 1
    args = split_args(args_text)
    name = re.fullmatch(r'"([A-Z0-9_]+)"', args[0])
    if not name:
      continue
    line = text.count("\n", 0, m.start()) + 1
    options.append({"name": name.group(1), "kind": m.group(2), "args": args, "category": category,
                    "comment": comment, "line": line})
    comment = None
  return options


def read_template(su2):
  """Option name -> (comment lines above it, example line) from config_template.cfg."""
  examples = {}
  with open(os.path.join(su2, "config_template.cfg"), encoding="utf-8") as f:
    lines = f.read().splitlines()
  notes = []
  for line in lines:
    m = re.match(r"^\s*%?\s*([A-Z][A-Z0-9_]+)\s*=(.*)$", line)
    if m and not line.lstrip().startswith("%%"):
      examples.setdefault(m.group(1), (list(notes), line.strip()))
      notes = []
    elif line.lstrip().startswith("%") and not line.lstrip().startswith("% ---") and line.strip("% -"):
      notes.append(line.strip().lstrip("%").strip())
    else:
      notes = []
  return examples


def read_testcases(su2):
  """Option name -> list of (cfg path relative to the SU2 root, value) from TestCases/**/*.cfg."""
  usage = {}
  setting = re.compile(r"^\s*([A-Z][A-Z0-9_]+)\s*=\s*([^%]*)")
  for path in sorted(glob.glob(os.path.join(su2, "TestCases", "**", "*.cfg"), recursive=True)):
    rel = os.path.relpath(path, su2)
    with open(path, encoding="utf-8", errors="replace") as f:
      for line in f:
        m = setting.match(line)
        if m:
          usage.setdefault(m.group(1), []).append((rel, " ".join(m.group(2).split())))
  return usage


def pick_examples(uses, count=3):
  """Up to count uses with different values, so that the examples show the range of settings."""
  picked, values = [], set()
  for rel, value in uses:
    if value and value not in values:
      picked.append((rel, value))
      values.add(value)
    if len(picked) == count:
      break
  return picked


def read_page_mentions(site, known, skip):
  """Option name -> list of (title, url) of the user guide pages and tutorials that name the option."""
  mentions = {}
  files = sorted(glob.glob(os.path.join(site, DOCS, "*.md"))) + \
          sorted(glob.glob(os.path.join(site, "_tutorials", "**", "*.md"), recursive=True))
  for path in files:
    base = os.path.splitext(os.path.basename(path))[0]
    if base.startswith(REF_PREFIX) or base in skip or "unfinished" in path:
      continue
    with open(path, encoding="utf-8", errors="replace") as f:
      text = f.read()
    title = re.search(r"^title:\s*(.+)$", text, re.M)
    link = re.search(r"^permalink:\s*(\S+)$", text, re.M)
    if not (title and link):
      continue
    title = title.group(1).strip().strip("\"'")
    for name in set(re.findall(r"\b[A-Z][A-Z0-9_]{2,}\b", text)) & known:
      in_guide = path.startswith(os.path.join(site, DOCS))
      mentions.setdefault(name, []).append((title, link.group(1).strip("\"'"), in_guide))
  for name in mentions:
    # User guide pages first, then tutorials.
    mentions[name].sort(key=lambda m: (not m[2], m[0]))
  return mentions


def describe_type(opt, maps):
  """Return (type text, list of allowed values, default text)."""
  kind, args = opt["kind"], opt["args"]
  values, default = [], ""
  if kind in ("Enum", "EnumList"):
    map_name = args[2] if kind == "Enum" else args[3]
    values = [k for k, _ in maps.get(map_name, [])]
    if kind == "Enum" and len(args) > 3:
      default = next((k for k, v in maps.get(map_name, []) if v == args[3]), args[3])
    return ("one of the values below" if kind == "Enum" else "list of the values below"), values, default
  if kind == "Bool" and len(args) > 2:
    default = {"true": "YES", "false": "NO"}.get(args[2], args[2])
  elif kind in ("Double", "String", "UnsignedLong", "UnsignedShort", "Long") and len(args) > 2:
    default = args[2].strip('"')
  elif kind in ("DoubleArray", "UShortArray", "ULongArray") and len(args) > 3:
    default = args[3]
  return TYPE_NAMES.get(kind, "special format, see config_template.cfg"), values, default


def normalize(text):
  return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()


def same_text(a, b):
  """True if b says nothing more than a (same words, or contained in a)."""
  return normalize(a) == normalize(b) or normalize(b) in normalize(a)


def write_reference(su2, site, options, maps, commit, repo, skip):
  examples = read_template(su2)
  usage = read_testcases(su2)
  mentions = read_page_mentions(site, {o["name"] for o in options}, skip)
  stamp = "Generated from SU2 commit `{}` on {}.".format(commit, datetime.date.today().isoformat())
  blob = "https://github.com/{}/blob/{}/".format(repo, commit)
  source = blob + "Common/src/CConfig.cpp"

  categories = {}
  for opt in options:
    categories.setdefault(opt["category"], []).append(opt)

  pages = []
  for cat, opts in categories.items():
    slug = "{}-{}".format(REF_PREFIX, slugify(cat))
    for opt in opts:
      opt["page"] = slug
    out = ["---", 'title: "{}"'.format(cat.replace('"', "'")), "permalink: /{}/{}/".format(URL, slug), "---",
           "", "{% raw %}", "", "> {} Do not edit this page. Edit the comment above the option in "
           "[`Common/src/CConfig.cpp`]({}) instead.".format(stamp, source), ""]
    for opt in opts:
      desc, opt_text, def_text = clean_comment(opt["comment"]) if opt["comment"] else ("", "", "")
      notes, tmpl_line = examples.get(opt["name"], ([], ""))
      notes = " ".join(notes)
      type_text, values, default = describe_type(opt, maps)
      uses = usage.get(opt["name"], [])
      pages_of = mentions.get(opt["name"], [])
      opt.update({"desc_words": len(desc.split()), "has_notes": bool(notes), "uses": len(uses),
                  "mentions": len(pages_of)})

      out.append("### {} {{#{}}}".format(opt["name"], opt["name"]))
      out.append("")
      if desc:
        out.append(desc)
      elif notes:
        out.append(notes + " *(from config_template.cfg)*")
      else:
        out.append("*No description in the code yet.*")
      out.append("")
      if desc and notes and not same_text(desc, notes):
        out.append("> **More in config_template.cfg:** " + notes)
        out.append("")
      out.append("- **Type:** {}".format(type_text))
      if values:
        out.append("- **Values:** " + ", ".join("`{}`".format(v) for v in values))
      elif opt_text and not opt_text.startswith(tuple(maps)):
        out.append("- **Values:** {}".format(opt_text))
      if default or def_text:
        out.append("- **Default:** `{}`".format(default if default else def_text))
      if tmpl_line:
        out.append("- **In config_template.cfg:** `{}`".format(tmpl_line))
      else:
        out.append("- **Not in config_template.cfg.**")
      if pages_of:
        out.append("- **Mentioned in:** " + ", ".join(
          "[{}]({})".format(title, url) for title, url, _ in pages_of[:6]) +
          (" and {} more".format(len(pages_of) - 6) if len(pages_of) > 6 else ""))
      if uses:
        files = len({rel for rel, _ in uses})
        out.append("- **Used in {} test case file{}**, for example: ".format(files, "s" if files > 1 else "") +
                   ", ".join("[`{}`]({}{}) (`{}`)".format(os.path.basename(rel), blob, rel, value or "empty")
                             for rel, value in pick_examples(uses)))
      out.append("- **Code:** [`CConfig.cpp` line {}]({}#L{})".format(opt["line"], source, opt["line"]))
      out.append("")
    out += ["{% endraw %}", ""]
    with open(os.path.join(site, DOCS, slug + ".md"), "w", encoding="utf-8") as f:
      f.write("\n".join(out))
    pages.append((slug, cat, len(opts)))

  missing_desc = sum(1 for o in options if not o["desc_words"])
  missing_both = sum(1 for o in options if not o["desc_words"] and not o["has_notes"])
  missing_tmpl = sum(1 for o in options if o["name"] not in examples)
  out = ["---", "title: Option Reference", "permalink: /{}/{}/".format(URL, REF_PREFIX), "---", "",
         "{% raw %}", "",
         "This reference lists every option of the SU2 configuration file. It is made from the source code, "
         "so it always matches the code. {}".format(stamp), "",
         "Each option shows the description in the code, the text of `config_template.cfg` when it adds "
         "something, the pages of this website that mention it, and real test case files that use it.", "",
         "- Options: {}".format(len(options)),
         "- Options without a description in the code: {} ({} without any text at all)".format(
           missing_desc, missing_both),
         "- Options not in `config_template.cfg`: {}".format(missing_tmpl),
         "- The [Docs Check](../Docs-Check/) page lists the options that need better documentation.", "",
         "## Groups", ""]
  out += ["- [{}](../{}/) ({} options)".format(cat, slug, n) for slug, cat, n in pages]
  out += ["", "## All options A-Z", ""]
  for opt in sorted(options, key=lambda o: o["name"]):
    out.append("[`{0}`](../{1}/#{0}) ".format(opt["name"], opt["page"]))
  out += ["", "{% endraw %}", ""]
  with open(os.path.join(site, DOCS, REF_PREFIX + ".md"), "w", encoding="utf-8") as f:
    f.write("\n".join(out))
  return [REF_PREFIX] + [p[0] for p in pages]


def read_code_names(su2):
  """Upper case identifiers of the C++ code (enums, constants), to tell them apart from options."""
  names = set()
  for folder in ("Common", "SU2_CFD"):
    for path in glob.glob(os.path.join(su2, folder, "**", "*.[ch]pp"), recursive=True):
      with open(path, encoding="utf-8", errors="replace") as f:
        names.update(re.findall(r"\b[A-Z][A-Z0-9_]{2,}\b", f.read()))
  return names


def check_docs(site, known, code_names=frozenset()):
  """Find options in the docs, tutorials and V&V pages that the code does not know.
  Only code blocks that look like an SU2 configuration file (at least half known options)
  and inline code such as `OPTION= VALUE` are checked."""
  found = []
  setting = re.compile(r"^\s*%?\s*([A-Z][A-Z0-9_]{2,})\s*=")
  files = sorted(glob.glob(os.path.join(site, DOCS, "*.md")) +
                 glob.glob(os.path.join(site, "_tutorials", "**", "*.md"), recursive=True) +
                 glob.glob(os.path.join(site, "_vandv", "*.md")))
  for path in files:
    if os.path.basename(path).startswith(REF_PREFIX):
      continue
    with open(path, encoding="utf-8", errors="replace") as f:
      lines = f.read().splitlines()
    hits, block, code = set(), None, False
    for num, line in enumerate(lines, 1):
      if line.lstrip().startswith("```"):
        if block is None:
          block = []
          # Blocks of source code (C++, Python, shell) are not configuration files.
          code = line.strip("` \t").lower() in ("cpp", "c++", "c", "hpp", "python", "py", "bash", "sh", "shell")
        else:
          names = [(n, m.group(1)) for n, m in block if m]
          if not code and names and 2 * sum(name in known for _, name in names) >= len(names):
            hits.update((n, name) for n, name in names)
          block = None
        continue
      if block is not None:
        block.append((num, setting.match(line)))
      else:
        # Inline code such as `MAXNVAR = 4` can be C++; only names unknown to the code are reported.
        hits.update((num, name) for name in re.findall(r"`\s*([A-Z][A-Z0-9_]{2,})\s*=", line)
                    if name not in code_names)
    for num, name in sorted(hits):
      if name not in known and name not in NOT_OPTIONS:
        found.append((os.path.relpath(path, site), num, name, lines[num - 1].strip()))
  return found


# Output classes, in the order of the page, with the solvers that use them.
OUTPUT_CLASSES = [
  ("COutput", "All solvers"), ("CFVMOutput", "All finite volume solvers"),
  ("CFlowOutput", "All flow solvers (compressible, incompressible, NEMO), incl. turbulence and species"),
  ("CFlowCompOutput", "Compressible flow"), ("CFlowIncOutput", "Incompressible flow"),
  ("CNEMOCompOutput", "Non-equilibrium flow (NEMO)"), ("CFlowCompFEMOutput", "Finite element flow (DG-FEM)"),
  ("CTurboOutput", "Turbomachinery"), ("CHeatOutput", "Heat equation"), ("CElasticityOutput", "Structural analysis"),
  ("CAdjFlowOutput", "All adjoint flow solvers"), ("CAdjFlowCompOutput", "Compressible adjoint"),
  ("CAdjFlowIncOutput", "Incompressible adjoint"), ("CAdjHeatOutput", "Heat adjoint"),
  ("CAdjElasticityOutput", "Structural adjoint"), ("CMultizoneOutput", "Multizone problems"),
  ("CMeshOutput", "Mesh conversion (SU2_DEF, SU2_SOL)"), ("CBaselineOutput", "Solution conversion (SU2_SOL)"),
]


def literal(arg):
  """The text of a C++ string literal, or None if the argument is an expression."""
  m = re.fullmatch(r'"((?:[^"\\]|\\.)*)"', arg.strip())
  return m.group(1) if m else None


def read_outputs(su2):
  """Output class -> list of (kind, name, group, description, line) of the fields it registers."""
  calls = {"AddHistoryOutput": ("history", 3, 4), "AddHistoryOutputPerSurface": ("history, per surface", 3, None),
           "AddVolumeOutput": ("volume", 2, 3)}
  token = re.compile(r"\b(" + "|".join(calls) + r")\s*\(")
  fields = {}
  for path in sorted(glob.glob(os.path.join(su2, "SU2_CFD", "src", "output", "*.cpp"))):
    cls = os.path.splitext(os.path.basename(path))[0]
    with open(path, encoding="utf-8", errors="replace") as f:
      text = f.read()
    seen = set()
    for m in token.finditer(text):
      body, _ = call_body(text, m.end() - 1)
      args = split_args(body)
      kind, igroup, idesc = calls[m.group(1)]
      name = literal(args[0])
      shown = name if name is not None else " ".join(args[0].split())
      if (kind, shown) in seen:
        continue
      seen.add((kind, shown))
      group = literal(args[igroup]) if len(args) > igroup else None
      desc = literal(args[idesc]) if idesc is not None and len(args) > idesc else None
      fields.setdefault(cls, []).append({"kind": kind, "name": shown, "literal": name is not None,
                                         "group": group if group is not None else "",
                                         "desc": desc or "", "line": text.count("\n", 0, m.start()) + 1})
  return fields


def cell(text):
  return text.replace("|", "\\|").replace("\n", " ")


def write_output_page(site, fields, commit, blob):
  name = "Output-Fields-Reference"
  out = ["---", "title: Output Fields Reference", "permalink: /{}/{}/".format(URL, name), "---", "", "{% raw %}", "",
         "This page lists the screen, history and volume output fields that each output class of SU2 can write. "
         "It is made from the `AddHistoryOutput` and `AddVolumeOutput` calls in `SU2_CFD/src/output/`. "
         "Generated from SU2 commit `{}`.".format(commit), "",
         "- A field is written only when the case needs it (for example turbulence fields only with a turbulence "
         "model). For the exact list of your case, run a dry run: `SU2_CFD -d your_config.cfg`. It prints "
         "\"Available screen/history output fields\" and \"Available volume output fields\".",
         "- How to select fields: [Custom Output](../Custom-Output/). How to add one: "
         "[Adding Output Fields](../Adding-Output-Fields/).",
         "- Names in *italics* are built at run time (for example one field per species or per marker).", ""]
  order = [c for c, _ in OUTPUT_CLASSES if c in fields] + sorted(c for c in fields if c not in dict(OUTPUT_CLASSES))
  used = dict(OUTPUT_CLASSES)
  for cls in order:
    path = "SU2_CFD/src/output/{}.cpp".format(cls)
    out += ["## {} {{#{}}}".format(cls, cls), "", "Used by: {}. Source: [`{}`]({}{}).".format(
      used.get(cls, "see the code"), path, blob, path), ""]
    for kind in ("history", "history, per surface", "volume"):
      rows = [f for f in fields[cls] if f["kind"] == kind]
      if not rows:
        continue
      out += ["**{} fields ({})**".format(kind.capitalize(), len(rows)), "", "| Field | Group | Description |",
              "|---|---|---|"]
      for f in rows:
        shown = "`{}`".format(f["name"]) if f["literal"] else "*{}*".format(cell(f["name"]))
        out.append("| [{}]({}{}#L{}) | {} | {} |".format(shown, blob, path, f["line"], cell(f["group"]),
                                                         cell(f["desc"]) or "*no description*"))
      out.append("")
  out += ["{% endraw %}", ""]
  with open(os.path.join(site, DOCS, name + ".md"), "w", encoding="utf-8") as f:
    f.write("\n".join(out))
  return name


def find_stale_pages(su2, pages_dir):
  """Developer pages whose named code files changed since the commit in their "Last checked" line."""
  checked = re.compile(r"Last checked against SU2 \w+ `([0-9a-f]{7,40})`")
  named = re.compile(r"`((?:Common|SU2_CFD|SU2_DEF|SU2_DOT|SU2_GEO|SU2_SOL|SU2_PY|UnitTests|TestCases|externals|"
                     r"meson_scripts|Docs|\.github)/[^`\s:]+|config_template\.cfg|meson\.build|meson_options\.txt|"
                     r"doxyfile)")
  report = []
  for path in sorted(glob.glob(os.path.join(pages_dir, "*.md"))):
    with open(path, encoding="utf-8") as f:
      text = f.read()
    page = os.path.splitext(os.path.basename(path))[0]
    m = checked.search(text)
    if not m:
      report.append((page, None, "no \"Last checked against\" line", []))
      continue
    files = sorted({p for p in named.findall(text) if os.path.exists(os.path.join(su2, p))})
    try:
      out = subprocess.check_output(["git", "-C", su2, "log", "--format=%h", "--name-only", m.group(1) + "..HEAD",
                                     "--"] + files, text=True, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.CalledProcessError):
      report.append((page, m.group(1), "commit not found in the history", []))
      continue
    counts = {}
    for line in out.splitlines():
      if line in files:
        counts[line] = counts.get(line, 0) + 1
    report.append((page, m.group(1), "", sorted(counts.items(), key=lambda c: (-c[1], c[0]))))
  return report


def doc_coverage(su2):
  """Classes in the headers of Common and SU2_CFD: total, with a Doxygen comment, and in a Doxygen group.
  Approximate: found with a regular expression, not with a C++ parser."""
  cls = re.compile(r"^\s*(?:template\s*<[^;{]*>\s*)?class\s+(?:\w+\s+)?([A-Z]\w*)\b[^;]*?(?:\{|:[^:])", re.M)
  stats = {}
  for folder in ("Common/include", "SU2_CFD/include"):
    for path in sorted(glob.glob(os.path.join(su2, folder, "**", "*.hpp"), recursive=True)):
      with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
      rel = os.path.relpath(path, su2)
      seen = set()
      for m in cls.finditer(text):
        # Skip repeated declarations (specializations) and all-capital helper types of template code.
        if m.group(1) in seen or m.group(1).isupper():
          continue
        seen.add(m.group(1))
        before = text[:m.start()].rstrip()
        # The comment must end just before the class (template lines are part of the match).
        comment = before[before.rfind("/*"):] if before.endswith("*/") else ""
        key = os.path.dirname(rel)
        total, documented, grouped, missing = stats.get(key, (0, 0, 0, []))
        stats[key] = (total + 1, documented + bool(comment), grouped + ("\\ingroup" in comment),
                      missing + ([m.group(1)] if not comment else []))
  return stats


def read_doxygen_log(path):
  """Doxygen warnings: total, by kind and by folder."""
  kinds, folders, total = {}, {}, 0
  if not path or not os.path.isfile(path):
    return None
  with open(path, encoding="utf-8", errors="replace") as f:
    for line in f:
      m = re.match(r"(.*?):\d+: warning: (.*)", line)
      if not m:
        continue
      total += 1
      kind = re.sub(r"'[^']*'|\"[^\"]*\"|\b\w+::\w+(::\w+)*\b", "X", m.group(2))
      kind = re.sub(r"\(.*", "", kind).strip()[:80]
      kinds[kind] = kinds.get(kind, 0) + 1
      parts = re.split(r"[\\/]", m.group(1))
      top = next((i for i, p in enumerate(parts) if p in ("Common", "SU2_CFD", "SU2_DEF", "SU2_DOT", "SU2_GEO",
                                                          "SU2_SOL", "Docs")), None)
      folder = "/".join(parts[top:top + 2]) if top is not None else "other"
      folders[folder] = folders.get(folder, 0) + 1
  return total, kinds, folders


def quality_sections(stale, coverage, doxy):
  """Extra sections of the Docs Check page: stale pages and code documentation coverage."""
  out = ["", "## Developer pages that may be out of date", "",
         "Each developer page ends with \"Last checked against SU2 develop `<commit>`\". These pages name code "
         "files that changed since that commit. Read the changes, update the page if needed, then update the "
         "commit in its last line.", ""]
  flagged = [r for r in stale if r[2] or r[3]]
  out.append("Pages to review: {} of {}.".format(len(flagged), len(stale)))
  out.append("")
  for page, commit, problem, changed in flagged:
    out.append("- [{}](../{}/){}".format(page.replace("-", " "), page, ": " + problem if problem else ""))
    for path, n in changed:
      out.append("  - `{}`: {} commit{}".format(path, n, "s" if n > 1 else ""))
  total = sum(v[0] for v in coverage.values())
  documented = sum(v[1] for v in coverage.values())
  grouped = sum(v[2] for v in coverage.values())
  out += ["", "## Code documentation coverage", "",
          "Classes declared in the headers of `Common/include` and `SU2_CFD/include` (counted with a simple "
          "text search, so the numbers are approximate). `Docs/docmain.hpp` asks that every class has a Doxygen "
          "comment with `\\brief` and `\\ingroup`.", "",
          "- Classes: {}".format(total),
          "- With a Doxygen comment: {} ({:.0f}%)".format(documented, 100.0 * documented / max(total, 1)),
          "- In a Doxygen group (`\\ingroup`): {} ({:.0f}%)".format(grouped, 100.0 * grouped / max(total, 1)), "",
          "| Folder | Classes | With comment | In a group |", "|---|---|---|---|"]
  for folder, (t, d, g, _) in sorted(coverage.items()):
    out.append("| `{}` | {} | {} | {} |".format(folder, t, d, g))
  out += ["", '<details markdown="1"><summary>Classes without a Doxygen comment</summary>', ""]
  for folder, (_, _, _, missing) in sorted(coverage.items()):
    if missing:
      out.append("- `{}`: {}".format(folder, ", ".join("`{}`".format(c) for c in missing)))
  out += ["", "</details>", ""]
  if doxy:
    total, kinds, folders = doxy
    out += ["", "## Doxygen warnings", "", "Doxygen found {} problems in the code comments.".format(total), "",
            "| Kind | Count |", "|---|---|"]
    out += ["| {} | {} |".format(cell(k), n) for k, n in sorted(kinds.items(), key=lambda k: -k[1])[:12]]
    out += ["", "| Folder | Count |", "|---|---|"]
    out += ["| `{}` | {} |".format(k, n) for k, n in sorted(folders.items(), key=lambda k: -k[1])]
  return out


def write_check_page(site, found, commit, options, extra=()):
  out = ["---", "title: Docs Check", "permalink: /{}/Docs-Check/".format(URL), "---", "", "{% raw %}", "",
         "This page is made at each build. It lists what needs work in the documentation of SU2 commit `{}`."
         .format(commit), "",
         "## Options named in the docs that do not exist", "",
         "The option was renamed or removed, or it has a typo. Fix the page in the website repository.", "",
         "Found: {}".format(len(found)), ""]
  last = None
  for path, num, name, line in found:
    if path != last:
      out += ["", "### `{}`".format(path), ""]
      last = path
    out.append("- line {}: `{}`".format(num, name))

  def link(opt):
    return "[`{0}`](../{1}/#{0})".format(opt["name"], opt["page"])

  checks = [
    ("Options without any description",
     "No description in `Common/src/CConfig.cpp` and no comment in `config_template.cfg`.",
     lambda o: not o["desc_words"] and not o["has_notes"]),
    ("Options with a description only in config_template.cfg",
     "Move or copy the text into the comment above the option in `Common/src/CConfig.cpp`.",
     lambda o: not o["desc_words"] and o["has_notes"]),
    ("Options with a very short description",
     "The comment in `Common/src/CConfig.cpp` has 6 words or fewer.",
     lambda o: 0 < o["desc_words"] <= 6),
    ("Options with no example and no explanation",
     "No test case uses the option and no page of this website names it.",
     lambda o: not o["uses"] and not o["mentions"]),
  ]
  for title, text, test in checks:
    selected = [o for o in options if test(o)]
    out += ["", "## {} ({})".format(title, len(selected)), "", text, ""]
    groups = {}
    for opt in selected:
      groups.setdefault(opt["category"], []).append(link(opt))
    for cat, links in groups.items():
      out.append("- **{}:** {}".format(cat, ", ".join(links)))
  out += list(extra)
  out += ["", "{% endraw %}", ""]
  with open(os.path.join(site, DOCS, "Docs-Check.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(out))


def copy_code_notes(su2, site):
  """Copy the Markdown pages kept in the code folders (e.g. Docs/Aeroacoustics.md) into the user guide."""
  folders = ["Docs", "Common", "SU2_CFD", "SU2_DEF", "SU2_DOT", "SU2_GEO", "SU2_SOL", "SU2_PY"]
  names = []
  for folder in folders:
    for path in sorted(glob.glob(os.path.join(su2, folder, "**", "*.md"), recursive=True)):
      if os.path.join("Docs", "website") in path:
        continue
      with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
      heading = re.search(r"^#\s+(.+)$", text, re.M)
      title = heading.group(1).strip() if heading else os.path.basename(path)
      if heading:
        text = text[:heading.start()] + text[heading.end():]
      name = "Code-Notes-" + slugify(os.path.splitext(os.path.basename(path))[0])
      rel = os.path.relpath(path, su2)
      page = ["---", 'title: "{}"'.format(title.replace('"', "'")), "permalink: /{}/{}/".format(URL, name), "---",
              "", "> This page is kept in the SU2 code repository: `{}`.".format(rel), "",
              "{% raw %}", text.strip(), "{% endraw %}", ""]
      with open(os.path.join(site, DOCS, name + ".md"), "w", encoding="utf-8") as f:
        f.write("\n".join(page))
      names.append(name)
  return names


def add_to_section(nav, title, pages):
  """Append pages to a section of the menu, unless they are listed already."""
  sec = re.search(r"- title: {}\n  docs_v7:\n((?:  - .*\n)+)".format(re.escape(title)), nav)
  if not sec:
    return nav
  added = "".join("  - {}\n".format(p) for p in pages if "  - {}\n".format(p) not in sec.group(1))
  return nav[:sec.end(1)] + added + nav[sec.end(1):]


def patch_nav(site, ref_pages, extra_dev_pages, user_pages):
  path = os.path.join(site, NAV)
  with open(path, encoding="utf-8") as f:
    nav = f.read()
  section = "- title: Option Reference\n  docs_v7:\n" + "".join("  - {}\n".format(p) for p in ref_pages) + "\n"
  nav = nav.replace("- title: Developer Docs", section + "- title: Developer Docs", 1)
  # The "Adding-..." pages get their own section, after the developer docs.
  howto = sorted((p for p in extra_dev_pages if p.startswith("Adding-")),
                 key=lambda p: (HOWTO_ORDER.index(p) if p in HOWTO_ORDER else len(HOWTO_ORDER), p))
  if howto:
    dev = re.search(r"- title: Developer Docs\n  docs_v7:\n(?:  - .*\n)+", nav)
    block = "- title: How to Add...\n  docs_v7:\n" + "".join("  - {}\n".format(p) for p in howto)
    nav = nav[:dev.end()] + block + nav[dev.end():] if dev else nav
  nav = add_to_section(nav, "Developer Docs", [p for p in extra_dev_pages if p not in howto])
  nav = add_to_section(nav, "Users Guide", user_pages)
  with open(path, "w", encoding="utf-8") as f:
    f.write(nav)


UPSTREAM_SITE = "https://github.com/su2code/su2code.github.io.git"

DEBUG_JS = """// Debug only: marks in red the pages and menu sections changed on this test site.
(function () {
  var pages = %s;
  var sections = %s;
  function changed(path) {
    return pages.some(function (p) { return path.replace(/index\\.html$/, "").endsWith(p); });
  }
  function red(el) {
    el.style.color = "#c00";
    el.querySelectorAll("a").forEach(function (a) { a.style.color = "#c00"; });
  }
  // Menu items of changed pages.
  document.querySelectorAll("a.list-group-item").forEach(function (a) {
    if (changed(new URL(a.href, location.href).pathname)) { red(a); a.style.fontWeight = "bold"; a.dataset.changed = "1"; }
  });
  // Menu sections that are new or contain a changed page.
  document.querySelectorAll(".panel").forEach(function (panel) {
    var title = panel.querySelector(".panel-title");
    if (title && (sections.indexOf(title.textContent.trim()) >= 0 || panel.querySelector("a[data-changed]"))) {
      red(title);
    }
  });
  if (changed(location.pathname)) {
    document.querySelectorAll("h1").forEach(function (h) { h.style.color = "#c00"; });
  }
})();
"""


def changed_site_pages(site):
  """Permalinks of the website pages that differ from the official website (su2code/su2code.github.io)."""
  try:
    subprocess.check_call(["git", "-C", site, "fetch", "-q", "--depth=50", UPSTREAM_SITE, "master"],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
      base = subprocess.check_output(["git", "-C", site, "merge-base", "FETCH_HEAD", "HEAD"], text=True).strip()
    except subprocess.CalledProcessError:
      base = "FETCH_HEAD"
    files = subprocess.check_output(["git", "-C", site, "diff", "--name-only", base, "HEAD"], text=True).split()
  except (OSError, subprocess.CalledProcessError):
    print("Could not compare the website with the official one; only the generated pages are marked.")
    return []
  urls = []
  for rel in files:
    path = os.path.join(site, rel)
    if not rel.endswith(".md") or not os.path.isfile(path):
      continue
    with open(path, encoding="utf-8", errors="replace") as f:
      link = re.search(r"^permalink:\s*(\S+)$", f.read(), re.M)
    if link:
      urls.append(link.group(1).strip("\"'"))
  return urls


def mark_changes(site, urls, sections):
  """Debug: add a script to every page that colours the changed pages and menu sections in red."""
  os.makedirs(os.path.join(site, "js"), exist_ok=True)
  with open(os.path.join(site, "js", "debug-changes.js"), "w", encoding="utf-8") as f:
    f.write(DEBUG_JS % (json.dumps(sorted(set(urls))), json.dumps(sections)))
  layout = os.path.join(site, "_layouts", "default.html")
  with open(layout, encoding="utf-8") as f:
    text = f.read()
  tag = '<script src="{{ "/js/debug-changes.js" | relative_url }}"></script>\n'
  if tag not in text:
    with open(layout, "w", encoding="utf-8") as f:
      f.write(text.replace("</body>", tag + "</body>", 1))


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--su2", required=True, help="root of the SU2 source tree")
  parser.add_argument("--site", required=True, help="checkout of su2code.github.io, changed in place")
  parser.add_argument("--light", action="store_true", help="remove the download folders (test site)")
  parser.add_argument("--doxygen-log", help="file with the Doxygen warnings (WARN_LOGFILE)")
  parser.add_argument("--mark-changes", action="store_true",
                      help="debug: show in red the pages that differ from the official website")
  args = parser.parse_args()

  try:
    commit = subprocess.check_output(["git", "-C", args.su2, "rev-parse", "--short=10", "HEAD"], text=True).strip()
  except (OSError, subprocess.CalledProcessError):
    commit = "unknown"

  # Pages kept in this repository replace or extend the ones of the website.
  pages_dir = os.path.join(args.su2, "Docs", "website", "pages")
  dev_pages = []
  for path in sorted(glob.glob(os.path.join(pages_dir, "*.md"))):
    shutil.copy(path, os.path.join(args.site, DOCS))
    dev_pages.append(os.path.splitext(os.path.basename(path))[0])
  for folder in ("img", "files"):
    src = os.path.join(pages_dir, folder)
    if os.path.isdir(src):
      shutil.copytree(src, os.path.join(args.site, "docs_files", "website"), dirs_exist_ok=True)

  # Links to the code point to the repository that is built (a fork during tests).
  repo = os.environ.get("GITHUB_REPOSITORY", "su2code/SU2")

  user_pages = copy_code_notes(args.su2, args.site)
  maps = read_maps(args.su2)
  options = read_options(args.su2)
  ref_pages = write_reference(args.su2, args.site, options, maps, commit, repo, set(dev_pages))
  found = check_docs(args.site, {o["name"] for o in options}, read_code_names(args.su2))
  output_page = write_output_page(args.site, read_outputs(args.su2), commit,
                                  "https://github.com/{}/blob/{}/".format(repo, commit))
  extra = quality_sections(find_stale_pages(args.su2, pages_dir), doc_coverage(args.su2),
                           read_doxygen_log(args.doxygen_log))
  write_check_page(args.site, found, commit, options, extra)
  patch_nav(args.site, ref_pages, dev_pages + ["Docs-Check"], user_pages + [output_page])

  if args.mark_changes:
    generated = ref_pages + dev_pages + user_pages + [output_page, "Docs-Check"]
    urls = ["/{}/{}/".format(URL, p) for p in generated] + changed_site_pages(args.site)
    mark_changes(args.site, urls, ["Option Reference", "How to Add..."])

  if args.light:
    for folder in HEAVY_FOLDERS:
      shutil.rmtree(os.path.join(args.site, folder), ignore_errors=True)

  print("Options: {}, groups: {}, unknown options named in the docs: {}".format(
    len(options), len(ref_pages) - 1, len(found)))
  for path, num, name, _ in found:
    print("  {}:{}: {}".format(path, num, name))
  return 0


if __name__ == "__main__":
  sys.exit(main())
