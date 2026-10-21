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
import os
import re
import shutil
import subprocess
import sys

DOCS = "_docs_v7"
URL = "docs_v7"
NAV = os.path.join("_data", "docs_v7.yml")
REF_PREFIX = "Option-Reference"

# Folders of the website that are only downloads; the test site links to the real site instead.
HEAVY_FOLDERS = ["documents", "download", "emails"]

# Upper case words followed by "=" in the docs that are not configuration options
# (keywords of the SU2 mesh format, etc.).
NOT_OPTIONS = {"NDIME", "NELEM", "NPOIN", "NMARK", "NZONE", "IZONE", "MARKER_TAG", "MARKER_ELEMS", "FFD_NBOX",
               "FFD_NLEVEL", "FFD_TAG", "FFD_LEVEL", "FFD_DEGREE_I", "FFD_DEGREE_J", "FFD_DEGREE_K",
               "FFD_PARENTS", "FFD_CHILDREN", "FFD_CORNER_POINTS", "FFD_CONTROL_POINTS", "FFD_SURFACE_POINTS",
               "NMARK_FFD", "NPERIODIC", "PERIODIC_INDEX", "NELEM_BOUND", "VALUE", "OPTION"}

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

  token = re.compile(r"/\*!(.*?)\*/|\badd(\w*)Option\s*\(", re.S)
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
      else:
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


def write_reference(su2, site, options, maps, commit):
  examples = read_template(su2)
  stamp = "Generated from SU2 commit `{}` on {}.".format(commit, datetime.date.today().isoformat())
  source = "https://github.com/su2code/SU2/blob/{}/Common/src/CConfig.cpp".format(commit)

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
      type_text, values, default = describe_type(opt, maps)
      out.append("### {} {{#{}}}".format(opt["name"], opt["name"]))
      out.append("")
      out.append(desc if desc else "*No description in the code yet.*")
      out.append("")
      out.append("- **Type:** {}".format(type_text))
      if values:
        out.append("- **Values:** " + ", ".join("`{}`".format(v) for v in values))
      elif opt_text and not opt_text.startswith(tuple(maps)):
        out.append("- **Values:** {}".format(opt_text))
      if default or def_text:
        out.append("- **Default:** `{}`".format(default if default else def_text))
      if opt["name"] in examples:
        notes, line = examples[opt["name"]]
        out.append("- **In config_template.cfg:** `{}`".format(line))
      else:
        out.append("- **Not in config_template.cfg.**")
      out.append("- **Code:** [`CConfig.cpp` line {}]({}#L{})".format(opt["line"], source, opt["line"]))
      out.append("")
    out += ["{% endraw %}", ""]
    with open(os.path.join(site, DOCS, slug + ".md"), "w", encoding="utf-8") as f:
      f.write("\n".join(out))
    pages.append((slug, cat, len(opts)))

  missing_desc = sum(1 for o in options if not o["comment"])
  missing_tmpl = sum(1 for o in options if o["name"] not in examples)
  out = ["---", "title: Option Reference", "permalink: /{}/{}/".format(URL, REF_PREFIX), "---", "",
         "{% raw %}", "",
         "This reference lists every option of the SU2 configuration file. It is made from the source code, "
         "so it always matches the code. {}".format(stamp), "",
         "- Options: {}".format(len(options)),
         "- Options without a description in the code: {}".format(missing_desc),
         "- Options not in `config_template.cfg`: {}".format(missing_tmpl), "",
         "## Groups", ""]
  out += ["- [{}](../{}/) ({} options)".format(cat, slug, n) for slug, cat, n in pages]
  out += ["", "## All options A-Z", ""]
  for opt in sorted(options, key=lambda o: o["name"]):
    out.append("[`{0}`](../{1}/#{0}) ".format(opt["name"], opt["page"]))
  out += ["", "{% endraw %}", ""]
  with open(os.path.join(site, DOCS, REF_PREFIX + ".md"), "w", encoding="utf-8") as f:
    f.write("\n".join(out))
  return [REF_PREFIX] + [p[0] for p in pages]


def check_docs(site, known):
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
    hits, block = set(), None
    for num, line in enumerate(lines, 1):
      if line.lstrip().startswith("```"):
        if block is None:
          block = []
        else:
          names = [(n, m.group(1)) for n, m in block if m]
          if names and 2 * sum(name in known for _, name in names) >= len(names):
            hits.update((n, name) for n, name in names)
          block = None
        continue
      if block is not None:
        block.append((num, setting.match(line)))
      else:
        hits.update((num, name) for name in re.findall(r"`\s*([A-Z][A-Z0-9_]{2,})\s*=", line))
    for num, name in sorted(hits):
      if name not in known and name not in NOT_OPTIONS:
        found.append((os.path.relpath(path, site), num, name, lines[num - 1].strip()))
  return found


def write_check_page(site, found, commit):
  out = ["---", "title: Docs Check", "permalink: /{}/Docs-Check/".format(URL), "---", "", "{% raw %}", "",
         "These pages name options that do not exist in SU2 commit `{}`. "
         "The option was renamed or removed, or it has a typo.".format(commit), "",
         "Found: {}".format(len(found)), ""]
  last = None
  for path, num, name, line in found:
    if path != last:
      out += ["", "### `{}`".format(path), ""]
      last = path
    out.append("- line {}: `{}`".format(num, name))
  out += ["", "{% endraw %}", ""]
  with open(os.path.join(site, DOCS, "Docs-Check.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(out))


def patch_nav(site, ref_pages, extra_dev_pages):
  path = os.path.join(site, NAV)
  with open(path, encoding="utf-8") as f:
    nav = f.read()
  section = "- title: Option Reference\n  docs_v7:\n" + "".join("  - {}\n".format(p) for p in ref_pages) + "\n"
  nav = nav.replace("- title: Developer Docs", section + "- title: Developer Docs", 1)
  dev = re.search(r"- title: Developer Docs\n  docs_v7:\n((?:  - .*\n)+)", nav)
  if dev:
    listed = dev.group(1)
    added = "".join("  - {}\n".format(p) for p in extra_dev_pages if "  - {}\n".format(p) not in listed)
    nav = nav[:dev.end(1)] + added + nav[dev.end(1):]
  with open(path, "w", encoding="utf-8") as f:
    f.write(nav)


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--su2", required=True, help="root of the SU2 source tree")
  parser.add_argument("--site", required=True, help="checkout of su2code.github.io, changed in place")
  parser.add_argument("--light", action="store_true", help="remove the download folders (test site)")
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

  maps = read_maps(args.su2)
  options = read_options(args.su2)
  ref_pages = write_reference(args.su2, args.site, options, maps, commit)
  found = check_docs(args.site, {o["name"] for o in options})
  write_check_page(args.site, found, commit)
  patch_nav(args.site, ref_pages, dev_pages + ["Docs-Check"])

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
