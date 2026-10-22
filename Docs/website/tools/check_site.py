#!/usr/bin/env python3
"""Check generated SU2 documentation without requesting external websites."""

import argparse
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import unquote, urljoin, urlsplit

import build_site

COLLECTIONS = ("docs_v7", "tutorials", "vandv")


class HTML(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.ids, self.links = set(), []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        self.ids.update(attrs[k] for k in ("id", "name") if k in attrs)
        self.links.extend(attrs[k].strip() for k in ("href", "src") if attrs.get(k))


def local_target(root, page, href, baseurl):
    base = "/" + baseurl.strip("/") if baseurl.strip("/") else ""
    link = urlsplit(urljoin("https://docs.invalid" + base + "/" + page, href))
    if link.scheme not in ("http", "https") or link.netloc != "docs.invalid":
        return None
    path = unquote(link.path)
    if base and path != base and not path.startswith(base + "/"):
        return root / "__outside_baseurl__", path, ""
    path = path[len(base):].lstrip("/")
    target = (root / path).resolve()
    if root != target and root not in target.parents:
        return root / "__outside_site__", path, ""
    if target.is_dir():
        target /= "index.html"
    return target, path, unquote(link.fragment)


def check_links(root, baseurl):
    root = Path(root).resolve()
    cache, problems, checked = {}, set(), 0

    def parse(path):
        if path not in cache:
            cache[path] = HTML(path.read_text(encoding="utf-8", errors="replace"))
        return cache[path]

    for collection in COLLECTIONS:
        for page in sorted((root / collection).rglob("*.html")):
            relative = page.relative_to(root).as_posix()
            for href in parse(page).links:
                result = local_target(root, relative, href, baseurl)
                if result is None:
                    continue
                target, path, fragment = result
                checked += 1
                reason = ""
                if not target.is_file():
                    reason = "missing file"
                elif fragment and target.suffix == ".html" and fragment not in parse(target).ids:
                    reason = "missing anchor"
                if reason:
                    problems.add("{}: {}{} ({})".format(
                        relative, path, "#" + fragment if fragment else "", reason))
    return problems, checked



def check_enum_values(source, options, maps):
    choices = {o["name"]: set(build_site.describe_type(o, maps)[1])
               for o in options if o["kind"] in ("Enum", "EnumList")}
    known = {o["name"] for o in options}
    problems = set()
    for collection in COLLECTIONS:
        for page in (source / ("_" + collection)).rglob("*.md"):
            text = page.read_text(encoding="utf-8")
            for block in re.finditer(r"^```([^\n]*)\n(.*?)^```", text, re.M | re.S):
                if block[1].strip().lower() in ("cpp", "c++", "c", "hpp", "python", "py", "bash", "sh", "shell"):
                    continue
                settings = re.findall(r"^\s*([A-Z][A-Z0-9_]+)\s*=\s*([^%\n]*)", block[2], re.M)
                if not settings or 2 * sum(name in known for name, _ in settings) < len(settings):
                    continue
                for name, value in settings:
                    value = value.strip()
                    if not choices.get(name) or not re.fullmatch(r"[A-Za-z0-9_, ()]+", value):
                        continue
                    keywords = re.findall(r"[A-Za-z0-9_]+", value)
                    for keyword in set(keywords) - choices[name]:
                        if keyword.startswith("__") and keyword.endswith("__"):
                            continue
                        problems.add("{}: {}= {} (unknown keyword {})".format(
                            page.relative_to(source).as_posix(), name, value, keyword))
    return problems


def check_source(su2, source, baseline):
    options = build_site.read_options(su2)
    known = {o["name"] for o in options}
    problems = ["{}:{}: unknown option {}".format(p, n, name)
                for p, n, name, _ in build_site.check_docs(str(source), known, build_site.read_code_names(su2))]
    template = {name for _, name in build_site.check_template(su2, known)}
    problems += ["config_template.cfg: unknown option " + name
                 for name in sorted(template - set(baseline["template_options"]))]
    missing = {o["name"] for o in options
               if not build_site.clean_comment(o["comment"] or "")[0]}
    problems += ["New option without a source description: " + name
                 for name in sorted(missing - set(baseline["missing_descriptions"]))]
    problems += sorted(check_enum_values(source, options, build_site.read_maps(su2)) - set(baseline["enum_values"]))
    return problems


def check_search(root, source, su2, baseurl):
    data = json.loads((root / "search.json").read_text())
    problems, urls = [], set()
    for entry in data:
        if not all(isinstance(entry.get(k), str) and entry[k] for k in ("title", "url", "text", "section")):
            problems.append("Search entry has missing or empty fields: " + str(entry.get("url")))
            continue
        url = entry["url"]
        if url in urls:
            problems.append("Duplicate search URL: " + url)
        urls.add(url)
        target = local_target(root, "index.html", url, baseurl)
        if target is None or not target[0].is_file():
            problems.append("Search URL has no local page: " + url)
    base = "/" + baseurl.strip("/") if baseurl.strip("/") else ""
    for collection in COLLECTIONS:
        for page in (source / ("_" + collection)).rglob("*.md"):
            text = page.read_text(encoding="utf-8")
            front = re.match(r"\A---\s*\n(.*?)\n---", text, re.S)
            if not front or not re.search(r"^title:", front[1], re.M):
                continue
            permalink = re.search(r"^permalink:\s*(\S+)", front[1], re.M)
            if permalink and base + "/" + permalink[1].strip("\"'/") + "/" not in urls:
                problems.append("Page missing from search: " + str(page.relative_to(source)))
    indexed = " ".join(entry.get("text", "") for entry in data).lower()
    fields = {f["name"] for fs in build_site.read_outputs(su2).values() for f in fs if f["literal"]}
    problems += ["Output field missing from search: " + name for name in sorted(fields) if name.lower() not in indexed]
    return problems, len(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--su2", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path, help="staged website after generation")
    parser.add_argument("--site", required=True, type=Path, help="assembled HTML website")
    parser.add_argument("--baseurl", default="/SU2")
    args = parser.parse_args()
    baseline = json.loads(Path(__file__).with_name("check-baseline.json").read_text())
    problems = check_source(args.su2, args.source, baseline)
    links, count = check_links(args.site, args.baseurl)
    problems += sorted(links - set(baseline["links"]))
    search, entries = check_search(args.site.resolve(), args.source, args.su2, args.baseurl)
    problems += search
    report = {"new_problems": problems, "references_checked": count, "search_entries": entries,
              "existing_link_problems": sorted(links & set(baseline["links"]))}
    (args.site / "docs-check.json").write_text(json.dumps(report, indent=2) + "\n")
    for problem in problems:
        print(problem)
    print("Checked {} local references and {} search entries; {} new problems, {} existing link problems.".format(
        count, entries, len(problems), len(report["existing_link_problems"])))
    return bool(problems)


if __name__ == "__main__":
    raise SystemExit(main())
