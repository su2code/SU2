#!/usr/bin/env python3
"""Build and check the combined SU2 website in an isolated output directory."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

TOOLS = Path(__file__).resolve().parent


def output(command, **kwargs):
    return subprocess.check_output(command, text=True, **kwargs).strip()


def revision(path):
    return {"commit": output(["git", "-C", str(path), "rev-parse", "HEAD"]),
            "dirty": bool(output(["git", "-C", str(path), "status", "--porcelain"]))}


def prepare_directory(destination, sources):
    destination = destination.resolve()
    for source in sources:
        source = source.resolve()
        if destination == source or destination in source.parents or source in destination.parents:
            raise ValueError("Build directory must be separate from both source trees")
    marker = destination / ".su2-docs-build"
    if destination.exists() and any(destination.iterdir()) and not marker.is_file():
        raise ValueError("Refusing to reuse a nonempty directory without .su2-docs-build")
    destination.mkdir(parents=True, exist_ok=True)
    marker.touch()
    for name in ("source", "site", "doxygen"):
        path = destination / name
        if path.is_symlink():
            raise ValueError("Build output must not be a symlink: " + str(path))
        if path.exists():
            shutil.rmtree(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--su2", type=Path, default=TOOLS.parents[2])
    parser.add_argument("--site", type=Path, required=True, help="website source checkout (not modified)")
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--baseurl", default="/SU2")
    parser.add_argument("--mark-changes", action="store_true", help="keep red preview highlighting")
    args = parser.parse_args()
    su2, website, work = args.su2.resolve(), args.site.resolve(), args.build_dir.resolve()
    if not (su2 / "Common/src/CConfig.cpp").is_file() or not (website / "_config.yml").is_file():
        parser.error("--su2 and --site must point to their respective source checkouts")
    expected = (website / ".ruby-version").read_text().strip()
    actual = output(["ruby", "-e", "print RUBY_VERSION"])
    if actual != expected:
        parser.error("Use Ruby {} from the website .ruby-version (found {})".format(expected, actual))
    for executable in ("bundle", "doxygen", "dot"):
        if not shutil.which(executable):
            parser.error("Required tool not found: " + executable)
    sources = {"su2": revision(su2), "website": revision(website)}
    epoch = max(int(output(["git", "-C", str(p), "show", "-s", "--format=%ct", "HEAD"])) for p in (su2, website))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", SOURCE_DATE_EPOCH=str(epoch), BUNDLE_FROZEN="true")
    subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", str(TOOLS)], env=env, check=True)
    prepare_directory(work, (su2, website))
    source, site = work / "source", work / "site"
    shutil.copytree(website, source, ignore=shutil.ignore_patterns(
        ".git", ".github", "_site", ".jekyll-cache", ".jekyll-metadata", ".bundle", "vendor",
        "documents", "download", "emails", "__pycache__"))
    env["BUNDLE_GEMFILE"] = str(source / "Gemfile")
    subprocess.run(["bundle", "check"], cwd=source, env=env, check=True)
    flags = {"OUTPUT_DIRECTORY": str(work / "doxygen"), "PROJECT_NUMBER": sources["su2"]["commit"][:10],
             "SOURCE_BROWSER": "NO", "CALL_GRAPH": "NO", "CALLER_GRAPH": "NO", "COLLABORATION_GRAPH": "NO",
             "INCLUDE_GRAPH": "NO", "INCLUDED_BY_GRAPH": "NO", "DIRECTORY_GRAPH": "NO", "DOT_IMAGE_FORMAT": "svg",
             "WARN_LOGFILE": str(work / "doxygen_warnings.txt"), "EXCLUDE": "./Docs/website"}
    config = (su2 / "doxyfile").read_text() + "\n" + "\n".join(k + ' = "' + v + '"' for k, v in flags.items())
    print("Building Doxygen; log: " + str(work / "doxygen.log"), flush=True)
    with (work / "doxygen.log").open("w") as log:
        subprocess.run(["doxygen", "-"], cwd=su2, input=config + "\n", text=True, stdout=log,
                       stderr=subprocess.STDOUT, env=env, check=True)
    command = [sys.executable, str(TOOLS / "build_site.py"), "--su2", str(su2), "--site", str(source),
               "--light", "--doxygen-log", str(work / "doxygen_warnings.txt")]
    if args.mark_changes:
        command += ["--mark-changes", "--compare-site", str(website)]
    subprocess.run(command, env=env, check=True)
    subprocess.run(["bundle", "exec", "jekyll", "build", "--baseurl", args.baseurl,
                    "--destination", str(site)], cwd=source, env=env, check=True)
    subprocess.run([sys.executable, str(TOOLS / "rebase_links.py"), str(site), args.baseurl], env=env, check=True)
    shutil.copytree(work / "doxygen/html", site / "doxygen")
    sources.update({"baseurl": args.baseurl, "mark_changes": args.mark_changes,
                    "source_date_epoch": epoch, "light": True,
                    "gemfile_lock_sha256": hashlib.sha256((source / "Gemfile.lock").read_bytes()).hexdigest(),
                    "versions": {"python": sys.version.split()[0], "ruby": actual,
                                 "bundler": output(["bundle", "--version"], cwd=source, env=env),
                                 "jekyll": output(["bundle", "exec", "jekyll", "--version"], cwd=source, env=env),
                                 "doxygen": output(["doxygen", "--version"]),
                                 "graphviz": subprocess.run(["dot", "-V"], capture_output=True, text=True, check=True).stderr.strip()}})
    (site / "build-info.json").write_text(json.dumps(sources, indent=2) + "\n")
    subprocess.run([sys.executable, str(TOOLS / "check_site.py"), "--su2", str(su2), "--source", str(source),
                    "--site", str(site), "--baseurl", args.baseurl], env=env, check=True)
    print("Checked website: " + str(site))


if __name__ == "__main__":
    main()
