import ast
import os
import re
import sys

target_key = ".test_vals"
skip_key = ".test_vals_"

target = sys.argv[1]

with open(target) as f:
    scriptlines = f.readlines()

with open(sys.argv[2]) as f:
    loglines = f.readlines()


def update_filediff_references(target, scriptlines, loglines):
    """Rewrite on-disk reference files for failed run_filediff() tests.

    run_filediff() (TestCase.py) has two comparison modes:

    - Strict (tol_file_percent == 0.0): prints a unified diff (fromfile/tofile
      headers starting with "--- "/"+++ ").
    - Tolerance-based (tol_file_percent > 0.0): compares word-by-word and, on
      failure, prints "ERROR: File entries '<a>' and '<b>' in line N, word M
      differ" - the exact same message whether the mismatch is a genuine
      type/string difference or just a numeric delta over the tolerance, so
      the two must be told apart by checking if both entries parse as floats.
      "ERROR: Number of lines/words ... differ" indicates a structural
      difference and is never auto-applied.

    Either way, on failure the full new content of the test file is printed
    as a Python list repr, which is the exact new file content and can be
    written back as the new gold reference (unlike a diff hunk, which may
    only cover part of the file).
    """
    # Map test tag -> cfg_dir/reference_file(_aarch64), by scanning the
    # TestCase(...) blocks in the script.
    test_info = {}
    current_tag = None
    for line in scriptlines:
        m = re.search(r"TestCase\(\s*['\"]([^'\"]+)['\"]\s*\)", line)
        if m:
            current_tag = m.group(1)
            test_info[current_tag] = {}
            continue
        if current_tag is None:
            continue
        m = re.search(r"\.cfg_dir\s*=\s*['\"]([^'\"]+)['\"]", line)
        if m:
            test_info[current_tag]['cfg_dir'] = m.group(1)
            continue
        m = re.search(r"\.reference_file_aarch64\s*=\s*['\"]([^'\"]*)['\"]", line)
        if m:
            test_info[current_tag]['reference_file_aarch64'] = m.group(1)
            continue
        m = re.search(r"\.reference_file\s*=\s*['\"]([^'\"]*)['\"]", line)
        if m:
            test_info[current_tag]['reference_file'] = m.group(1)

    # Log lines may be prefixed by CI timestamps (e.g.
    # "2026-07-04T01:24:49.2142723Z "), so markers are searched for anywhere
    # in the line rather than anchored to its start.
    i = 0
    while i < len(loglines):
        m = re.search(r"==================== Start Test: (\S+) ====================", loglines[i])
        if not m:
            i += 1
            continue
        tag = m.group(1)
        i += 1
        block = []
        while i < len(loglines) and f"==================== End Test: {tag}" not in loglines[i]:
            block.append(loglines[i])
            i += 1

        if not any("Diff duration:" in line for line in block):
            continue  # Not a run_filediff() test.

        fromfile = None
        for line in block:
            m_from = re.search(r"--- (\S.*?)\t", line)
            if m_from:
                fromfile = m_from.group(1).strip()
                break

        if fromfile is None:
            # No unified-diff header: either the test passed, or it failed
            # via the tolerance-based word-by-word comparison.
            block_text = "".join(block)
            m_err = re.search(
                r"ERROR: File entries '([^']*)' and '([^']*)' in line \d+, word \d+ differ",
                block_text)
            is_pure_numeric_drift = False
            if m_err:
                try:
                    float(m_err.group(1))
                    float(m_err.group(2))
                    is_pure_numeric_drift = True
                except ValueError:
                    is_pure_numeric_drift = False

            if not is_pure_numeric_drift:
                if "ERROR" in block_text:
                    print(f"WARNING: {tag} failed with a non-numeric filediff error "
                          f"(missing reference file, mismatched line/word count, or a "
                          f"type/string mismatch); needs manual investigation, not just "
                          f"a value refresh.")
                continue

            info = test_info.get(tag, {})
            cpu_arch = ""
            for line in block:
                m_arch = re.search(r"CPU architecture:\s*(\S+)", line)
                if m_arch:
                    cpu_arch = m_arch.group(1)
                    break
            if "aarch64" in cpu_arch and info.get('reference_file_aarch64'):
                fromfile = info['reference_file_aarch64']
            else:
                fromfile = info.get('reference_file')

        if not fromfile:
            print(f"WARNING: could not determine reference file for test {tag} in "
                  f"{target}; skipping reference file update.")
            continue

        if tag not in test_info or 'cfg_dir' not in test_info[tag]:
            print(f"WARNING: could not find cfg_dir for test {tag} in {target}; "
                  f"skipping reference file update.")
            continue

        literal_line = None
        for line in reversed(block):
            stripped = line.strip()
            start, end = stripped.find("["), stripped.rfind("]")
            if start != -1 and end > start:
                literal_line = stripped[start:end + 1]
                break

        if literal_line is None:
            print(f"WARNING: filediff failure for {tag} detected but no new file "
                  f"contents found in log; skipping.")
            continue

        try:
            new_lines = ast.literal_eval(literal_line)
        except (ValueError, SyntaxError):
            print(f"WARNING: could not parse new file contents for {tag}; skipping.")
            continue

        ref_path = os.path.join(os.path.dirname(target), test_info[tag]['cfg_dir'], fromfile)
        with open(ref_path, "w") as f:
            f.writelines(new_lines)
        print(f"Updated reference file for {tag}: {ref_path}")


update_filediff_references(target, scriptlines, loglines)

# Jump to place in script where test values begin.
for cursor, line in enumerate(scriptlines):
    if target_key in line:
        break

for i, line in enumerate(loglines):
    if not "test_vals (stored): " in line: continue

    stored = line.split("test_vals (stored): ")[1].strip()
    computed = loglines[i+1].split("sim_vals (computed): ")[1].strip()

    if stored != computed:
        if not computed:
            # Empty computed values mean the test crashed/errored rather than
            # drifted numerically (e.g. "iteration could not be found").
            # Overwriting stored test_vals with an empty list would silently
            # destroy the reference data, so skip it and warn instead.
            print(f"WARNING: empty computed values for {scriptlines[cursor].split('=')[0].strip()} "
                  f"in {target} (log line {i+1}); leaving stored test_vals unchanged.")
        else:
            newline = scriptlines[cursor].split("=")[0]+"= ["+computed+"]\n"
            scriptlines[cursor] = newline

    for offset, line in enumerate(scriptlines[cursor+1:]):
        if target_key in line and not skip_key in line:
            if not "#" in line.split(target_key)[0]:
                break
    cursor += offset+1

with open(target,"w") as f:
    f.writelines(scriptlines)
