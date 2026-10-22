#!/usr/bin/env python3

## \file rebase_links.py
#  \brief Make the absolute links of a built SU2 website work under a sub-path
#         (e.g. https://<user>.github.io/SU2/). Links to the download folders that
#         the test site does not carry go to the real website.
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

import os
import re
import sys

REAL_SITE = "https://su2code.github.io"
REMOTE_FOLDERS = ("documents/", "download/", "emails/")


def main():
  if len(sys.argv) != 3:
    print("usage: rebase_links.py <built site folder> <base path, e.g. /SU2>")
    return 1
  root = sys.argv[1]
  base = "/" + sys.argv[2].strip("/") if sys.argv[2].strip("/") else ""
  link = re.compile(r'((?:href|src|action)\s*=\s*["\'])/(?!/)')

  def fix(m, text):
    rest = text[m.end():]
    if base and rest.startswith(base.lstrip("/") + "/"):
      return m.group(0)
    if rest.startswith(REMOTE_FOLDERS):
      return m.group(1) + REAL_SITE + "/"
    return m.group(1) + base + "/"

  changed = 0
  for folder, _, files in os.walk(root):
    for name in files:
      if not name.endswith((".html", ".xml")):
        continue
      path = os.path.join(folder, name)
      with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
      new = link.sub(lambda m: fix(m, text), text)
      if new != text:
        with open(path, "w", encoding="utf-8") as f:
          f.write(new)
        changed += 1
  print("Rebased links in {} files".format(changed))
  return 0


if __name__ == "__main__":
  sys.exit(main())
