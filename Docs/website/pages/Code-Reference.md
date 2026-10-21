---
title: Code Reference (Doxygen)
permalink: /docs_v7/Code-Reference/
---

The code reference is made by [Doxygen](https://www.doxygen.nl) from the comments in the SU2 source code.
It is rebuilt each time the code changes, so it always matches the code.

**[Open the code reference]({{ "/doxygen/index.html" | relative_url }})**

Good places to start:

- [Class list]({{ "/doxygen/annotated.html" | relative_url }}): all classes, with their parents and children.
- [File list]({{ "/doxygen/files.html" | relative_url }}): all source files.
- `CDriver`: runs the simulation and creates all the other objects.
- `CConfig`: reads the configuration file. The [Option Reference](../Option-Reference/) is made from it.
- `CSolver`, `CNumerics`, `CVariable`, `CGeometry`: the solvers, the numerical schemes, the stored variables and the mesh.

To build it on your computer, run this in the SU2 folder. The result is in `Docs/html/index.html`.

```
doxygen doxyfile
```
