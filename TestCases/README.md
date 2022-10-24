--------------------------------
  SU2: The Open-Source CFD Code
--------------------------------

Computational analysis tools have revolutionized the way we design aerospace systems, but most established codes are proprietary, unavailable, or prohibitively expensive for many users. The SU2 team is changing this, making computational analysis and design freely available as open-source software and involving everyone in its creation and development.

-----------------------------------------------------------
  SU2 TEST CASES
-----------------------------------------------------------

SU2 provides an extensive collection of common test cases. The test cases contain meshes and configuration files that can be run to ensure different components of the suite are working properly. The directory structure is organized by governing equations and their specific cases.

The files found within these directories serve two purposes:

1. A subset of the available cases are used for regression testing internally by the development team. The configuration files and meshes are automatically downloaded and executed with the latest version of SU2 at regular intervals. Any changes in the output of the specified test cases are immediately reported to the developers. These cases are controlled by the Python scripts in the root of the test cases directory, e.g., serial_regression.py.

2. The entire suite of test cases is provided to the community as a way to get started with SU2 and its many configuration options, including settings that the developers consider to be good starting points. Often, you will find a test case that is similar to your problem of interest, and the available configuration files can be modified to suit your needs.

Note that, while many of the cases are used for regression testing, this test case suite is provided **without any guarantees on performance or expected results** (check out the [tutorials](https://su2code.github.io/tutorials/home/) for more thoroughly maintained cases). Keep in mind that small changes to the configuration options can often result in large changes to the output. We encourage the community to experiment with these test cases, or even try (and share!) other interesting validation cases with SU2! 

Happy testing!

- The SU2 Development Team
