#!/usr/bin/env python

## \file resolver.py
#  \brief Utility for custom metric sensors in mesh adaptation.
#  \author B. Munguía
#  \version 8.4.0 "Harrier"
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


class CustomSensorRegistry:
    """Registry mapping sensor names to per-node callables for mesh adaptation.

    Usage::

        registry = CustomSensorRegistry()
        registry["MACH"] = my_sensor  # fn(driver, iPoint) -> float

        driver = pysu2.CSinglezoneDriver("case.cfg", 1, comm)
        registry.initialize(driver)   # resolve indices once after construction

        for inner_iter in range(N):
            driver.Preprocess(inner_iter)
            driver.Run()
            registry.populate(driver) # fill sensor slots before ComputeMetricField
            driver.Postprocess()
            driver.Update()
            driver.Output(inner_iter)
    """

    def __init__(self):
        self._sensors = {}  # name -> fn(driver, iPoint) -> float
        self._indices = {}  # name -> iSensor (cached by initialize)

    def register(self, name, fn):
        """Register a sensor function ``fn(driver, iPoint) -> float``.

        Args:
            name: Sensor name exactly as listed in METRIC_SENSOR in the config.
            fn:   Called once per node per iteration. Must cover all nodes
                  (domain + halo) so the Green-Gauss gradient has correct
                  values at MPI boundaries.
        """
        self._sensors[name] = fn

    def __setitem__(self, name, fn):
        self.register(name, fn)

    def __getitem__(self, name):
        return self._sensors[name]

    def initialize(self, driver):
        """Resolve and cache sensor indices. Call once after driver construction.

        Raises:
            ValueError: if a registered name is not listed in METRIC_SENSOR.
        """
        for name in self._sensors:
            idx = driver.GetMetricSensorIndex(name)
            if idx < 0:
                raise ValueError(
                    "Custom sensor '{}' not found in the driver. "
                    "Ensure it is listed in METRIC_SENSOR in the config.".format(name)
                )
            self._indices[name] = idx

    def populate(self, driver):
        """Evaluate all sensors and push values to the driver. Call before Postprocess()."""
        nNodes = driver.GetNumberNodes()
        for name, fn in self._sensors.items():
            iSensor = self._indices[name]
            for iPoint in range(nNodes):
                driver.SetSensorAdapt(iPoint, iSensor, fn(driver, iPoint))
