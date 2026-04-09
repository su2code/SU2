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

from __future__ import annotations

from collections.abc import Callable
from typing import Any

SU2Driver = Any
SensorFn = Callable[[SU2Driver, int], float]


class CustomSensorRegistry:
    """Registry mapping sensor names to callables for mesh adaptation."""

    def __init__(self, sensors: dict[str, SensorFn] | None = None) -> None:
        """
        Args:
            sensors: Optional mapping of sensor name to callable. Each callable
                must have the signature ``fn(driver, iPoint) -> float`` and names
                must match entries in ``METRIC_SENSOR`` in the config exactly.
        """
        self._sensors: dict[str, SensorFn] = (
            dict(sensors) if sensors is not None else {}
        )
        self._indices: dict[str, int] = {}

    def __setitem__(self, name: str, fn: SensorFn) -> None:
        """Register or replace the sensor callable for ``name``."""
        self._sensors[name] = fn

    def __getitem__(self, name: str) -> SensorFn:
        """Return the sensor callable registered under ``name``."""
        return self._sensors[name]

    def initialize(self, driver: SU2Driver) -> None:
        """Resolve and cache sensor indices. Call once after driver construction.

        Args:
            driver: Active SU2 driver instance (e.g. ``CSinglezoneDriver``).

        Raises:
            ValueError: If a registered sensor name is not listed in
                ``METRIC_SENSOR`` in the config.
        """
        for name in self._sensors:
            idx = driver.GetMetricSensorIndex(name)
            if idx < 0:
                raise ValueError(
                    "Custom sensor '{}' not found in the driver. "
                    "Ensure it is listed in METRIC_SENSOR in the config.".format(name)
                )
            self._indices[name] = idx

    def populate(self, driver: SU2Driver) -> None:
        """Evaluate all sensors and push values to the driver.

        Call after ``Run()`` and before ``Postprocess()`` on each iteration.

        Args:
            driver: Active SU2 driver instance (e.g. ``CSinglezoneDriver``).
        """
        nNodes = driver.GetNumberNodes()
        for name, fn in self._sensors.items():
            iSensor = self._indices[name]
            for iPoint in range(nNodes):
                driver.SetSensorAdapt(iPoint, iSensor, fn(driver, iPoint))
