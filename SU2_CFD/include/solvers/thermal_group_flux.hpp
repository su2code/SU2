/*!
 * \file thermal_group_flux.hpp
 * \brief Thermal residuals and temperature derivatives at grouped boundary vertices.
 * \version 8.6.0 "Harrier"
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 * Distributed under the GNU Lesser General Public License, version 2.1 or later.
 */
#pragma once

#include <cstddef>

namespace ThermalGroupFlux {

template <class T, class Group>
bool Isothermal(bool eligible, std::size_t count, const T& conductivity, const T& wallTemperature,
                Group group, T& flux, T& derivative) {
  if (!eligible) return false;
  flux = 0.0;
  derivative = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    T area = 0.0, distance = 0.0, temperature = 0.0;
    group(i, area, distance, temperature);
    if (!(distance > 0.0)) return false;
    flux += conductivity * (wallTemperature - temperature) / distance * area;
    derivative += conductivity * area / distance;
  }
  return true;
}

template <class T>
void HeatTransfer(const T& coefficient, const T& farTemperature, const T& wallTemperature, const T& area,
                  T& flux, T& residualDerivative) {
  flux = coefficient * (farTemperature - wallTemperature) * area;
  residualDerivative = coefficient * area;
}

}  // namespace ThermalGroupFlux
