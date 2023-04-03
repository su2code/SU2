/*!
 * \file CTurbomachineryOutput.hpp
 * \brief  Headers of the turbomachinery output.
 * \author F. Palacios, T. Economon, M. Colonno, J. Kelly
 * \version 7.5.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../../Common/include/option_structure.hpp"

class CGeometry;
class CConfig;
class CSolver;
class CIntegration;

/*!
 * \class CTurbomachineryOutput
 * \brief Class for writing the turbomachinery output
 * \author J. Kelly
 */
class CTurbomachineryOutput{

    unsigned short nSpanWiseSections, nMarkerTurboPerf;

    su2double **TotalStaticEfficiency,
        **TotalTotalEfficiency,
        **KineticEnergyLoss,
        **TRadius,
        **TotalPressureLoss,
        **MassFlowIn,
        **MassFlowOut,
        **FlowAngleIn,
        **FlowAngleIn_BC,
        **FlowAngleOut,
        **EulerianWork,
        **TotalEnthalpyIn,
        **TotalEnthalpyIn_BC,
        **EntropyIn,
        **EntropyOut,
        **EntropyIn_BC,
        **PressureRatio,
        **TotalTemperatureIn,
        **EnthalpyOut,
        ***MachIn,
        ***MachOut,
        **VelocityOutIs,
        **DensityIn,
        **PressureIn,
        ***TurboVelocityIn,
        **DensityOut,
        **PressureOut,
        ***TurboVelocityOut,
        **EnthalpyOutIs,
        **EntropyGen,
        **AbsFlowAngleIn,
        **TotalEnthalpyOut,
        **RothalpyIn,
        **RothalpyOut,
        **TotalEnthalpyOutIs,
        **AbsFlowAngleOut,
        **PressureOut_BC,
        **TemperatureIn,
        **TemperatureOut,
        **TotalPressureIn,
        **TotalPressureOut,
        **TotalTemperatureOut,
        **EnthalpyIn,
        **TurbIntensityIn,
        **Turb2LamViscRatioIn,
        **TurbIntensityOut,
        **Turb2LamViscRatioOut,
        **NuFactorIn,
        **NuFactorOut;

    protected:
        int rank,
        size;
    public:
    /*!
   * \brief Constructor of the class.
   */
  CTurbomachineryOutput(CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbomachineryOutput(void);

};