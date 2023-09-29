/*!
 * \file prop_defs.hpp
 * \brief Headers of the CEulerSolver class
 * \author Chandukrishna Y., T. N. Venkatesh and Josy P. Pullockara
 * \version 8.0.0 "Harrier"
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

#define BEM_MAXR      50
#define BEM_MAXALF   100
#define BEM_MAX_ITER  20

typedef struct
{
  int nblades;
  float dia;
  float rhub;
  float ang0_75;
} propeller_geom_struct;

typedef struct
{
  int nblades;
  su2double dia;
  su2double rhub;
  su2double ang0_75;
} dpropeller_geom_struct;

typedef struct
{
  int nalf;
  int nrad;
  float r1[BEM_MAXR];
  float chord[BEM_MAXR];
  float setangle[BEM_MAXR];
  float alf[BEM_MAXALF][BEM_MAXR];
  float cl_arr[BEM_MAXALF][BEM_MAXR];
  float cd_arr[BEM_MAXALF][BEM_MAXR];
} propeller_section_struct;

typedef struct
{
  int nalf;
  int nrad;
  su2double r1[BEM_MAXR];
  su2double chord[BEM_MAXR];
  su2double setangle[BEM_MAXR];
  su2double alf[BEM_MAXALF][BEM_MAXR];
  su2double cl_arr[BEM_MAXALF][BEM_MAXR];
  su2double cd_arr[BEM_MAXALF][BEM_MAXR];
} dpropeller_section_struct;
