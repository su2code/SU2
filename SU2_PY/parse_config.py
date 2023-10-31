#!/usr/bin/env python

## \file parse_config.py
#  \brief Builds a worksheet of all SU2.cpp options
#  \author A. Aranake, F. Palacios
#  \version 8.0.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

# make print(*args) function available in PY2.6+, does'nt work on PY < 2.6
from __future__ import print_function

import os, sys, xlwt

# note: requires xlwt for spreadsheet output
# http://pypi.python.org/pypi/xlwt


class config_option:
    option_name = ""
    option_type = ""
    option_category = ""
    option_values = []
    option_default = ""
    option_description = ""

    def __init__(self, name, otype, category, values, default, description):
        self.option_name = name
        self.option_type = otype
        self.option_category = category
        self.option_values = values
        self.option_default = default
        self.option_description = description

    def print_data(self):
        print("Option Name: %s " % self.option_name)
        print("Option Type: %s " % self.option_type)
        print("Option Category: %s " % self.option_category)
        print("Option values: ", self.option_values)
        print("Option default: %s" % self.option_default)
        print("Option description: %s " % self.option_description)
        print("")


def parse_config(config_cpp, config_hpp):

    # List of option types
    option_types = [
        "AddEnumOption",
        "AddMathProblem",
        "AddSpecialOption",
        "AddScalarOption",
        "AddMarkerOption",
        "AddMarkerPeriodic",
        "AddMarkerDirichlet",
        "AddMarkerInlet",
        "AddMarkerOutlet",
        "AddMarkerDisplacement",
        "AddMarkerLoad",
        "AddMarkerFlowLoad",
        "AddArrayOption",
        "AddListOption",
        "AddConvectOption",
        "AddEnumListOption",
        "AddDVParamOption",
    ]

    # Build a dictionary of enum options from hpp file
    enum_options = {}
    f = open(config_hpp, "r")
    while 1:  # Find beginning of enum definitions
        s = f.readline()
        if s.find("BEGIN_CONFIG_ENUMS") > -1:
            break
    while 1:
        s = f.readline()
        if s.find("END_CONFIG_ENUMS") > -1:
            break  # Reached end

        if s.find("CCreateMap") > -1:
            dict_key = (s.split("=")[0]).split(">")[1].strip()
            dict_val = []
            while 1:
                s2 = f.readline()
                thisval = s2.split('"')[1]
                dict_val.append(thisval)
                if s2.find(";") > -1:
                    break
            enum_options[dict_key] = dict_val
    f.close()

    # Temporary: For now, build a list of all of the schemes
    scheme_list = enum_options["Upwind_Map"]
    scheme_list.extend(enum_options["Centered_Map"][1:])
    # Read the Options section of config_structure.cpp into a list of strings
    lines = []
    f = open(config_cpp, "r")
    while 1:
        s = f.readline()
        if s.find("BEGIN_CONFIG_OPTIONS") > -1:
            break
    while 1:
        s = f.readline()
        # Check if we've reached the end
        if s.find("END_CONFIG_OPTIONS") > -1:
            break
        lines.append(s)
    f.close()

    option_list = []
    present_category = "None"
    # ----- Main text parsing loop -----
    for j, line in enumerate(lines):

        # Check for a category description
        if line.find("CONFIG_CATEGORY") > -1:
            present_category = line.split(":")[1].strip().strip("*/").strip()
            print(present_category)

        # Check for an option type
        for option_type in option_types:
            if line.find(option_type) > -1:  # Found an option
                # Get option name
                name = line.split('"')[1]

                # Permitted values
                values = ["YES", "NO"]
                if option_type == "AddEnumOption":
                    try:
                        enum_mapname = line.split(",")[2].strip()
                        values = enum_options[enum_mapname]
                    except KeyError:
                        print("KeyError, key=%s" % enum_mapname)
                        print("enum_options: ", enum_options)
                        sys.exit(1)
                    except TypeError:
                        print("TypeError, key=%s" % enum_mapname)
                        print("enum_options: ", enum_options)
                        sys.exit(1)
                elif option_type == "AddMathProblem":
                    values = ["DIRECT", "CONTINUOUS_ADJOINT", "LINEARIZED"]
                elif option_type == "AddScalarOption":
                    values = ["A scalar constant"]
                elif option_type in (
                    "AddMarkerOption",
                    "AddMarkerPeriodic",
                    "AddMarkerDirichlet",
                    "AddMarkerInlet",
                    "AddMarkerOutlet",
                    "AddMarkerDisplacement",
                    "AddMarkerLoad",
                    "AddMarkerFlowLoad",
                ):
                    values = ["Valid marker name from grid file"]
                elif option_type == "AddArrayOption":
                    values = ["Array"]
                elif option_type == "AddListOption":
                    values = ["List"]
                elif option_type == "AddConvectOption":
                    values = scheme_list
                    print("Convect Option: ", name)
                elif option_type == "AddEnumListOption":
                    values = ["Enum list"]
                elif option_type == "AddDVParamOption":
                    values = ["DV Param"]

                # A first pass at finding the default value (Check the last item in parenthesis)
                jdefault = j
                while lines[jdefault].find(";") == -1:
                    jdefault = jdefault + 1
                default = (
                    lines[jdefault]
                    .strip()
                    .strip(");")
                    .split(",")[-1]
                    .strip()
                    .strip('"')
                )

                # A whole bunch of corrections for what the default should be...
                if default.find('string("') == 0:
                    default = default.split('"')[1]
                if default == "default_vec_3d":
                    default = "(1.0, 100.0, 1.0)"
                if default == "default_vec_6d":
                    default = "( -1E15, -1E15, -1E15, 1E15, 1E15, 1E15 )"
                if option_type == "AddMathProblem":
                    default = "DIRECT"
                if option_type == "AddConvectOption":
                    default = "ROE-1ST_ORDER"
                if default == "RK_Alpha_Step":
                    default = "( 0.66667, 0.66667, 1.000000 )"
                if default == "RK_Beta_Step":
                    default = "( 1.00000, 0.00000, 0.00000 )"

                if default == "false":
                    default = "NO"
                elif default == "true":
                    default = "YES"

                # Check for a description tag
                description = "No description"
                if lines[j - 1].find("DESCRIPTION") > -1:
                    description = lines[j - 1].split(":")[1].strip().strip("*/").strip()

                # Add a new option
                option_list.append(
                    config_option(
                        name,
                        option_type[3:],
                        present_category,
                        values,
                        default,
                        description,
                    )
                )

                break

    return option_list


def print_all(option_list):
    # Dumps the option list to screen
    for option in option_list:
        option.print_data()


def make_spreadsheet(filename, option_list):

    wbk = xlwt.Workbook()

    sheet_name = ""
    jp = 0
    for j, opt in enumerate(option_list):
        jp = jp + 1
        if not sheet_name == opt.option_category:
            # Create new sheet for new category
            sheet_name = opt.option_category
            sheet = wbk.add_sheet(sheet_name[:31].replace("/", "-"))

            # Write spreadsheet header
            sheet.write(0, 0, "Option Name")
            sheet.write(0, 1, "Option Type")
            sheet.write(0, 2, "Option Category")
            sheet.write(0, 3, "Option Values")
            sheet.write(0, 4, "Option Default")
            sheet.write(0, 5, "Option Description")
            jp = 1

        sheet.write(jp, 0, opt.option_name)
        sheet.write(jp, 1, opt.option_type)
        sheet.write(jp, 2, opt.option_category)
        sheet.write(jp, 3, (",").join(opt.option_values))
        sheet.write(jp, 4, opt.option_default)
        sheet.write(jp, 5, opt.option_description)

    wbk.save(filename)


if __name__ == "__main__":

    # These variables should point to the configuration files
    su2_home = os.environ["SU2_HOME"]
    config_cpp = os.path.join(su2_home, "Common/src/config_structure.cpp")
    config_hpp = os.path.join(su2_home, "Common/include/option_structure.hpp")

    # Check that files exist
    if not os.path.isfile(config_cpp):
        sys.exit(
            "Could not find cpp file, please check that su2_basedir is set correctly in parse_config.py"
        )
    if not os.path.isfile(config_hpp):
        sys.exit(
            "Could not find hpp file, please check that su2_basedir is set correctly in parse_config.py"
        )

    # Run the parser
    option_list = parse_config(config_cpp, config_hpp)

    # Dump parsed data to screen
    print_all(option_list)

    # make_spreadsheet('out.xls',option_list)
