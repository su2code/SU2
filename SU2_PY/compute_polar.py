#!/usr/bin/env python

## \file Compute_polar.py
#  \brief Python script for performing polar sweep.
#  \author E Arad (based on T. Lukaczyk and  F. Palacios script)
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
#
#
# Several combinations of angles are possible:
# ------------------------------------------------
# 1. Polar-sweep in alpha per given phi                        ...... polarVar   = aoa
# 2. Polar-sweep in alpha per given beta (side slip angle)     ...... polarVar   = aoa
# 3. Polar-sweep in phi per given alpha                        ...... polarVar   = phi
# 4. Mach ramp  (single values for alpha, phi or both permitted)  ... polarVar   = MachRampNumbers
#
# Note: Seting a list of both phi and beta is impossible
#       For mach ramp you can specify alpha, phi (or both), but not a list of either of them

# make print(*args) function available in PY2.6+, does'nt work on PY < 2.6
from __future__ import print_function

# imports
import os, sys, shutil
from optparse import OptionParser

sys.path.append(os.environ["SU2_RUN"])
import SU2
import SU2.util.polarSweepLib as psl
import copy
import numpy as np


def main():
    # Command Line Options
    parser = OptionParser()
    parser.add_option(
        "-c",
        "--ctrl",
        dest="ctrlFile",
        help="reads polar control parameters from FILE (default:polarCtrl.in) ",
        metavar="FILE",
        default="polarCtrl.in",
    )
    parser.add_option(
        "-n",
        "--partitions",
        dest="partitions",
        default=2,
        help="number of PARTITIONS",
        metavar="PARTITIONS",
    )
    parser.add_option(
        "-i",
        "--iterations",
        dest="iterations",
        default=-1,
        help="number of ITERATIONS",
        metavar="ITERATIONS",
    )
    parser.add_option(
        "-d",
        "--dimension",
        dest="geomDim",
        default=2,
        help="Geometry dimension (2 or 3)",
        metavar="geomDim",
    )
    parser.add_option(
        "-w",
        "--Wind",
        action="store_true",
        dest="Wind",
        default=False,
        help=" Wind system (default is body system)",
    )
    parser.add_option(
        "-v",
        "--Verbose",
        action="store_true",
        dest="verbose",
        default=False,
        help=" Verbose printout (if activated)",
    )

    (options, args) = parser.parse_args()
    options.partitions = int(options.partitions)
    options.iterations = int(options.iterations)
    options.geomDim = int(options.geomDim)

    d2r = np.pi / 180
    #
    sweepOption = []
    sweepOption.append(" Polar sweep type: 1. Sweep in AOA per given roll angle")
    sweepOption.append(" Polar sweep type: 2. Sweep in AOA per given sideslip-angle")
    sweepOption.append(" Polar sweep type: 3. Sweep in phi per given AOA")
    sweepOption.append(
        " Polar sweep type: 4. Mach ramp (single- value AOA and  sideslip-angle"
    )

    #
    # --------------- now read the parameters control file and parse it

    fc = open(options.ctrlFile, "r")
    ctrl = fc.readlines()
    nc = np.size(ctrl)
    fc.close()

    print(str(nc) + " lines read from control file: " + options.ctrlFile)

    (
        PA,
        polarSweepType,
        velDirOption,
        nAlpha,
        nBeta,
        nPhi,
        nMach,
        alpha,
        beta,
        phi,
        MachList,
        polarVar,
    ) = psl.setPolaraType(ctrl, nc, options.verbose)

    if options.verbose:
        velDirOptionLegend = ["V(alpha,phi)", "V(alpha,beta)"]
        print(
            ">>>  Control file details: Pitch axis is "
            + PA
            + ". Polar sweep type is "
            + str(polarSweepType)
            + "; polarVar = "
            + polarVar
        )
        print(">>>  Velocity definiton: " + velDirOptionLegend[velDirOption - 1])
        print(
            ">>>  nAalpha = "
            + str(nAlpha)
            + "; nBeta = "
            + str(nBeta)
            + "; nPhi = "
            + str(nPhi)
            + "; nMach = "
            + str(nMach)
        )
    if polarSweepType < 4:
        nPolara = max(nAlpha, nPhi)
    else:
        nPolara = nMach

    # -------------Configuration base file ----------------------
    inputbaseFileString = "input base file"
    keyWordInputbaseFile = inputbaseFileString.lower()
    iBaseInputF = psl.parLocator(keyWordInputbaseFile, ctrl, nc, -1, options.verbose)
    bIFLine = ctrl[iBaseInputF]
    icol = bIFLine.index(":")
    sBIF = bIFLine[icol + 1 :]
    inputbaseFile = sBIF.strip(" ")
    inputbaseFile = inputbaseFile.strip("\n")

    print(" ")
    print(
        "--------------------------------------------------------------------------------------"
    )
    print(" ")
    print("Configuration file: " + inputbaseFile)
    print(
        "PolarSweepType = "
        + str(polarSweepType)
        + " Polar sweep in "
        + polarVar
        + " using "
        + str(nPolara)
        + " angles/Mach No "
    )
    print(" ")
    print(
        "--------------------------------------------------------------------------------------"
    )
    print(" ")

    if polarSweepType == 4:
        nPolara = 1  # prevent angles inner loop
    if options.geomDim not in [2, 3]:
        raise SystemExit("ERROR: dimension can be either 2 or 3 (-d parameter)  ")

    if options.Wind:
        outSystem = "Wind"
    else:
        outSystem = "Body"

    print(" ")
    print(
        "==============================================================================="
    )
    print(
        "   Polar sweep in "
        + str(options.geomDim)
        + "D ; output in "
        + outSystem
        + " system"
    )
    print(
        "==============================================================================="
    )
    print(" ")

    # load config, start state
    config = SU2.io.Config(inputbaseFile)
    state = SU2.io.State()
    # Set SU2 defaults units, if definitions are not included in the cfg file
    if "SYSTEM_MEASUREMENTS" not in config:
        config.SYSTEM_MEASUREMENTS = "SI"
    if config.SOLVER == "NAVIER_STOKES":
        if "REYNOLDS_LENGTH" not in config:
            config.REYNOLDS_LENGTH = 1.0

    # prepare config
    config.NUMBER_PART = options.partitions
    if options.iterations > 0:
        config.ITER = options.iterations
    config.NZONES = 1

    # find solution files if they exist
    state.find_files(config)

    # start results data
    results = SU2.util.bunch()

    if nMach == 0:
        if "MACH_NUMBER" in config:
            MachList.append(config.MACH_NUMBER)
        else:
            MachList.append(0.5)
        nMach = 1

    if nAlpha == 0:
        if "AOA" in config:
            alpha.append(config.AOA)
        else:
            alpha.append(0.0)
        nAlpha = 1

    if nPhi == 0:
        phi.append(0.0)
        nPhi = 1
        noPhi_in_CTRL = True
    else:
        noPhi_in_CTRL = False

    if nBeta == 0:
        if noPhi_in_CTRL:
            if "SIDESLIP_ANGLE" in config:
                beta.append(config.SIDESLIP_ANGLE)
            else:
                beta.append(0.0)
            nBeta = 1
        else:
            if polarSweepType < 4:  # alpha sweep with phi set
                tAlpha = [np.tan(d2r * x) for x in alpha]
                tPhi = [np.tan(d2r * x) for x in phi]
                tb = [x * y for y in tAlpha for x in tPhi]
                beta = [np.arctan(x) / d2r for x in tb]
                nBeta = np.size(beta)
            else:  # Mach ramp
                if "SIDESLIP_ANGLE" in config:
                    beta.append(config.SIDESLIP_ANGLE)
                else:
                    beta.append(0.0)
                nBeta = 1

    if options.verbose:
        print(">>> alpha: " + str(alpha))
        print(">>> beta:  " + str(beta))
        print(">>> phi:   " + str(phi))
        print(">>> Mach   " + str(MachList))

    results.AOA = alpha
    results.MACH = MachList
    results.SIDESLIP_ANGLE = beta

    if options.geomDim == 3:
        results.MOMENT_X = []
        results.MOMENT_Y = []
    if options.Wind:
        results.DRAG = []
        results.LIFT = []
        if options.geomDim == 3:
            results.SIDEFORCE = []
    else:
        results.FORCE_X = []
        results.FORCE_Y = []
        if options.geomDim == 3:
            results.FORCE_Z = []

    results.MOMENT_Z = []

    if polarSweepType == 4:
        outFile = "machRamp_aoa" + str(alpha[0]) + ".dat"
    else:
        outFile = "Polar_M" + str(MachList[0]) + ".dat"
    bufsize = 12
    #
    # ----------- Prepare output header ---------------
    #
    if config.SYSTEM_MEASUREMENTS == "SI":
        length_dimension = "m"
    else:
        length_dimension = "in"
    f = open(outFile, "w", bufsize)
    if options.verbose:
        print("Opening polar sweep file: " + outFile)
    f.write("% \n%  Main coefficients for a polar sweep \n% \n%  ")
    f.write(sweepOption[polarSweepType - 1])
    if polarSweepType == 1:
        satxt = "  ;  Roll angle = %7.2f " % (phi[0])
    elif polarSweepType == 2:
        satxt = "  ;  Sideslip angle = %7.2f " % (beta[0])
    elif polarSweepType == 3:
        satxt = " ;   AOA = %7.2f " % (alpha[0])
    elif polarSweepType == 4:
        satxt = "  ;  AOA = %7.2f Side slip angle = %7.2f ; " % (alpha[0], beta[0])
    f.write(satxt)
    f.write("\n% \n")
    f.write("% ================== Reference parameteres ======================\n%\n")
    XR = config.REF_ORIGIN_MOMENT_X
    YR = config.REF_ORIGIN_MOMENT_Y
    ZR = config.REF_ORIGIN_MOMENT_Z
    f.write("%  Reference point for moments   : [ ")
    line_text = "%s , %s , %s ]  [ %s ] " % (XR, YR, ZR, length_dimension)
    f.write(line_text)
    f.write("\n")
    f.write("%  Reference area and length: ")
    line_text = "Aref : %s Lref : %s  [ %s ] " % (
        config.REF_AREA,
        config.REF_AREA,
        length_dimension,
    )
    f.write(line_text)
    f.write("\n%  ")
    line_text = "Mach : %7.2f  ,  " % (config.MACH_NUMBER)
    f.write(line_text)
    if config.SOLVER == "NAVIER_STOKES":
        line_text = "Reynolds Number  :  %s   " % (config.REYNOLDS_NUMBER)
        f.write(line_text)
        line_text = "Reynolds length :   %s   [ %s ] " % (
            config.REYNOLDS_LENGTH,
            length_dimension,
        )
    else:
        line_text = "Physical problem : %s " % (config.SOLVER)
    f.write(line_text)
    f.write("\n%  ")
    rho = float(config.FREESTREAM_PRESSURE) / (
        float(config.GAS_CONSTANT) * float(config.FREESTREAM_TEMPERATURE)
    )
    line_text = "Reference pressure : %s  ," % (config.FREESTREAM_PRESSURE)
    f.write(line_text)
    line_text = "  Reference density : %7.4f  , " % (rho)
    f.write(line_text)
    line_text = "  Reference Temperature : %s " % (config.FREESTREAM_TEMPERATURE)
    f.write(line_text)
    f.write("\n%  ")
    line_text = "Constant specific heat ratio  : %s  ,   " % (config.GAMMA_VALUE)
    f.write(line_text)
    line_text = "Gas constant :  %s  " % (config.GAS_CONSTANT)
    f.write(line_text)
    f.write("\n%  ")
    line_text = "Grid file :  %s  " % (config.MESH_FILENAME)
    f.write(line_text)
    f.write("\n%  ")
    symmmetry_exists = False
    if "MARKER_SYM" in config:
        if config.MARKER_SYM != "NONE":
            symmmetry_exists = True
    if symmmetry_exists:
        line_text = "Symmetry surface :  yes"
    else:
        line_text = "Symmetry surface :  no"
    f.write(line_text)
    f.write("\n% \n")
    # -----------------   end reference parameter section --------------
    if options.Wind:
        f.write("% AOA,  Mach,       CL,         CD,            ")
        if options.geomDim == 3:
            f.write("CY,           ")
    else:
        if options.geomDim == 2:
            f.write("% AOA,  Mach,       CX,          CY,            ")
        else:
            f.write("% AOA,  Mach,       CX,          CZ,          CY,            ")

    if options.geomDim == 3:
        f.write("Cmx,         Cmz,           Cmy  \n")
    else:
        f.write("        Cmz \n")

    firstSweepPoint = True
    # iterate mach
    for MachNumber in MachList:

        # iterate angles
        for j in range(0, nPolara):
            if polarSweepType < 3:
                AngleAttack = alpha[j]
                SIDESLIP_ANGLE = beta[0]
            elif polarSweepType == 3:
                AngleAttack = alpha[0]
                SIDESLIP_ANGLE = beta[j]
            else:
                AngleAttack = alpha[0]
                SIDESLIP_ANGLE = beta[0]

            if options.verbose:
                print(
                    "Sweep step " + str(j) + ": Mach = " + str(MachNumber) + ", aoa = ",
                    str(AngleAttack) + ", beta = " + str(SIDESLIP_ANGLE),
                )

            # local config and state
            konfig = copy.deepcopy(config)
            # enable restart in polar sweep
            konfig.DISCARD_INFILES = "YES"
            ztate = copy.deepcopy(state)
            #
            # The eval functions below requires definition of various optimization
            # variables, though we are handling here only a direct solution.
            # So, if they are missing in the cfg file (and only then), some dummy values are
            # introduced here
            if "OBJECTIVE_FUNCTION" not in konfig:
                konfig.OBJECTIVE_FUNCTION = "DRAG"
            if "DV_KIND" not in konfig:
                konfig.DV_KIND = ["FFD_SETTING"]
            if "DV_PARAM" not in konfig:
                konfig.DV_PARAM = {"FFDTAG": ["1"], "PARAM": [[0.0, 0.5]], "SIZE": [1]}
            if "DEFINITION_DV" not in konfig:
                konfig.DEFINITION_DV = {
                    "FFDTAG": [[]],
                    "KIND": ["HICKS_HENNE"],
                    "MARKER": [["WING"]],
                    "PARAM": [[0.0, 0.05]],
                    "SCALE": [1.0],
                    "SIZE": [1],
                }
            if "OPT_OBJECTIVE" not in konfig:
                obj = {}
                obj["DRAG"] = {"SCALE": 1.0e-2, "OBJTYPE": "DEFAULT", "MARKER": "None"}
                konfig.OPT_OBJECTIVE = obj
            #
            # --------- end of dummy optimization variables definition section ---------
            #

            # set angle of attack and side-slip angle
            konfig.AOA = AngleAttack
            konfig.SIDESLIP_ANGLE = SIDESLIP_ANGLE
            konfig.MACH_NUMBER = MachNumber
            caseName = "DIRECT_M_" + str(MachNumber) + "_AOA_" + str(AngleAttack)
            print("Mach = ", konfig.MACH_NUMBER, "AOA = ", konfig.AOA)
            print("case :" + caseName)

            if firstSweepPoint:
                # if caseName exists copy the restart file from it for run continuation
                # Continue from previous sweep point if this is not he first
                if os.path.isdir(caseName):
                    command = "cp " + caseName + "/" + config.SOLUTION_FILENAME + " ."
                    if options.verbose:
                        print(command)
                    shutil.copy2(caseName + "/" + config.SOLUTION_FILENAME, os.getcwd())
                    konfig.RESTART_SOL = "YES"
                else:
                    konfig.RESTART_SOL = "NO"
                firstSweepPoint = False
            else:
                konfig.RESTART_SOL = "YES"
            if konfig.RESTART_SOL == "YES":
                ztate.FILES.DIRECT = config.SOLUTION_FILENAME
            # run su2
            if options.Wind:
                drag = SU2.eval.func("DRAG", konfig, ztate)
                lift = SU2.eval.func("LIFT", konfig, ztate)
                if options.geomDim == 3:
                    sideforce = SU2.eval.func("SIDEFORCE", konfig, ztate)
            else:
                force_x = SU2.eval.func("FORCE_X", konfig, ztate)
                force_y = SU2.eval.func("FORCE_Y", konfig, ztate)
                if options.geomDim == 3:
                    force_z = SU2.eval.func("FORCE_Z", konfig, ztate)

            momentz = SU2.eval.func("MOMENT_Z", konfig, ztate)
            if options.geomDim == 3:
                momentx = SU2.eval.func("MOMENT_X", konfig, ztate)
                momenty = SU2.eval.func("MOMENT_Y", konfig, ztate)

            # append results

            if options.Wind:
                results.DRAG.append(drag)
                results.LIFT.append(lift)
                if options.geomDim == 3:
                    results.SIDEFORCE.append(sideforce)
            else:
                results.FORCE_X.append(force_x)
                results.FORCE_Y.append(force_y)
                if options.geomDim == 3:
                    results.FORCE_Z.append(force_z)

            results.MOMENT_Z.append(momentz)
            if options.geomDim == 3:
                results.MOMENT_X.append(momentx)
                results.MOMENT_Y.append(momenty)

            output = "  " + str(AngleAttack) + ",   " + str(MachNumber) + ", "

            if options.Wind:
                output = output + str(lift) + ", " + str(drag)
                if options.geomDim == 3:
                    output = output + ", " + str(sideforce)
            else:
                if options.geomDim == 2:
                    output = output + str(force_x) + ", " + str(force_y)
                else:
                    output = (
                        output
                        + str(force_x)
                        + ", "
                        + str(force_z)
                        + ", "
                        + str(force_y)
                    )
            if options.geomDim == 3:
                output = output + ", " + str(momentx) + ", " + str(momentz) + ", "
                output = output + str(momenty) + " \n"
            else:
                output = output + ", " + str(momentz) + " \n"

            f.write(output)
            # save data
            SU2.io.save_data("results.pkl", results)
            shutil.copy2("results.pkl", "DIRECT")
            shutil.copy2(config.SOLUTION_FILENAME, "DIRECT")

            if os.path.isdir(caseName):
                command = (
                    "cat "
                    + caseName
                    + "/history_direct.dat DIRECT/history_direct.dat > tmp && mv tmp "
                    + "DIRECT/history_direct.dat"
                )
                if options.verbose:
                    print(command)
                os.system(command)
                shutil.rmtree(caseName)

            command = "cp -p -R DIRECT " + caseName
            if options.verbose:
                print(command)
            shutil.copytree("DIRECT", caseName)

    # Close open file
    f.close()
    if os.path.isdir("DIRECT"):
        shutil.rmtree("DIRECT")
    if os.path.isfile(config.SOLUTION_FILENAME):
        os.remove(config.SOLUTION_FILENAME)
    if os.path.isfile("results.pkl"):
        os.remove("results.pkl")
    print("Post sweep cleanup completed")

    #         sys.exit(0)

    # ----------------------------------------------------------#

    #: for each angle

    # plotting
    # plt.figure()
    # plt.plot( results.MACH_NUMBER, results.AOA , results.LIFT , results.DRAG )
    # plt.show()


if __name__ == "__main__":
    main()
