# \file polarSweepLib.py
#  \brief Functions library for compute_polar.py script.
#  \author E Arad
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

from __future__ import division, print_function, absolute_import
from numpy import *


def parLocator(keyWord, b, n, iDoNot, verbose):

    # ---- -- locate the relevant line in base input file
    # --- do not select line iDoNot (unless it is -1)
    #
    keyWord = keyWord.lower()
    iFocus = -1
    icol = -1
    for i in range(1, n):

        lineString = str(b[i]).lower()
        # check if : exist in line
        try:
            icol = lineString.index(":")
        except ValueError:
            pass  # do nothing
        if icol > -1:
            # verify that this line was not commented out
            try:
                ii = lineString[: icol - 1].index("#")
                pass  # do nothing
            except ValueError:
                # This line wasn't commented out
                try:
                    ii = lineString.index(keyWord)
                    if i != iDoNot:
                        # string.index and not string.find is used here, since index raises
                        # exception when search is failed
                        if verbose:
                            print("parLocator: " + str(i) + " found:  " + str(b[i]))
                        iFocus = i
                        break
                    else:
                        iFocus = -1
                except ValueError:
                    pass  # do nothing

    if iFocus == -1:
        if verbose:
            print("parLocator: Keyword ->" + str(keyWord) + "<-  not found")
    return iFocus


def stringLocator(keyWord, b, n, verbose):

    # ---- -- locate the relevant line in a file
    #
    keyWord = keyWord.lower()
    iFocus = -1
    for i in range(1, n):
        lineString = str(b[i]).lower()
        try:
            ii = lineString.index(keyWord)
            if verbose:
                print("parLocator: " + str(i) + " found:  " + str(b[i]))
            iFocus = i
            break
        except ValueError:
            pass  # do nothing

    if iFocus == -1:
        if verbose:
            print("parLocator: Keyword ->" + str(keyWord) + "<-  not found")

    return iFocus


def readList(dataFile, iLine, verbose):

    from numpy import size

    #
    # ----read list from file to a local float list
    #
    listDataLine = dataFile[iLine]
    icol = listDataLine.index(":")
    Data = listDataLine[icol + 1 :]
    lData = Data.split(",")
    nData = size(lData)

    if verbose:
        print("readList nData = " + str(nData))
    fData = map(float, lData)
    return list(fData), nData


def readParameter(dataFile, nLines, keyWord, iDoNot, verbose):

    from numpy import size

    #
    # ----read a parameter from a file-list
    #
    keyWord = keyWord.lower()
    ipar = parLocator(keyWord, dataFile, nLines, iDoNot, verbose)
    if ipar == -1:
        if verbose:
            print(
                " failed to locate " + keyWord + " in base input file; Set value to 1"
            )
        paVal = 1
    else:
        paLine = dataFile[ipar]
        icol = paLine.index(":")
        try:
            iComment = paLine.index("#")
            paVal = paLine[icol + 1 : iComment - 1].lower()
        except ValueError:
            paVal = paLine[icol + 1 :].lower()

    if verbose:
        if ipar != -1:
            print(keyWord + " = " + paVal)

    return paVal, ipar


def setContribution(dataFile, nLines, keyWord, iDoNot, verbose):

    from numpy import size
    import string

    #
    # default values
    #
    nameText = ""
    removeContribution = False
    #
    # ----Determine if a given amily contribute to force
    #
    # Start by locating lines setting contribution
    #
    keyWord = keyWord.lower()
    ipar = parLocator(keyWord, dataFile, nLines, iDoNot, verbose)
    if ipar == -1:
        if verbose:
            print(
                " failed to locate " + keyWord + " in base input file; Set value to 1"
            )
        paVal = 1
    else:
        paLine = dataFile[ipar]
        icol = paLine.index(":")

        # Now identify the first part of this line
        firstPart = paLine[0:icol]

        # now find out where the standard text ends
        iBF = firstPart.lower().index("family") + 6
        nameText = string.join(firstPart[iBF:].split(), "")

        # component name located. Now check about its contribution
        try:
            iComment = paLine.index("#")
            secondPart = paLine[icol + 1 : iComment - 1].lower()
        except ValueError:
            secondPart = paLine[icol + 1 :].lower()

        # find the second colon of this line
        icol2 = secondPart.index(":")
        try:
            iComment = secondPart.index("#")
            yesNoText = secondPart[icol2 + 1 : iComment - 1].lower()
        except ValueError:
            yesNoText = secondPart[icol2 + 1 :].lower()

        try:
            noFound = yesNoText.lower().index("no")
            removeContribution = True
        except ValueError:
            removeContribution = False

    if verbose:
        if ipar != -1:
            print(
                " part: "
                + nameText
                + " remove contribution: "
                + str(removeContribution)
            )

    return nameText, removeContribution, ipar


def setPolaraType(ctrl, nc, verbose):

    # scan the control file and determine polara type and angles
    # Determine pitch direction from control file
    # ---------------------------------------------------

    keyWordPitchAxis = "pitch axis"
    iPA = parLocator(keyWordPitchAxis, ctrl, nc, -1, verbose)

    if iPA == -1:
        PA = "z"  # This is the default
    else:
        paLine = ctrl[iPA]
        icol = paLine.index(":")
        paVal = paLine[icol + 1 :].lower()
        zFound = "z" in paVal
        if zFound:
            PA = "z"
        else:
            yFound = "y" in paVal
            if yFound:
                PA = "y"
            else:
                raise SystemExit(
                    "ERROR in control file: only Y or Z can be given for control keyWord  ->"
                    + keyWordPitchAxis
                    + "<-"
                )

    if verbose:
        print("Pitch axis is " + PA.upper())
    #
    # angles definitions:
    #   alpha ... angle of attack
    #   beta  ... side-slip angle
    #   phi   ... roll angle
    #
    # Note: Actually alpha here is the angle of rotation about the above-defined pitch axis
    #       Thus, by replacing the pitch-axis, all that is said here about alpha is actually for beta
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
    #
    # Now let us find out which angles are specified in the control file, to figure out polarSweepType and polarVar
    #
    keyWordListAOA = "angles of attack"
    iListAOA = parLocator(keyWordListAOA, ctrl, nc, -1, verbose)
    keyWordListPhi = "roll angles"
    iListPhi = parLocator(keyWordListPhi, ctrl, nc, -1, verbose)
    keyWordListBeta = "side slip angle"
    iListBeta = parLocator(keyWordListBeta, ctrl, nc, -1, verbose)
    keyWordListMRN = "mach ramp numbers"
    iListMRN = parLocator(keyWordListMRN, ctrl, nc, -1, verbose)

    #
    # Check first if this is a Mach ramp session
    #
    if iListMRN > -1:
        polarSweepType = 4
        # This is a Mach rmp session
        polarVar = "MachRampNumbers"
        MachList, nMach = readList(ctrl, iListMRN, verbose)
        #
        # Now check if any angle was specified
        #
        if iListBeta == -1:
            nBeta = 0
            beta = []
            velDirOption = 1
            # Velocity dirction vector v(alpha,phi). May be overwritten below
        else:
            beta, nBeta = readList(ctrl, iListBeta, verbose)
            velDirOption = 2
            # Velocity dirction vector v(alpha,beta)
            if nBeta > 1:
                raise SystemExit(
                    "ERROR in control file: >>>>>>>> nBeta > 1 in a Mach Ramp session <<<<<<<<"
                )

        if iListAOA == -1:
            if velDirOption == 2:
                alpha = [0.0]
                nAalpha = 1
            else:
                alpha = []
                nAalpha = 0
                velDirOption = 0
                # No specification of Velocity dirction vector. May be overwritten below
        else:
            alpha, nAalpha = readList(ctrl, iListAOA, verbose)
            if nAalpha > 1:
                raise SystemExit(
                    "ERROR in control file: >>>>>>>>  nAlpha > 1 in a Mach Ramp session  <<<<<<<<"
                )

        if iListPhi == -1:
            if velDirOption != 1:
                phi = []
                nPhi = 0
            else:
                phi = [0.0]
                nPhi = 1
        else:
            phi, nPhi = readList(ctrl, iListPhi, verbose)
            if nPhi > 1:
                raise SystemExit(
                    "ERROR in control file:  >>>>>>>> nPhi > 1 in a Mach Ramp session  <<<<<<<<"
                )
            if velDirOption == 0:
                #    if phi is specified, then this is a alpha,phi case, with alpha = 0
                velDirOption = 1
                alpha = [0.0]
                nAalpha = 1

        if nPhi + nBeta >= 2:
            raise SystemExit(
                "ERROR in control file:  >>>>>>>> Both phi and Beta specified  (in a Mach Ramp session)  <<<<<<<<"
            )

    else:
        #
        # this is not a mach ramp
        MachList = []
        nMach = 0
        #
        # So, this is a polar-sweep and not mach ramp.
        # Now find out polarSweepType (1,2,or 3)
        #

        if iListPhi == -1:
            if iListBeta == -1:
                polarSweepType = 1
                polarVar = "aoa"
                # phi/beta not found. Polar sweep in alpha for phi=beta=0
                velDirOption = 1
                #  Velocity dirction vector v(alpha,phi).
                phi = [0.0]
                nPhi = 1
                nBeta = 0
                beta = []
            else:
                #         beta was found in control file, phi is not there; check how about alpha
                nPhi = 0
                phi = []
                polarSweepType = 2
                polarVar = "aoa"
                velDirOption = 2
                #  Velocity dirction vector v(alpha,beta).
                beta, nBeta = readList(ctrl, iListBeta, verbose)
                if nBeta > 1:
                    raise SystemExit(
                        "ERROR in control file: nBeta > 1. For polar sweep in beta exchange pitch-axis and use aoa"
                    )

            if iListAOA == -1:
                raise SystemExit(
                    "ERROR in control file: phi and alpha are missing. Polar sweep not defined"
                )

            alpha, nAalpha = readList(ctrl, iListAOA, verbose)

        else:
            #
            #      phi was found in control file, so beta must not be there
            #
            if iListBeta > -1:
                raise SystemExit(
                    "ERROR in control file: both phi and beta specified.  Polar sweep not defined "
                )

            nBeta = 0
            beta = []
            velDirOption = 1
            #  Velocity dirction vector v(alpha,phi).
            #
            #     Check now if alpha appears
            #
            if iListAOA == -1:
                #
                #     phi found in control file, but alpha is missing, so it is a polar-sweep in phi with alpha=0
                #
                polarSweepType = 3
                polarVar = "phi"
                alpha = [0.0]
                nAalpha = 1

            else:
                #
                #     Both alpha and phi found in control file. Find out which one is a list
                #
                alpha, nAalpha = readList(ctrl, iListAOA, verbose)

            phi, nPhi = readList(ctrl, iListPhi, verbose)

            if nAalpha == 1:
                if nPhi > 1:
                    polarSweepType = 3
                    polarVar = "phi"
                else:
                    polarSweepType = 1
                    polarVar = "aoa"

                nBeta = 0
                beta = []
            else:
                #
                # -----that is nAlpha > 1
                #
                if nPhi > 1:
                    raise SystemExit(
                        "ERROR in control file: read lists in both alpha and phi.  Polar sweep not defined "
                    )

                polarSweepType = 1
                polarVar = "aoa"
                nBeta = 0
                beta = []

    #
    # Here we end the long if cycle, refrring to Mach ramp or angle sweep

    # -------------------------------------------------------------------------------------------
    if verbose:
        if polarSweepType == 1:
            print(
                "Sweep type: "
                + str(polarSweepType)
                + " in alpha. nAalpha = "
                + str(nAalpha)
                + " phi = "
                + str(phi)
            )
        elif polarSweepType == 2:
            print(
                "Sweep type: "
                + str(polarSweepType)
                + " in alpha. nAalpha = "
                + str(nAalpha)
                + " beta = "
                + str(beta)
            )
        elif polarSweepType == 3:
            print(
                "Sweep type: "
                + str(polarSweepType)
                + " in phi. nPhi = "
                + str(nPhi)
                + " alpha = ",
                str(alpha),
            )
        elif polarSweepType == 4:
            print(
                "Sweep type: " + str(polarSweepType) + " in Mach. nMach = " + str(nMach)
            )

    return (
        PA,
        polarSweepType,
        velDirOption,
        nAalpha,
        nBeta,
        nPhi,
        nMach,
        alpha,
        beta,
        phi,
        MachList,
        polarVar,
    )


def setVelDir(velDirOption, PA, alphar, phir, betar):

    # set the velocity direction
    from numpy import sin, cos, tan, size

    #
    # Check for alpha and if we are dealing with values greater than 88deg (near 90)
    # In such cases cases change to  sin cos formualtion
    #

    a88 = 1.5359  #   88 degrees
    if velDirOption == 2:

        if PA == "z":
            if alphar < a88:
                dv1 = [cos(betar) for x in alphar]
                dv2 = tan(alphar) * cos(betar)
                dv3 = [sin(betar) for x in alphar]
            else:
                dv1 = cos(betar) * cos(alphar)
                dv2 = sin(alphar) * cos(betar)
                dv3 = sin(betar) * cos(alphar)
        else:
            if alphar < a88:
                dv1 = [cos(betar) for x in alphar]
                dv2 = [sin(betar) for x in alphar]
                dv3 = tan(alphar) * cos(betar)
            else:
                dv1 = cos(betar) * cos(alphar)
                dv2 = sin(betar) * cos(alphar)
                dv3 = sin(alphar) * cos(betar)
    else:

        if size(alphar) > size(phir):
            dummyVec = alphar
        else:
            dummyVec = phir

        if PA == "z":
            if alphar < a88:
                dv1 = [1.0 for x in dummyVec]
                dv2 = tan(alphar) * cos(phir)
                dv3 = tan(alphar) * sin(phir)
            else:
                if size(alphar) > 1:
                    dv1 = cos(alphar)
                else:
                    dv1 = [cos(alphar) for x in dummyVec]

                dv2 = sin(alphar) * cos(phir)
                dv3 = sin(alphar) * sin(phir)
        else:
            if alphar < a88:
                dv1 = [1.0 for x in dummyVec]
                dv2 = tan(alphar) * sin(phir)
                dv3 = tan(alphar) * cos(phir)
            else:
                if size(alphar) > 1:
                    dv1 = cos(alphar)
                else:
                    dv1 = [cos(alphar) for x in dummyVec]
                dv2 = sin(alphar) * sin(phir)
                dv3 = sin(alphar) * cos(phir)

    return dv1, dv2, dv3


def processAddAngle(addRunStr, nPolara, parAngle, angleEqualCriterion):
    #
    # ---------------------------------------------------------------------
    # Process the list of interactively added angles.
    # Note that parAngle can receive also MachList
    # --------------------------------------------------------------------

    # ------------ create a list out of entered angles

    addRunList = addRunStr.split(",")
    fAddRunListRaw = map(float, addRunList)
    fAddRunList = sort(fAddRunListRaw)
    nAddRun = size(fAddRunList)

    # ------- By default, do not compute the cases in input file, since they were computed already

    computeCase = [False for j in range(0, nPolara + nAddRun)]
    rerunCase = [False for j in range(0, nPolara + nAddRun)]

    # ----- Now, for each  new angle/Mach verify if it is a rerun or inserted new value

    closestAngle = [0 for i in range(0, nAddRun)]
    for i in range(0, nAddRun):

        diff = [abs(fAddRunList[i] - x) for x in parAngle]
        iClose = diff.index(min(diff))
        closestAngle[i] = parAngle[iClose]
        if min(diff) < angleEqualCriterion:
            computeCase[iClose] = True
            rerunCase[iClose] = True

        else:

            if fAddRunList[i] > closestAngle[i]:
                ii = iClose + 1
            else:
                ii = iClose

            tmpAng1 = parAngle[:ii]
            tmpAng1.append(fAddRunList[i])
            tmpAng2 = parAngle[ii:]
            tmpAng1.extend(tmpAng2)

            parAngle = tmpAng1
            nPolara = nPolara + 1
            computeCase[ii] = True

    return nPolara, parAngle, computeCase, rerunCase


def updatedControlFile(ctrl, nc, parAngle, ctrlFile, verbose):

    # generate a modified control file for case with addRun options

    import os

    #
    # -- get a proper list of updated parameter-angle
    st1 = str(parAngle)
    updatedAngleList = st1[1:-1]
    # Now let us find out which angles are specified in the control file, to figure out polarSweepType and polarVar
    #
    keyWordListAOA = "angles of attack"
    iListAOA = parLocator(keyWordListAOA, ctrl, nc, -1, verbose)
    keyWordListPhi = "roll angles"
    iListPhi = parLocator(keyWordListPhi, ctrl, nc, -1, verbose)
    keyWordListBeta = "side slip angle"
    iListBeta = parLocator(keyWordListBeta, ctrl, nc, -1, verbose)
    keyWordListMRN = "mach ramp numbers"
    iListMRN = parLocator(keyWordListMRN, ctrl, nc, -1, verbose)
    #
    # Check first if this is a Mach ramp session
    #
    if iListMRN > -1:
        polarSweepType = 4
        # This is a Mach rmp session
        polarVar = "MachRampNumbers"
        MachList, nMach = readList(ctrl, iListMRN, verbose)
        #
        # Now check if any angle was specified
        #
        if iListBeta == -1:
            nBeta = 0
            beta = []
            velDirOption = 1
            # Velocity dirction vector v(alpha,phi). May be overwritten below
        else:
            beta, nBeta = readList(ctrl, iListBeta, verbose)
            velDirOption = 2
            # Velocity dirction vector v(alpha,beta)
            if nBeta > 1:
                raise SystemExit(
                    "ERROR in control file: >>>>>>>> nBeta > 1 in a Mach Ramp session <<<<<<<<"
                )

        if iListAOA == -1:
            if velDirOption == 2:
                alpha = [0.0]
                nAalpha = 1
            else:
                alpha = []
                nAalpha = 0
                velDirOption = 0
                # No specification of Velocity dirction vector. May be overwritten below
        else:
            alpha, nAalpha = readList(ctrl, iListAOA, verbose)
            if nAalpha > 1:
                raise SystemExit(
                    "ERROR in control file: >>>>>>>>  nAlpha > 1 in a Mach Ramp session  <<<<<<<<"
                )

            if iListPhi == -1:
                if velDirOption != 1:
                    phi = []
                    nPhi = 0
                else:
                    phi = [0.0]
                    nPhi = 1
            else:
                phi, nPhi = readList(ctrl, iListPhi, verbose)
                if nPhi > 1:
                    raise SystemExit(
                        "ERROR in control file:  >>>>>>>> nPhi > 1 in a Mach Ramp session  <<<<<<<<"
                    )
                if velDirOption == 0:
                    #    if phi is specified, then this is a alpha,phi case, with alpha = 0
                    velDirOption = 1
                    alpha = [0.0]
                    nAalpha = 1

        if nPhi + nBeta >= 2:
            raise SystemExit(
                "ERROR in control file:  >>>>>>>> Both phi and Beta specified  (in a Mach Ramp session)  <<<<<<<<"
            )

        ctrl[iListMRN] = " Mach ramp numbers :  " + updatedAngleList + "\n"

    else:
        #
        # this is not a mach ramp
        MachList = []
        nMach = 0
        if iListPhi == -1:
            if iListBeta == -1:
                polarSweepType = 1
                polarVar = "aoa"
                # phi/beta not found. Polar sweep in alpha for phi=beta=0
                phi = [0.0]
                nPhi = 1
                nBeta = 0
                beta = []
            else:
                #         beta was found in control file, phi is not there; check how about alpha
                nPhi = 0
                phi = []
                polarSweepType = 2
                polarVar = "aoa"
                beta, nBeta = readList(ctrl, iListBeta, verbose)
                if nBeta > 1:
                    raise SystemExit(
                        "ERROR in control file: nBeta > 1. For polar sweep in beta exchange pitch-axis and use aoa"
                    )

            if iListAOA == -1:
                raise SystemExit(
                    "ERROR in control file: phi and alpha are missing. Polar sweep not defined"
                )

            alpha, nAalpha = readList(ctrl, iListAOA, verbose)
            ctrl[iListAOA] = " angles of attack :  " + updatedAngleList + "\n"

        else:
            #      phi was found in control file, so beta must not be there
            if iListBeta > -1:
                raise SystemExit(
                    "ERROR in control file: both phi and beta specified.  Polar sweep not defined "
                )

            nBeta = 0
            beta = []
            #     Check now if alpha appears
            if iListAOA == -1:
                #     phi found in control file, but alpha is missing, so it is a polar-sweep in phi with alpha=0
                polarSweepType = 3
                polarVar = "phi"
                alpha = [0.0]
                nAalpha = 1

            else:
                #
                #     Both alpha and phi found in control file. Find out which one is a list
                #
                alpha, nAalpha = readList(ctrl, iListAOA, verbose)
                if nAalpha > 1:
                    ctrl[iListAOA] = " angles of attack :  " + updatedAngleList + "\n"

            phi, nPhi = readList(ctrl, iListPhi, verbose)
            if nPhi > 1:
                ctrl[iListPhi] = " roll angles :  " + updatedAngleList + "\n"

    # Prepare a backup of control file

    shutil.copy2(ctrlFile, ctrlFile + ".bck")
    #
    # --- Write down the updated file
    fc = open(ctrlFile, "w")
    fc.writelines(ctrl)
    fc.close()

    print(
        "More cases were added. Original ctrl file saved at "
        + ctrlFile
        + ".bck File "
        + ctrlFile
        + " updated"
    )

    return


def retrievePhysicalData(b, n, polarSweepType, verbose):

    # scan the control file and retrieve physical data parameters and their location
    # Included are Mach and reynolds number (non-dim group)
    #              Pref, rho_ref, Tref (ref group)
    # ---------------------------------------------------
    #
    # physical data, needed for Mach ramp
    #
    keyWord = "mach for coefficients"
    MachNumCoef, iparMcoeff = readParameter(b, n, keyWord, -1, verbose)
    keyWord = "mach"
    MachNum, iprMach = readParameter(
        b, n, keyWord, iparMcoeff, verbose
    )  # look for Mach, but avoid Mach for coefficients
    keyWord = "reynolds length (in meter)"
    ReNumRefLength, iprDRe = readParameter(b, n, keyWord, -1, verbose)
    keyWord = "reynolds"
    ReNum, iprRe = readParameter(b, n, keyWord, iprDRe, verbose)
    sNonDimNum = [MachNum, MachNumCoef, ReNum, ReNumRefLength]
    nonDimNum = map(float, sNonDimNum)
    nonDimNumLoc = [iprMach, iparMcoeff, iprRe, iprDRe]
    #
    # the next set of parameters might, or might not appear in base input file
    # If they appear, they should be updated in a Mach ramp. All 3 of them are needed.
    #
    refParNo = 0
    keyWord = "Reference pressure (in Pa)"
    pRef, iprPr = readParameter(b, n, keyWord, -1, verbose)
    if iprPr > -1:
        refParNo = refParNo + 1

    keyWord = "Reference density (in kg/m^3)"
    rhoRef, iprRho = readParameter(b, n, keyWord, -1, verbose)
    if iprRho > -1:
        refParNo = refParNo + 1

    keyWord = "Reference temperature (in K)"
    TRef, iprT = readParameter(b, n, keyWord, -1, verbose)
    if iprT > -1:
        refParNo = refParNo + 1

    if refParNo == 3:
        refParExist = True
    elif refParNo == 0:
        refParExist = False
    else:
        if polarSweepType == 4:
            raise SystemExit(
                "ERROR in control file: in Mach ramp, base file should include (Pr,rho_r,Tr) or none of them"
            )

    if refParExist:
        sRefPar = [pRef, rhoRef, TRef]
        refPar = map(float, sRefPar)
        refParLoc = [iprPr, iprRho, iprT]
    #
    # Thermodynamic properties
    #
    keyWord = "Constant specific heat ratio"
    gamma, iprGamma = readParameter(b, n, keyWord, -1, verbose)
    keyWord = "Gas constant (J/(kg K))"
    rGas, iprGasC = readParameter(b, n, keyWord, -1, verbose)
    keyWord = "Free stream temperature (in K)"
    TFreeS, iprTFreeS = readParameter(b, n, keyWord, -1, verbose)

    sThermoPar = [gamma, rGas, TFreeS]
    thermoPar = map(float, sThermoPar)
    thermoParLoc = [iprGamma, iprGasC, iprTFreeS]

    if verbose:
        print("base case parameters of Mach ramp")
        print("---------------------------------")
        print(" M = " + sNonDimNum[0] + " Reynolds = " + sNonDimNum[2])
        if refParExist:
            print(
                " Pref = "
                + str(refPar[0])
                + " rhor = "
                + str(refPar[1])
                + " Tr = "
                + str(refPar[2])
            )
            print(
                " gamma = "
                + str(thermoPar[0])
                + " Gas Const = "
                + str(thermoPar[1])
                + " T_freeStream = "
                + str(thermoPar[2])
            )

    return (
        nonDimNum,
        nonDimNumLoc,
        refParExist,
        refPar,
        refParLoc,
        thermoPar,
        thermoParLoc,
    )


def fMachIsentropic(Mach, Gamma):

    # Isentropic relation of Mach
    # ---------------------------------------------------
    #
    fMach = 1.0 + (Gamma - 1.0) / 2.0 * Mach * Mach
    return fMach


#
#
def extractUy(filename, outFile, inDepVar, depVar, verbose):

    import os
    import sys

    # --------------- read the  file

    fc = open(filename, "r")
    data = fc.readlines()
    nc = size(data)
    fc.close()
    print(str(nc) + " lines were written from file " + filename + ". File closed")

    # --------------Retreive the variables names in the Tecplot file

    ivb = stringLocator("VARIABLES", data, nc, verbose)
    if ivb == -1:
        raise SystemExit("ERROR: failed to trace VARIABLES list in input file")

    izo = stringLocator("ZONE", data, nc, verbose)
    if izo == -1:
        raise SystemExit("ERROR: failed to trace ZONE list in input file")

    izo = izo - 1  #  last variables line
    print("list of variables traced between lines " + str(ivb) + " and ", str(izo))

    varListLines = data[ivb:izo]
    nV = len(varListLines)
    varList = []
    iX = -1
    iY = -1
    for i in range(0, nV):
        i1 = varListLines[i].index('"') + 1
        i2 = varListLines[i].rindex('"')
        varList.append(varListLines[i][i1:i2])
        if iX == -1:
            try:
                ifound = varList[i].index(inDepVar)
                iX = i
            except ValueError:
                pass  # do nothing
        if iY == -1:
            try:
                ifound = varList[i].index(depVar)
                iY = i
            except ValueError:
                pass  # do nothing

    print(
        "inDepVar: "
        + inDepVar
        + " : "
        + str(iX + 1)
        + " . DepVar: "
        + depVar
        + " : "
        + str(iY + 1)
        + " of "
        + str(nV)
        + " variables"
    )

    # find out how many nodes

    inodes = stringLocator("Nodes", data, nc, verbose)
    if inodes == -1:
        raise SystemExit("ERROR: failed to trace nodes in input file")

    i1 = data[inodes].index("=") + 1
    i2 = data[inodes].index(",")
    Nodes = int(data[inodes][i1:i2])
    print("Nodes = ", str(Nodes))
    #
    # now map the whole matrix
    #
    i1 = inodes + 3
    i2 = i1 + Nodes
    X = []
    Y = []
    for i in range(i1 + 1, i2):
        ff = map(float, data[i][1:-1].split(" "))
        X.append(ff[iX])
        Y.append(ff[iY])

    nP = len(X)

    # ------       sorting by X

    ind = lexsort((Y, X))
    Xs = take(X, ind)
    Ys = take(Y, ind)
    # write down to a simple 2-columns file
    foc = open(outFile, "w")
    fileHeader = "      " + inDepVar + "                     " + depVar
    foc.write("% " + fileHeader + " \n% -----------------------------------\n%\n")
    for i in range(0, nP):
        Line1 = "   %10.5f    %14.5g  " % (Xs[i], Ys[i]) + "  \n"
        foc.write(Line1)
    foc.close()


#    numpy.plot(Xs,Ys,"-b")


#
def loadArray(Fin, nCol):
    #
    # load a polar-sweep file as an array
    #
    f = open(Fin, "r")
    b = f.readlines()
    n = size(b)
    f.close()
    #
    data = []
    nd = 0
    for i in range(0, n):
        sline = (
            b[i]
            .replace("   ", " ")
            .replace("  ", " ")
            .replace("  ", " ")
            .replace("  ", " ")
            .strip()
            .split(" ")
        )
        sv = size(sline)
        if sv == nCol:
            try:
                dd = map(float, sline)
                data.append(dd)
                nd = nd + 1
            except ValueError:
                pass  # do nothing

    return data, nd


def locateSteps(d, nd, nCol):
    #
    # read polarsweep files and identify steps
    #
    eps = 0.001
    nColD = nCol - 2  # cxbase and quality are not checked
    a = array(d)
    dx = diff(a[:, 0], n=1, axis=0)
    nStairs = []
    for ic in range(1, nColD):
        dy = diff(a[:, ic], n=1, axis=0)
        dydx = dy / dx
        adydx = abs(dydx)
        madydx = adydx.mean(axis=0)
        mmxadydx = max(adydx)
        mmnadydx = min(adydx)
        madydx2 = (mmxadydx + mmnadydx) / 2
        iic = where(adydx < eps * madydx)
        iic2 = where(adydx < eps * madydx2)
        nStairs.append(size(iic) + size(iic2))

    nStM = max(nStairs)
    if nStM > 0:
        fst = open("stairs", "w")
        fst.write("%  \n%     Polara stairs report  \n%  \n")
        fst.write(
            "%  Note that the identified number might be something between the correct number of stairs \n"
        )
        fst.write(
            "%   and 2X this number since since 2 criteria are added in the search script \n \n "
        )
        Headers = ["CX ", "CY ", "CZ ", "Cmx", "Cmy", "Cmz"]
        for ic in range(0, nColD - 1):
            refLine1 = Headers[ic] + ": Number of stairs identified: %i \n " % (
                nStairs[ic]
            )
            fst.write(refLine1)

        fst.close()
    return nStairs, nStM


def find_index(ar, eps):
    #
    # locate array components that are > eps
    #
    ia = []
    for i, v in enumerate(ar):
        if v > eps:
            ia.append(i)
    return ia


def testComponentSum(cbdOutput, verbose):
    #
    # check cbd summation
    #

    coeffNames = ["Cfx", "Cfy", "Cfz", "Cmx", "Cmy", "Cmz"]
    try:
        fd = open(cbdOutput, "r")
        d = fd.readlines()
        nd = size(d)
        fd.close()
        if verbose:
            print("CBD file " + cbdOutput + " loaded by testComponentSum")

    except IOError:
        raise SystemExit("testComponentSum: Failed to find file " + cbdOutput)

    # now read the numerical values from the cdb file

    data, nd = loadArray(cbdOutput, 6)
    # transpose the array
    td = zip(*data)
    # now check correct som for each variable
    eps = 0.01
    errorA = []
    sumD = []
    for i in range(0, 6):
        fsumD = sum(td[i][:-1])
        sumD.append(fsumD)
        if abs(td[i][nd - 1]) > eps:
            error = abs((fsumD - td[i][nd - 1]) / td[i][nd - 1])
        else:
            error = 0
        errorA.append(error)

    iErr = find_index(errorA, 0.005)
    nER = size(iErr)
    if nER > 0:
        print("testComponentSum: Error is components sumation in file " + cbdOutput)
        for i in range(0, nER):
            print(
                "Error found in "
                + coeffNames[iErr[i]]
                + " Error = "
                + str(100 * errorA[iErr[i]])
                + " %"
            )

        corrDataLine = "  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  " % (
            sumD[0],
            sumD[1],
            sumD[2],
            sumD[3],
            sumD[4],
            sumD[5],
        )
    else:
        corrDataLine = " "
    return nER, corrDataLine


def retreiveNumPar(ctrl, nc, keyWord, parType, verbose):
    # get the parameter from the control file. Set it to unity if not found
    # parType: 1  -> integer   2 -> float
    #
    ipar = parLocator(keyWord, ctrl, nc, -1, verbose)
    if ipar == -1:
        #  default value
        if parType == 1:
            parVal = 1
        else:
            parVal = 1.0
    else:
        PARLine = ctrl[ipar]
        icol = PARLine.index(":")
        if parType == 1:
            parVal = int(PARLine[icol + 1 :])
        else:
            parVal = float(PARLine[icol + 1 :])
    return parVal


# -----------------------------------------------------


def loadData(filename, delim):
    # read a 2D data from a file, separated by delim
    # (may be , (comma) or ' ' (space )
    #
    # do array(dout)  (in calling) to obtain result as an array (numpy imported)
    # dout=loadData(filename,delim)
    # v=array(dout)

    import csv

    #    import numpy

    data = []
    with open(filename, "rb") as f:
        # -avoid NULL  error

        reader = csv.reader(
            (
                line.replace("\0", "").replace("   ", " ").replace("  ", " ").strip()
                for line in f
            ),
            delimiter=delim,
        )

        data = [" "]
        for row in reader:
            try:

                data.append(map(float, row))
            except ValueError:
                print("Line doesnt match map float: ")
                print(row)
    # check square matrix
    N1 = len(data[0])
    dout = []
    for i in range(1, len(data)):
        if N1 <= 1:
            N1 = len(data[i])
        if len(data[i]) != N1:
            print("WARNING: Line " + str(i) + ": size does not match. Skipped")
        else:
            dout.append(data[i])

    # adout=array(dout)

    return dout
