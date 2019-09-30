/*
 * Example C++ program to write a binary data file for Tecplot.
 * This example does the following:
 *
 *   1.  Opens a datafile called "icosahedrons.plt"
 *   2.  Assigns values for x, y, z and p.
 *   3.  Writes out 6 polyhedral zones.
 *   4.  Close the datafile.
 *
 * All zones share connectivity (facemap) and variables x and p.
 */

#include "TECIO.h"

#include <cmath>

#define PHI 1.61803f // 0.5 + sqrt(1.25). Used in the definition of icosahedra
#define PI 3.14159f
#define NULL 0

int main(int argc, char** argv)
{
    INTEGER4 debug = 1;
    INTEGER4 fileIsDouble = 0;
    INTEGER4 fileType = 0;
    INTEGER4 fileFormat = 0; // 0 == PLT, 1 == SZPLT; Only PLT is currently supported
    INTEGER4 i;

    i = TECINI142((char*)"Icosahedral zones with shared connectivity",
        (char*)"X Y Z P",
        (char*)"icosahedrons.plt",
        (char*)".",
        &fileFormat,
        &fileType,
        &debug,
        &fileIsDouble);

    INTEGER4 zoneType = 7;     // FEPolyhedron
    INTEGER4 numNodes = 12;    // Number of nodes in the zone.
    INTEGER4 numElems = 1;     // Number of icosahedral elements.
    INTEGER4 numFaces = 20;    // Number of faces in the icosahedron.
    INTEGER4 iCellMax = 0;     // Not Used
    INTEGER4 jCellMax = 0;     // Not Used
    INTEGER4 kCellMax = 0;     // Not Used
    double   solTime = 0.0;
    INTEGER4 strandID = 0;     // Static Zone
    INTEGER4 parentZn = 0;     // No Parent Zone
    INTEGER4 isBlock = 1;      // Block
    INTEGER4 nFConns = 0;
    INTEGER4 fNMode = 0;
    INTEGER4 numFaceNodes = 3 * numFaces; // All faces are triangular, so 3 nodes per face
    INTEGER4 numBFaces = 0;
    INTEGER4 numBConnections = 0;
    INTEGER4 shrConn = 0;

    i = TECZNE142((char*)"Icosahedral Zone",
        &zoneType,
        &numNodes,
        &numElems,
        &numFaces,
        &iCellMax,
        &jCellMax,
        &kCellMax,
        &solTime,
        &strandID,
        &parentZn,
        &isBlock,
        &nFConns,
        &fNMode,
        &numFaceNodes,
        &numBFaces,
        &numBConnections,
        NULL,
        NULL,  // All nodal variables
        NULL,
        &shrConn);

    float x[] = { 4, 4, 4, 4, 3, 5, 3, 5, 4 - PHI, 4 + PHI, 4 - PHI, 4 + PHI };
    float y[] = { -1, 1, -1, 1, -PHI, -PHI, PHI, PHI, 0, 0, 0, 0 };
    float z[] = { -PHI, -PHI, PHI, PHI, 0, 0, 0, 0, -1, -1, 1, 1 };
    float p[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
    INTEGER4 dIsDouble = 0;

    i = TECDAT142(&numNodes, &x[0], &dIsDouble);
    i = TECDAT142(&numNodes, &y[0], &dIsDouble);
    i = TECDAT142(&numNodes, &z[0], &dIsDouble);
    i = TECDAT142(&numNodes, &p[0], &dIsDouble);

    INTEGER4 faceNodes[] = {
        1, 2, 10,
        2, 8, 10,
        8, 12, 10,
        12, 6, 10,
        6, 1, 10,
        1, 9, 2,
        2, 9, 7,
        2, 7, 8,
        8, 7, 4,
        8, 4, 12,
        12, 4, 3,
        12, 3, 6,
        6, 3, 5,
        6, 5, 1,
        1, 5, 9,
        9, 5, 11,
        7, 9, 11,
        4, 7, 11,
        3, 4, 11,
        5, 3, 11 };
    INTEGER4 faceLeftElems[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    INTEGER4 faceRightElems[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    INTEGER4 faceNodeCounts[] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

    i = TECPOLYFACE142(
        &numFaces,
        faceNodeCounts,
        &faceNodes[0],
        &faceLeftElems[0],
        &faceRightElems[0]);

    // Rotate in the X-Z plane and create more zones that share connectivity, y, and p with the first zone.
    shrConn = 1;
    int shareVar[4] = { 0, 1, 0, 1 };
    float x1[12];
    float z1[12];
    for (int angle = 60; angle < 360; angle += 60)
    {
        i = TECZNE142((char*)"Icosahedral Zone",
            &zoneType,
            &numNodes,
            &numElems,
            &numFaces,
            &iCellMax,
            &jCellMax,
            &kCellMax,
            &solTime,
            &strandID,
            &parentZn,
            &isBlock,
            &nFConns,
            &fNMode,
            &numFaceNodes,
            &numBFaces,
            &numBConnections,
            NULL,
            NULL,  // All nodal variables
            shareVar,
            &shrConn);
        for (int n = 0; n < 12; ++n)
        {
            x1[n] = (float)(x[n] * cos(angle * PI / 180.0) - z[n] * sin(angle * PI / 180.0));
            z1[n] = (float)(x[n] * sin(angle * PI / 180.0) + z[n] * cos(angle * PI / 180.0));
        }

        i = TECDAT142(&numNodes, &x1[0], &dIsDouble);
        i = TECDAT142(&numNodes, &z1[0], &dIsDouble);
    }

    i = TECEND142();

    return 0;
}

