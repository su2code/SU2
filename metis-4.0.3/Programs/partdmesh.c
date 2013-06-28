/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * partdmesh.c
 *
 * This file reads in the element node connectivity array of a mesh and 
 * partitions both the elements and the nodes using KMETIS on the dual graph.
 *
 * Started 9/29/97
 * George
 *
 * $Id: partdmesh.c,v 1.1 1998/11/27 17:59:38 karypis Exp $
 *
 */

#include <metis.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  int i, j, ne, nn, etype, numflag=0, nparts, edgecut;
  idxtype *elmnts, *epart, *npart;
  timer IOTmr, DUALTmr;
  char etypestr[4][5] = {"TRI", "TET", "HEX", "QUAD"};
  GraphType graph;

  if (argc != 3) {
    printf("Usage: %s <meshfile> <nparts>\n",argv[0]);
    exit(0);
  }

  nparts = atoi(argv[2]);
  if (nparts < 2) {
    printf("nparts must be greater than one.\n");
    exit(0);
  }
   
  cleartimer(IOTmr);
  cleartimer(DUALTmr);

  starttimer(IOTmr);
  elmnts = ReadMesh(argv[1], &ne, &nn, &etype);
  stoptimer(IOTmr);

  epart = idxmalloc(ne, "main: epart");
  npart = idxmalloc(nn, "main: npart");

  printf("**********************************************************************\n");
  printf("%s", METISTITLE);
  printf("Mesh Information ----------------------------------------------------\n");
  printf("  Name: %s, #Elements: %d, #Nodes: %d, Etype: %s\n\n", argv[1], ne, nn, etypestr[etype-1]);
  printf("Partitioning Dual Graph... ------------------------------------------\n");


  starttimer(DUALTmr);
  METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
  stoptimer(DUALTmr);

  printf("  %d-way Edge-Cut: %7d, Balance: %5.2f\n", nparts, edgecut, ComputeElementBalance(ne, nparts, epart));

  starttimer(IOTmr);
  WriteMeshPartition(argv[1], nparts, ne, epart, nn, npart);
  stoptimer(IOTmr);


  printf("\nTiming Information --------------------------------------------------\n");
  printf("  I/O:          \t\t %7.3f\n", gettimer(IOTmr));
  printf("  Partitioning: \t\t %7.3f\n", gettimer(DUALTmr));
  printf("**********************************************************************\n");

/*
  graph.nvtxs = nn;
  graph.xadj = idxmalloc(nn+1, "xadj");
  graph.vwgt = idxsmalloc(nn, 1, "vwgt");
  graph.adjncy = idxmalloc(20*nn, "adjncy");
  graph.adjwgt = idxsmalloc(20*nn, 1, "adjncy");

  METIS_MeshToNodal(&ne, &nn, elmnts, &etype, &numflag, graph.xadj, graph.adjncy);

  ComputePartitionInfo(&graph, nparts, npart);

  GKfree(&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
*/

  GKfree(&elmnts, &epart, &npart, LTERM);

}


