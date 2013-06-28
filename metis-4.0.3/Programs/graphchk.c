/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * graphchk.c
 *
 * This file checks the validity of a graph
 *
 * Started 8/28/94
 * George
 *
 * $Id: graphchk.c,v 1.1 1998/11/27 17:59:33 karypis Exp $
 *
 */

#include <metis.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  GraphType graph;
  char filename[256];
  int wgtflag;

  if (argc != 2) {
    printf("Usage: %s <GraphFile>\n", argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);

  ReadGraph(&graph, filename, &wgtflag);
  if (graph.nvtxs == 0) {
    printf("Empty graph!\n");
    exit(0);
  }

  printf("**********************************************************************\n");
  printf("%s", METISTITLE);
  printf("Graph Information ---------------------------------------------------\n");
  printf("  Name: %s, #Vertices: %d, #Edges: %d\n\n", filename, graph.nvtxs, graph.nedges/2);
  printf("Checking Graph... ---------------------------------------------------\n");

  if (CheckGraph(&graph))
    printf("   The format of the graph is correct!\n");
  else
    printf("   The format of the graph is incorrect!\n");

  printf("\n**********************************************************************\n");


  GKfree(&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
}  


