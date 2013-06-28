/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * io.c
 *
 * This file contains routines related to I/O
 *
 * Started 8/28/94
 * George
 *
 * $Id: io.c,v 1.1 1998/11/27 17:59:34 karypis Exp $
 *
 */

#include <metis.h>



/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
void ReadGraph(GraphType *graph, char *filename, int *wgtflag)
{
  int i, j, k, l, fmt, readew, readvw, ncon, edge, ewgt;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt;
  char *line, *oldstr, *newstr;
  FILE *fpin;

  InitGraph(graph);

  line = (char *)malloc(sizeof(char)*(MAXLINE+1));

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  while (fgets(line, MAXLINE, fpin) && line[0] == '%');

  if (feof(fpin)) {
    graph->nvtxs = 0;
    free(line);
    return;
  }

  fmt = ncon = 0;
  sscanf(line, "%d %d %d %d", &(graph->nvtxs), &(graph->nedges), &fmt, &ncon);

  readew = (fmt%10 > 0);
  readvw = ((fmt/10)%10 > 0);
  if (fmt >= 100) {
    printf("Cannot read this type of file format!");
    exit(0);
  }


  *wgtflag = 0;
  if (readew)
    *wgtflag += 1;
  if (readvw)
    *wgtflag += 2;

  if (ncon > 0 && !readvw) {
    printf("------------------------------------------------------------------------------\n");
    printf("***  I detected an error in your input file  ***\n\n");
    printf("You specified ncon=%d, but the fmt parameter does not specify vertex weights\n", ncon);
    printf("Make sure that the fmt parameter is set to either 10 or 11.\n");
    printf("------------------------------------------------------------------------------\n");
    exit(0);
  }

  graph->nedges *=2;
  ncon = graph->ncon = (ncon == 0 ? 1 : ncon);

  /*printf("%d %d %d %d %d [%d %d]\n", fmt, fmt%10, (fmt/10)%10, ncon, graph->ncon, readew, readvw);*/

  if (graph->nvtxs > MAXIDX) 
    errexit("\nThe matrix is too big: %d [%d %d]\n", graph->nvtxs, MAXIDX, sizeof(idxtype));

  xadj = graph->xadj = idxsmalloc(graph->nvtxs+1, 0, "ReadGraph: xadj");
  adjncy = graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: adjncy");

  vwgt = graph->vwgt = (readvw ? idxmalloc(ncon*graph->nvtxs, "ReadGraph: vwgt") : NULL);
  adjwgt = graph->adjwgt = (readew ? idxmalloc(graph->nedges, "ReadGraph: adjwgt") : NULL);

  /* Start reading the graph file */
  for (xadj[0]=0, k=0, i=0; i<graph->nvtxs; i++) {
    while ((oldstr = fgets(line, MAXLINE, fpin)) && line[0] == '%'); /* skip lines with '#' */
    if (oldstr == NULL && feof(fpin)) {
      printf("------------------------------------------------------------------------------\n");
      printf("Premature end of input file when reading the adjacency list for vertex %d\n", i+1);
      printf("------------------------------------------------------------------------------\n");
      exit(0);
    }
    oldstr = line;
    newstr = NULL;

    if (strlen(line) == MAXLINE) 
      errexit("\nBuffer for fgets not big enough!\n");

    if (readvw) {
      for (l=0; l<ncon; l++) {
        vwgt[i*ncon+l] = (int)strtol(oldstr, &newstr, 10);
        oldstr = newstr;
      }
    }

    for (;;) {
      edge = (int)strtol(oldstr, &newstr, 10) -1;
      oldstr = newstr;

      if (readew) {
        ewgt = (int)strtol(oldstr, &newstr, 10);
        oldstr = newstr;
      }

      if (edge < 0)
        break;

      adjncy[k] = edge;
      if (readew) 
        adjwgt[k] = ewgt;
      k++;
    } 
    xadj[i+1] = k;
  }

  fclose(fpin);

  if (k != graph->nedges) {
    printf("------------------------------------------------------------------------------\n");
    printf("***  I detected an error in your input file  ***\n\n");
    printf("In the first line of the file, you specified that the graph contained\n%d edges. However, I only found %d edges in the file.\n", graph->nedges/2, k/2);
    if (2*k == graph->nedges) {
      printf("\n *> I detected that you specified twice the number of edges that you have in\n");
      printf("    the file. Remember that the number of edges specified in the first line\n");
      printf("    counts each edge between vertices v and u only once.\n\n");
    }
    printf("Please specify the correct number of edges in the first line of the file.\n");
    printf("------------------------------------------------------------------------------\n");
    exit(0);
  }

  free(line);
}



/*************************************************************************
* This function writes out the partition vector
**************************************************************************/
void WritePartition(char *fname, idxtype *part, int n, int nparts)
{
  FILE *fpout;
  int i;
  char filename[256];

  sprintf(filename,"%s.part.%d",fname, nparts);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<n; i++)
    fprintf(fpout,"%d\n",part[i]);

  fclose(fpout);

}


/*************************************************************************
* This function writes out the partition vectors for a mesh
**************************************************************************/
void WriteMeshPartition(char *fname, int nparts, int ne, idxtype *epart, int nn, idxtype *npart)
{
  FILE *fpout;
  int i;
  char filename[256];

  sprintf(filename,"%s.epart.%d",fname, nparts);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<ne; i++)
    fprintf(fpout,"%d\n", epart[i]);

  fclose(fpout);

  sprintf(filename,"%s.npart.%d",fname, nparts);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<nn; i++)
    fprintf(fpout,"%d\n", npart[i]);

  fclose(fpout);


}



/*************************************************************************
* This function writes out the partition vector
**************************************************************************/
void WritePermutation(char *fname, idxtype *iperm, int n)
{
  FILE *fpout;
  int i;
  char filename[256];

  sprintf(filename,"%s.iperm",fname);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the permutation file: %s", filename);

  for (i=0; i<n; i++)
    fprintf(fpout,"%d\n", iperm[i]);

  fclose(fpout);

}



/*************************************************************************
* This function checks if a graph is valid
**************************************************************************/
int CheckGraph(GraphType *graph)
{
  int i, j, k, l, nvtxs, err=0;
  idxtype *xadj, *adjncy, *adjwgt;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;


  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];

      if (i == k) {
        printf("Vertex %d contains a self-loop (i.e., diagonal entry in the matrix)!\n", i);
        err++;
      }
      else {
        for (l=xadj[k]; l<xadj[k+1]; l++) {
          if (adjncy[l] == i) {
            if (adjwgt != NULL && adjwgt[l] != adjwgt[j]) {
              printf("Edges (%d %d) and (%d %d) do not have the same weight! %d %d\n", i,k,k,i, adjwgt[l], adjwgt[adjncy[j]]);
              err++;
            }
            break;
          }
        }
        if (l == xadj[k+1]) {
          printf("Missing edge: (%d %d)!\n", k, i);
          err++;
        }
      }
    }
  }

  if (err > 0) 
    printf("A total of %d errors exist in the input file. Correct them, and run again!\n", err);

  return (err == 0 ? 1 : 0);
}


/*************************************************************************
* This function reads the element node array of a mesh
**************************************************************************/
idxtype *ReadMesh(char *filename, int *ne, int *nn, int *etype)
{
  int i, j, k, esize;
  idxtype *elmnts;
  FILE *fpin;

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  if (fscanf(fpin, "%d %d", ne, etype) != 2) {
    printf("Header line of input file does not contain two numbers.\n");
    exit(0);
  }

  switch (*etype) {
    case 1:
      esize = 3;
      break;
    case 2:
      esize = 4;
      break;
    case 3:
      esize = 8;
      break;
    case 4:
      esize = 4;
      break;
    default:
      errexit("Unknown mesh-element type: %d\n", *etype);
  }

  elmnts = idxmalloc(esize*(*ne), "ReadMesh: elmnts");

  for (j=esize*(*ne), i=0; i<j; i++) {
    if (fscanf(fpin, "%d", elmnts+i) != 1) {
      printf("Missing node number %d for element %d\n", i%esize+1, i/esize);
      exit(0);
    }
    elmnts[i]--;
  }

  fclose(fpin);

  *nn = elmnts[idxamax(j, elmnts)]+1;

  return elmnts;
}


/*************************************************************************
* This function writes a graphs into a file 
**************************************************************************/
void WriteGraph(char *filename, int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, j;
  FILE *fpout;

  if ((fpout = fopen(filename, "w")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fprintf(fpout, "%d %d", nvtxs, xadj[nvtxs]/2);

  for (i=0; i<nvtxs; i++) {
    fprintf(fpout, "\n");
    for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fpout, " %d", adjncy[j]+1);
  }

  fclose(fpout);
}


/*************************************************************************
* This function writes a graphs into a file 
**************************************************************************/
void WriteMocGraph(GraphType *graph)
{
  int i, j, nvtxs, ncon;
  idxtype *xadj, *adjncy;
  float *nvwgt;
  char filename[256];
  FILE *fpout;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  nvwgt = graph->nvwgt;

  sprintf(filename, "moc.graph.%d.%d", nvtxs, ncon);

  if ((fpout = fopen(filename, "w")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fprintf(fpout, "%d %d 10 1 %d", nvtxs, xadj[nvtxs]/2, ncon);

  for (i=0; i<nvtxs; i++) {
    fprintf(fpout, "\n");
    for (j=0; j<ncon; j++)
      fprintf(fpout, "%d ", (int)((float)10e6*nvwgt[i*ncon+j]));

    for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fpout, " %d", adjncy[j]+1);
  }

  fclose(fpout);
}
