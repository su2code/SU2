/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mtest.c
 *
 * This file is a comprehensive tester for all the graph partitioning/ordering
 * routines of METIS
 *
 * Started 9/18/98
 * George
 *
 * $Id: mtest.c,v 1.8 1998/09/20 17:36:14 karypis Exp $
 *
 */

#include <metis.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  int i, nparts, options[10];
  idxtype *part;
  float lbvec[MAXNCON], rubvec[MAXNCON];
  GraphType graph;
  int numflag = 0, wgtflag = 0, edgecut;


  if (argc != 2) {
    printf("Usage: %s <GraphFile>\n",argv[0]);
    exit(0);
  }
    
  ReadGraph(&graph, argv[1], &wgtflag);

  printf("**********************************************************************\n");
  printf("%s", METISTITLE);
  printf("Graph Information ---------------------------------------------------\n");
  printf("  Name: %s, #Vertices: %d, #Edges: %d\n", argv[1], graph.nvtxs, graph.nedges/2);

  Test_PartGraph(graph.nvtxs, graph.xadj, graph.adjncy);

  Test_PartGraphV(graph.nvtxs, graph.xadj, graph.adjncy);

  Test_PartGraphmC(graph.nvtxs, graph.xadj, graph.adjncy);

  Test_ND(graph.nvtxs, graph.xadj, graph.adjncy);

  printf("\n---------------------------------------------------------------------\n");
  printf(" Testing completed.\n");
  printf("**********************************************************************\n");

  GKfree(&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
}  



/*************************************************************************
* This function tests the regular graph partitioning routines
**************************************************************************/
void Test_PartGraph(int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, j, jj, k, tstnum, rcode;
  int nparts, numflag, wgtflag, edgecut, options[10];
  idxtype *part, *vwgt, *adjwgt;
  float tpwgts[256];

  vwgt = idxmalloc(nvtxs, "vwgt");
  for (i=0; i<nvtxs; i++)
    vwgt[i] = RandomInRange(10);

  adjwgt = idxmalloc(xadj[nvtxs], "adjwgt");
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (i < k) {
        adjwgt[j] = 1+RandomInRange(5);
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          if (adjncy[jj] == i) {
            adjwgt[jj] = adjwgt[j];
            break;
          }
        }
      }
    }
  }


  part = idxmalloc(nvtxs, "part");

  tpwgts[0] = .1;
  tpwgts[1] = .2;
  tpwgts[2] = .3;
  tpwgts[3] = .1;
  tpwgts[4] = .05;
  tpwgts[5] = .25;


  /*===========================================================================*/
  printf("\nTesting METIS_PartGraphRecursive ------------------------------------\n  ");
  tstnum = 1;

/**/
  numflag = 0; wgtflag = 0; nparts = 20;
  options[0] = 0;
  METIS_PartGraphRecursive(&nvtxs, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
                           &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, NULL, NULL, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 20;
  options[0] = 0;
  METIS_PartGraphRecursive(&nvtxs, xadj, adjncy, NULL, adjwgt, &wgtflag, &numflag,
                           &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, NULL, adjwgt, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 2; nparts = 20;
  options[0] = 0;
  METIS_PartGraphRecursive(&nvtxs, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag,
                           &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, vwgt, NULL, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 0;
  METIS_PartGraphRecursive(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                           &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_PartGraphRecursive(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                           &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_PartGraphRecursive(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                           &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n");


 
  /*===========================================================================*/
  printf("\nTesting METIS_WPartGraphRecursive -----------------------------------\n  ");
  tstnum = 1;


/**/
  numflag = 0; wgtflag = 0; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphRecursive(&nvtxs, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
                           &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, NULL, NULL, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphRecursive(&nvtxs, xadj, adjncy, NULL, adjwgt, &wgtflag, &numflag,
                           &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, NULL, adjwgt, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 2; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphRecursive(&nvtxs, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag,
                           &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, vwgt, NULL, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphRecursive(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                           &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_WPartGraphRecursive(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                           &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_WPartGraphRecursive(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                           &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n");



  /*===========================================================================*/
  printf("\nTesting METIS_PartGraphKway -----------------------------------------\n  ");
  tstnum = 1;

/**/
  numflag = 0; wgtflag = 0; nparts = 20;
  options[0] = 0;
  METIS_PartGraphKway(&nvtxs, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
                      &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, NULL, NULL, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 20;
  options[0] = 0;
  METIS_PartGraphKway(&nvtxs, xadj, adjncy, NULL, adjwgt, &wgtflag, &numflag,
                      &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, NULL, adjwgt, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 2; nparts = 20;
  options[0] = 0;
  METIS_PartGraphKway(&nvtxs, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag,
                      &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, vwgt, NULL, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 0;
  METIS_PartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                      &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_PartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                      &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_PartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                      &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 2; options[4] = 0;
  METIS_PartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                      &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 3; options[4] = 0;
  METIS_PartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                      &nparts, options, &edgecut, part);

  if ((rcode = VerifyPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n");


 
  /*===========================================================================*/
  printf("\nTesting METIS_WPartGraphKway ----------------------------------------\n  ");
  tstnum = 1;


/**/
  numflag = 0; wgtflag = 0; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphKway(&nvtxs, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
                       &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, NULL, NULL, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphKway(&nvtxs, xadj, adjncy, NULL, adjwgt, &wgtflag, &numflag,
                       &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, NULL, adjwgt, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 2; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphKway(&nvtxs, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag,
                       &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, vwgt, NULL, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                       &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_WPartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                       &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_WPartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                       &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 2; options[4] = 0;
  METIS_WPartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                       &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 3; options[4] = 0;
  METIS_WPartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                       &nparts, tpwgts, options, &edgecut, part);

  if ((rcode = VerifyWPart(nvtxs, xadj, adjncy, vwgt, adjwgt, nparts, tpwgts, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n");

  GKfree(&vwgt, &adjwgt, &part, LTERM);
}



/*************************************************************************
* This function verifies that the partitioning was computed correctly
**************************************************************************/
int VerifyPart(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
               idxtype *adjwgt, int nparts, int edgecut, idxtype *part)
{
  int i, j, k, cut, vfree=0, efree=0, rcode=0;
  idxtype *pwgts;

  if (part[idxamax(nvtxs, part)] != nparts-1)
    return 1;  /* the total number of partitions is different than nparts */

  /* compute the cut and the pwgts */
  if (vwgt == NULL) {
    vfree = 1;
    vwgt = idxsmalloc(nvtxs, 1, "vwgt");
  }
  if (adjwgt == NULL) {
    efree = 1;
    adjwgt = idxsmalloc(xadj[nvtxs], 1, "adjwgt");
  }
  pwgts = idxsmalloc(nparts, 0, "pwgts");

  for (cut=0, i=0; i<nvtxs; i++) {
    pwgts[part[i]] += vwgt[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) 
      cut += (part[i] != part[adjncy[j]] ? adjwgt[j] : 0);
  }

  if (cut != 2*edgecut)
    rcode = 2;

  if (nparts*pwgts[idxamax(nparts, pwgts)] > 1.10*idxsum(nparts, pwgts))
    rcode = 3;

  if (vfree)
    free(vwgt);

  if (efree)
    free(adjwgt);

  free(pwgts);

  MALLOC_CHECK(NULL);

  return rcode;
}


/*************************************************************************
* This function verifies that the partitioning was computed correctly
**************************************************************************/
int VerifyWPart(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                idxtype *adjwgt, int nparts, float *tpwgts, int edgecut, idxtype *part)
{
  int i, j, k, tvwgt, cut, vfree=0, efree=0, rcode=0;
  idxtype *pwgts;

  if (part[idxamax(nvtxs, part)] != nparts-1) 
    return 1;  /* the total number of partitions is different than nparts */

  /* compute the cut and the pwgts */
  if (vwgt == NULL) {
    vfree = 1;
    vwgt = idxsmalloc(nvtxs, 1, "vwgt");
  }
  if (adjwgt == NULL) {
    efree = 1;
    adjwgt = idxsmalloc(xadj[nvtxs], 1, "adjwgt");
  }
  pwgts = idxsmalloc(nparts, 0, "pwgts");

  for (cut=0, i=0; i<nvtxs; i++) {
    pwgts[part[i]] += vwgt[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) 
      cut += (part[i] != part[adjncy[j]] ? adjwgt[j] : 0);
  }

  if (cut != 2*edgecut)
    rcode = 2;

  tvwgt = idxsum(nparts, pwgts);
  for (i=0; i<nparts; i++) {
    if (pwgts[i] > 1.10*tpwgts[i]*tvwgt) {
      rcode = 3;
      break;
    }
  }

  if (vfree)
    free(vwgt);

  if (efree)
    free(adjwgt);

  free(pwgts);

  MALLOC_CHECK(NULL);

  return rcode;
}



/*************************************************************************
* This function tests the regular graph partitioning routines
**************************************************************************/
void Test_PartGraphV(int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, j, jj, k, tstnum, rcode;
  int nparts, numflag, wgtflag, totalv, options[10];
  idxtype *part, *vwgt, *vsize;
  float tpwgts[256];

  vwgt = idxmalloc(nvtxs, "vwgt");
  for (i=0; i<nvtxs; i++)
    vwgt[i] = RandomInRange(10);

  vsize = idxmalloc(nvtxs, "vsize");
  for (i=0; i<nvtxs; i++)
    vsize[i] = RandomInRange(10);

  part = idxmalloc(nvtxs, "part");

  tpwgts[0] = .1;
  tpwgts[1] = .2;
  tpwgts[2] = .3;
  tpwgts[3] = .1;
  tpwgts[4] = .05;
  tpwgts[5] = .25;


  /*===========================================================================*/
  printf("\nTesting METIS_PartGraphVKway ----------------------------------------\n  ");
  tstnum = 1;

/**/
  numflag = 0; wgtflag = 0; nparts = 20;
  options[0] = 0;
  METIS_PartGraphVKway(&nvtxs, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
                       &nparts, options, &totalv, part);

  if ((rcode = VerifyPartV(nvtxs, xadj, adjncy, NULL, NULL, nparts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 20;
  options[0] = 0;
  METIS_PartGraphVKway(&nvtxs, xadj, adjncy, NULL, vsize, &wgtflag, &numflag,
                       &nparts, options, &totalv, part);

  if ((rcode = VerifyPartV(nvtxs, xadj, adjncy, NULL, vsize, nparts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 2; nparts = 20;
  options[0] = 0;
  METIS_PartGraphVKway(&nvtxs, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag,
                       &nparts, options, &totalv, part);

  if ((rcode = VerifyPartV(nvtxs, xadj, adjncy, vwgt, NULL, nparts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 0;
  METIS_PartGraphVKway(&nvtxs, xadj, adjncy, vwgt, vsize, &wgtflag, &numflag,
                       &nparts, options, &totalv, part);

  if ((rcode = VerifyPartV(nvtxs, xadj, adjncy, vwgt, vsize, nparts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_PartGraphVKway(&nvtxs, xadj, adjncy, vwgt, vsize, &wgtflag, &numflag,
                       &nparts, options, &totalv, part);

  if ((rcode = VerifyPartV(nvtxs, xadj, adjncy, vwgt, vsize, nparts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_PartGraphVKway(&nvtxs, xadj, adjncy, vwgt, vsize, &wgtflag, &numflag,
                       &nparts, options, &totalv, part);

  if ((rcode = VerifyPartV(nvtxs, xadj, adjncy, vwgt, vsize, nparts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 20;
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 3; options[4] = 0;
  METIS_PartGraphVKway(&nvtxs, xadj, adjncy, vwgt, vsize, &wgtflag, &numflag,
                       &nparts, options, &totalv, part);

  if ((rcode = VerifyPartV(nvtxs, xadj, adjncy, vwgt, vsize, nparts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n");


 
  /*===========================================================================*/
  printf("\nTesting METIS_WPartGraphVKway ---------------------------------------\n  ");
  tstnum = 1;


/**/
  numflag = 0; wgtflag = 0; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphVKway(&nvtxs, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
                        &nparts, tpwgts, options, &totalv, part);

  if ((rcode = VerifyWPartV(nvtxs, xadj, adjncy, NULL, NULL, nparts, tpwgts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphVKway(&nvtxs, xadj, adjncy, NULL, vsize, &wgtflag, &numflag,
                        &nparts, tpwgts, options, &totalv, part);

  if ((rcode = VerifyWPartV(nvtxs, xadj, adjncy, NULL, vsize, nparts, tpwgts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 2; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphVKway(&nvtxs, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag,
                        &nparts, tpwgts, options, &totalv, part);

  if ((rcode = VerifyWPartV(nvtxs, xadj, adjncy, vwgt, NULL, nparts, tpwgts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 0;
  METIS_WPartGraphVKway(&nvtxs, xadj, adjncy, vwgt, vsize, &wgtflag, &numflag,
                        &nparts, tpwgts, options, &totalv, part);

  if ((rcode = VerifyWPartV(nvtxs, xadj, adjncy, vwgt, vsize, nparts, tpwgts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_WPartGraphVKway(&nvtxs, xadj, adjncy, vwgt, vsize, &wgtflag, &numflag,
                        &nparts, tpwgts, options, &totalv, part);

  if ((rcode = VerifyWPartV(nvtxs, xadj, adjncy, vwgt, vsize, nparts, tpwgts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_WPartGraphVKway(&nvtxs, xadj, adjncy, vwgt, vsize, &wgtflag, &numflag,
                        &nparts, tpwgts, options, &totalv, part);

  if ((rcode = VerifyWPartV(nvtxs, xadj, adjncy, vwgt, vsize, nparts, tpwgts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 3; nparts = 6;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 3; options[4] = 0;
  METIS_WPartGraphVKway(&nvtxs, xadj, adjncy, vwgt, vsize, &wgtflag, &numflag,
                        &nparts, tpwgts, options, &totalv, part);

  if ((rcode = VerifyWPartV(nvtxs, xadj, adjncy, vwgt, vsize, nparts, tpwgts, totalv, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n");


  GKfree(&vwgt, &vsize, &part, LTERM);
}


/*************************************************************************
* This function verifies that the partitioning was computed correctly
**************************************************************************/
int VerifyPartV(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
               idxtype *vsize, int nparts, int totalv, idxtype *part)
{
  int i, j, k, ttlv, vfree=0, efree=0, rcode=0;
  idxtype *pwgts, *marker;

  if (part[idxamax(nvtxs, part)] != nparts-1)
    return 1;  /* the total number of partitions is different than nparts */

  /* compute the cut and the pwgts */
  if (vwgt == NULL) {
    vfree = 1;
    vwgt = idxsmalloc(nvtxs, 1, "vwgt");
  }
  if (vsize == NULL) {
    efree = 1;
    vsize = idxsmalloc(nvtxs, 1, "vsize");
  }
  pwgts = idxsmalloc(nparts, 0, "pwgts");

  marker = idxsmalloc(nparts, -1, "htable");
  for (ttlv=0, i=0; i<nvtxs; i++) {
    pwgts[part[i]] += vwgt[i];
    marker[part[i]] = i;
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (marker[part[adjncy[j]]] != i) {
        ttlv += vsize[i];
        marker[part[adjncy[j]]] = i;
      }
    }
  }

  if (ttlv != totalv)
    rcode = 2;

  if (nparts*pwgts[idxamax(nparts, pwgts)] > 1.05*idxsum(nparts, pwgts))
    rcode = 3;

  if (vfree)
    free(vwgt);

  if (efree)
    free(vsize);

  free(pwgts);
  free(marker);

  MALLOC_CHECK(NULL);

  return rcode;
}


/*************************************************************************
* This function verifies that the partitioning was computed correctly
**************************************************************************/
int VerifyWPartV(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                idxtype *vsize, int nparts, float *tpwgts, int totalv, idxtype *part)
{
  int i, j, k, tvwgt, ttlv, vfree=0, efree=0, rcode=0;
  idxtype *pwgts, *marker;

  if (part[idxamax(nvtxs, part)] != nparts-1) 
    return 1;  /* the total number of partitions is different than nparts */

  /* compute the cut and the pwgts */
  if (vwgt == NULL) {
    vfree = 1;
    vwgt = idxsmalloc(nvtxs, 1, "vwgt");
  }
  if (vsize == NULL) {
    efree = 1;
    vsize = idxsmalloc(nvtxs, 1, "vsize");
  }
  pwgts = idxsmalloc(nparts, 0, "pwgts");

  marker = idxsmalloc(nparts, -1, "htable");
  for (ttlv=0, i=0; i<nvtxs; i++) {
    pwgts[part[i]] += vwgt[i];
    marker[part[i]] = i;
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (marker[part[adjncy[j]]] != i) {
        ttlv += vsize[i];
        marker[part[adjncy[j]]] = i;
      }
    }
  }

  if (ttlv != totalv)
    rcode = 2;

  tvwgt = idxsum(nparts, pwgts);
  for (i=0; i<nparts; i++) {
    if (pwgts[i] > 1.05*tpwgts[i]*tvwgt) {
      rcode = 3;
      break;
    }
  }

  if (vfree)
    free(vwgt);

  if (efree)
    free(vsize);

  free(pwgts);
  free(marker);

  MALLOC_CHECK(NULL);

  return rcode;
}



/*************************************************************************
* This function tests the regular graph partitioning routines
**************************************************************************/
void Test_PartGraphmC(int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, j, jj, k, tstnum, rcode;
  int nparts, ncon, numflag, wgtflag, edgecut, options[10];
  idxtype *part, *vwgt, *adjwgt;
  float ubvec[MAXNCON];

  ncon = 3;

  vwgt = idxmalloc(nvtxs*ncon, "vwgt");
  for (i=0; i<ncon*nvtxs; i++)
    vwgt[i] = RandomInRange(10);

  adjwgt = idxmalloc(xadj[nvtxs], "adjwgt");
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (i < k) {
        adjwgt[j] = 1+RandomInRange(5);
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          if (adjncy[jj] == i) {
            adjwgt[jj] = adjwgt[j];
            break;
          }
        }
      }
    }
  }


  part = idxmalloc(nvtxs, "part");



  /*===========================================================================*/
  printf("\nTesting METIS_mCPartGraphRecursive ----------------------------------\n  ");
  tstnum = 1;

/**/
  numflag = 0; wgtflag = 0; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, NULL, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 4; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 5; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 6; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 7; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 8; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n  ");

/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 4; options[2] = 2; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 3; options[2] = 2; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                             &nparts, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n");



  printf("\nTesting METIS_mCPartGraphKway ---------------------------------------\n  ");
  tstnum = 1;

/**/
  numflag = 0; wgtflag = 0; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, NULL, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 4; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 5; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 6; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 7; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 8; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n  ");

/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 6; options[2] = 2; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.5; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 4; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 1.05; ubvec[1] = 1.5; ubvec[2] = 1.5;
  options[0] = 1; options[1] = 3; options[2] = 2; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; wgtflag = 1; nparts = 10;
  ubvec[0] = 2.05; ubvec[1] = 1.05; ubvec[2] = 1.05;
  options[0] = 1; options[1] = 4; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_mCPartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag,
                        &nparts, ubvec, options, &edgecut, part);

  if ((rcode = VerifyPartmC(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, nparts, ubvec, edgecut, part)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


  printf("\n");

  GKfree(&vwgt, &adjwgt, &part, LTERM);
}



/*************************************************************************
* This function verifies that the partitioning was computed correctly
**************************************************************************/
int VerifyPartmC(int nvtxs, int ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
               idxtype *adjwgt, int nparts, float *ubvec, int edgecut, idxtype *part)
{
  int i, j, k, cut, efree=0, rcode=0;
  idxtype *pwgts;
  float lb;

  if (part[idxamax(nvtxs, part)] != nparts-1)
    return 1;  /* the total number of partitions is different than nparts */

  if (adjwgt == NULL) {
    efree = 1;
    adjwgt = idxsmalloc(xadj[nvtxs], 1, "adjwgt");
  }
  pwgts = idxsmalloc(ncon*nparts, 0, "pwgts");

  for (cut=0, i=0; i<nvtxs; i++) {
    for (j=0; j<ncon; j++)
      pwgts[part[i]*ncon+j] += vwgt[i*ncon+j];
    for (j=xadj[i]; j<xadj[i+1]; j++) 
      cut += (part[i] != part[adjncy[j]] ? adjwgt[j] : 0);
  }

  if (cut != 2*edgecut)
    rcode = 2;

/*
printf("\n");
for (i=0; i<nparts; i++) {
  for (j=0; j<ncon; j++) 
    printf("%5d ", pwgts[i*ncon+j]);
  printf("\n");
}
printf("---------------------------------\n");
for (j=0; j<ncon; j++) 
  printf("%5d ", idxsum_strd(nparts, pwgts+j, ncon));
printf("\n---------------------------------\n");
for (j=0; j<ncon; j++) 
  printf("%5d ", pwgts[ncon*idxamax_strd(nparts, pwgts+j, ncon)+j]);
printf("\n%d %d\n", idxsum(ncon*nvtxs, vwgt), idxsum(ncon*nparts, pwgts));
*/

  for (i=0; i<ncon; i++) {
    lb = 1.0*nparts*pwgts[ncon*idxamax_strd(nparts, pwgts+i, ncon)+i]/(1.0*idxsum_strd(nparts, pwgts+i, ncon));
    /*printf("[%3.2f]", lb);*/
    if (lb > ubvec[i])
      rcode = 3;
  }


  if (efree)
    free(adjwgt);

  free(pwgts);

  MALLOC_CHECK(NULL);

  return rcode;
}


/*************************************************************************
* This function tests the regular graph partitioning routines
**************************************************************************/
void Test_ND(int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, j, jj, k, tstnum, rcode;
  int numflag, wgtflag, options[10];
  idxtype *perm, *iperm, *vwgt;

  vwgt = idxmalloc(nvtxs, "vwgt");
  for (i=0; i<nvtxs; i++)
    vwgt[i] = 1+RandomInRange(10);


  perm = idxmalloc(nvtxs, "perm");
  iperm = idxmalloc(nvtxs, "iperm");



  /*===========================================================================*/
  printf("\nTesting METIS_EdgeND ------------------------------------------------\n  ");
  tstnum = 1;

/**/
  numflag = 0; 
  options[0] = 0;
  METIS_EdgeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_EdgeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_EdgeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n");



  /*===========================================================================*/
  printf("\nTesting METIS_NodeND ------------------------------------------------\n  ");
  tstnum = 1;

/**/
  numflag = 0; 
  options[0] = 0;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  options[5] = 0; options[6] = 0; options[7] = 1;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0; 
  options[5] = 0; options[6] = 0; options[7] = 1;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 1; options[4] = 0;
  options[5] = 0; options[6] = 0; options[7] = 1;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 2; options[3] = 1; options[4] = 0;
  options[5] = 0; options[6] = 0; options[7] = 1;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 2; options[4] = 0;
  options[5] = 0; options[6] = 0; options[7] = 1;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 2; options[4] = 0;
  options[5] = 1; options[6] = 0; options[7] = 1;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 2; options[4] = 0;
  options[5] = 2; options[6] = 0; options[7] = 1;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 2; options[4] = 0;
  options[5] = 3; options[6] = 0; options[7] = 1;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 2; options[4] = 0;
  options[5] = 3; options[6] = 40; options[7] = 1;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n  ");

/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 2; options[4] = 0;
  options[5] = 3; options[6] = 20; options[7] = 1;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 2; options[4] = 0;
  options[5] = 3; options[6] = 20; options[7] = 2;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 2; options[4] = 0;
  options[5] = 0; options[6] = 0; options[7] = 2;
  METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n");



  /*===========================================================================*/
  printf("\nTesting METIS_NodeWND -----------------------------------------------\n  ");
  tstnum = 1;

/**/
  numflag = 0; 
  options[0] = 0;
  METIS_NodeWND(&nvtxs, xadj, adjncy, vwgt, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 1; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_NodeWND(&nvtxs, xadj, adjncy, vwgt, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 2; options[2] = 1; options[3] = 1; options[4] = 0; 
  METIS_NodeWND(&nvtxs, xadj, adjncy, vwgt, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 1; options[4] = 0;
  METIS_NodeWND(&nvtxs, xadj, adjncy, vwgt, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 2; options[3] = 1; options[4] = 0;
  METIS_NodeWND(&nvtxs, xadj, adjncy, vwgt, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);


/**/
  numflag = 0; 
  options[0] = 1; options[1] = 3; options[2] = 1; options[3] = 2; options[4] = 0;
  METIS_NodeWND(&nvtxs, xadj, adjncy, vwgt, &numflag, options, perm, iperm);

  if ((rcode = VerifyND(nvtxs, perm, iperm)) == 0)
    printf("[%d:ok]", tstnum++);
  else
    printf("[%d:err-%d]", tstnum++, rcode);
  fflush(stdout);

  printf("\n");


  GKfree(&vwgt, &perm, &iperm, LTERM);
}



/*************************************************************************
* This function verifies that the partitioning was computed correctly
**************************************************************************/
int VerifyND(int nvtxs, idxtype *perm, idxtype *iperm)
{
  int i, j, k, rcode=0;

  for (i=0; i<nvtxs; i++) {
    if (i != perm[iperm[i]])
      rcode = 1;
  }

  for (i=0; i<nvtxs; i++) {
    if (i != iperm[perm[i]])
      rcode = 2;
  }

  MALLOC_CHECK(NULL);

  return rcode;
}


