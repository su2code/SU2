/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * smbfactor.c
 *
 * This file performs the symbolic factorization of a matrix
 *
 * Started 8/1/97
 * George
 *
 * $Id: smbfactor.c,v 1.1 1998/11/27 17:59:40 karypis Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function sets up data structures for fill-in computations
**************************************************************************/
void ComputeFillIn(GraphType *graph, idxtype *iperm)
{
  int i, j, k, nvtxs, maxlnz, maxsub;
  idxtype *xadj, *adjncy;
  idxtype *perm, *xlnz, *xnzsub, *nzsub;
  double opc;

/*
  printf("\nSymbolic factorization... --------------------------------------------\n");
*/

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  maxsub = 4*xadj[nvtxs];

  /* Relabel the vertices so that it starts from 1 */
  k = xadj[nvtxs];
  for (i=0; i<k; i++)
    adjncy[i]++;
  for (i=0; i<nvtxs+1; i++)
    xadj[i]++;

  /* Allocate the required memory */
  perm = idxmalloc(nvtxs+1, "ComputeFillIn: perm");
  xlnz = idxmalloc(nvtxs+1, "ComputeFillIn: xlnz");
  xnzsub = idxmalloc(nvtxs+1, "ComputeFillIn: xnzsub");
  nzsub = idxmalloc(maxsub, "ComputeFillIn: nzsub");

  /* Construct perm from iperm and change the numbering of iperm */
  for (i=0; i<nvtxs; i++)
    perm[iperm[i]] = i;
  for (i=0; i<nvtxs; i++) {
    iperm[i]++;
    perm[i]++;
  }
  
  /*
   * Call sparspak routine.
   */
  if (smbfct(nvtxs, xadj, adjncy, perm, iperm, xlnz, &maxlnz, xnzsub, nzsub, &maxsub)) {
    free(nzsub);

    maxsub = 4*maxsub; 
    nzsub = idxmalloc(maxsub, "ComputeFillIn: nzsub");
    if (smbfct(nvtxs, xadj, adjncy, perm, iperm, xlnz, &maxlnz, xnzsub, nzsub, &maxsub)) 
      errexit("MAXSUB is too small!");
  }

  opc = 0;
  for (i=0; i<nvtxs; i++)
    xlnz[i]--;
  for (i=0; i<nvtxs; i++)
    opc += (xlnz[i+1]-xlnz[i])*(xlnz[i+1]-xlnz[i]) - (xlnz[i+1]-xlnz[i]);

  printf("  Nonzeros: %d, \tOperation Count: %6.4le\n", maxlnz, opc);


  GKfree(&perm, &xlnz, &xnzsub, &nzsub, LTERM);


  /* Relabel the vertices so that it starts from 0 */
  for (i=0; i<nvtxs; i++)
    iperm[i]--;
  for (i=0; i<nvtxs+1; i++)
    xadj[i]--;
  k = xadj[nvtxs];
  for (i=0; i<k; i++)
    adjncy[i]--;

}



/*************************************************************************
* This function sets up data structures for fill-in computations
**************************************************************************/
idxtype ComputeFillIn2(GraphType *graph, idxtype *iperm)
{
  int i, j, k, nvtxs, maxlnz, maxsub;
  idxtype *xadj, *adjncy;
  idxtype *perm, *xlnz, *xnzsub, *nzsub;
  double opc;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  maxsub = 4*xadj[nvtxs];

  /* Relabel the vertices so that it starts from 1 */
  k = xadj[nvtxs];
  for (i=0; i<k; i++)
    adjncy[i]++;
  for (i=0; i<nvtxs+1; i++)
    xadj[i]++;

  /* Allocate the required memory */
  perm = idxmalloc(nvtxs+1, "ComputeFillIn: perm");
  xlnz = idxmalloc(nvtxs+1, "ComputeFillIn: xlnz");
  xnzsub = idxmalloc(nvtxs+1, "ComputeFillIn: xnzsub");
  nzsub = idxmalloc(maxsub, "ComputeFillIn: nzsub");

  /* Construct perm from iperm and change the numbering of iperm */
  for (i=0; i<nvtxs; i++)
    perm[iperm[i]] = i;
  for (i=0; i<nvtxs; i++) {
    iperm[i]++;
    perm[i]++;
  }
  
  /*
   * Call sparspak routine.
   */
  if (smbfct(nvtxs, xadj, adjncy, perm, iperm, xlnz, &maxlnz, xnzsub, nzsub, &maxsub)) {
    free(nzsub);

    maxsub = 4*maxsub; 
    nzsub = idxmalloc(maxsub, "ComputeFillIn: nzsub");
    if (smbfct(nvtxs, xadj, adjncy, perm, iperm, xlnz, &maxlnz, xnzsub, nzsub, &maxsub)) 
      errexit("MAXSUB is too small!");
  }

  opc = 0;
  for (i=0; i<nvtxs; i++)
    xlnz[i]--;
  for (i=0; i<nvtxs; i++)
    opc += (xlnz[i+1]-xlnz[i])*(xlnz[i+1]-xlnz[i]) - (xlnz[i+1]-xlnz[i]);


  GKfree(&perm, &xlnz, &xnzsub, &nzsub, LTERM);


  /* Relabel the vertices so that it starts from 0 */
  for (i=0; i<nvtxs; i++)
    iperm[i]--;
  for (i=0; i<nvtxs+1; i++)
    xadj[i]--;
  k = xadj[nvtxs];
  for (i=0; i<k; i++)
    adjncy[i]--;

  return maxlnz;

}


/*****************************************************************          
**********     SMBFCT ..... SYMBOLIC FACTORIZATION       ********* 
******************************************************************          
*   PURPOSE - THIS ROUTINE PERFORMS SYMBOLIC FACTORIZATION               
*   ON A PERMUTED LINEAR SYSTEM AND IT ALSO SETS UP THE               
*   COMPRESSED DATA STRUCTURE FOR THE SYSTEM.                         
*
*   INPUT PARAMETERS -                                               
*      NEQNS - NUMBER OF EQUATIONS.                                 
*      (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.                   
*      (PERM, INVP) - THE PERMUTATION VECTOR AND ITS INVERSE.     
*
*   UPDATED PARAMETERS -                                         
*      MAXSUB - SIZE OF THE SUBSCRIPT ARRAY NZSUB.  ON RETURN,  
*             IT CONTAINS THE NUMBER OF SUBSCRIPTS USED        
*
*   OUTPUT PARAMETERS -                                       
*      XLNZ - INDEX INTO THE NONZERO STORAGE VECTOR LNZ.   
*      (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT VECTORS. 
*      MAXLNZ - THE NUMBER OF NONZEROS FOUND.             
*
*******************************************************************/
int smbfct(int neqns, idxtype *xadj, idxtype *adjncy, idxtype *perm, idxtype *invp, idxtype *xlnz, 
           int *maxlnz, idxtype *xnzsub, idxtype *nzsub, int *maxsub)
{
  /* Local variables */
  int node, rchm, mrgk, lmax, i, j, k, m, nabor, nzbeg, nzend;
  int kxsub, jstop, jstrt, mrkflg, inz, knz, flag;
  idxtype *mrglnk, *marker, *rchlnk;

  rchlnk = idxmalloc(neqns+1, "smbfct: rchlnk");
  marker = idxsmalloc(neqns+1, 0, "smbfct: marker");
  mrglnk = idxsmalloc(neqns+1, 0, "smbfct: mgrlnk");

  /* Parameter adjustments */
  --marker;
  --mrglnk;
  --rchlnk;
  --nzsub;
  --xnzsub;
  --xlnz;
  --invp;
  --perm;
  --adjncy;
  --xadj;

  /* Function Body */
  flag = 0;
  nzbeg = 1;
  nzend = 0;
  xlnz[1] = 1;

  /* FOR EACH COLUMN KNZ COUNTS THE NUMBER OF NONZEROS IN COLUMN K ACCUMULATED IN RCHLNK. */
  for (k = 1; k <= neqns; ++k) {
    knz = 0;
    mrgk = mrglnk[k];
    mrkflg = 0;
    marker[k] = k;
    if (mrgk != 0) 
      marker[k] = marker[mrgk];
    xnzsub[k] = nzend;
    node = perm[k];

    if (xadj[node] >= xadj[node+1]) {
      xlnz[k+1] = xlnz[k];
      continue;
    }

    /* USE RCHLNK TO LINK THROUGH THE STRUCTURE OF A(*,K) BELOW DIAGONAL */
    rchlnk[k] = neqns+1;
    for (j=xadj[node]; j<xadj[node+1]; j++) {
      nabor = invp[adjncy[j]];
      if (nabor <= k) 
        continue;
      rchm = k;

      do {
        m = rchm;
        rchm = rchlnk[m];
      } while (rchm <= nabor); 

      knz++;
      rchlnk[m] = nabor;
      rchlnk[nabor] = rchm;
      if (marker[nabor] != marker[k]) 
        mrkflg = 1;
    }

    /* TEST FOR MASS SYMBOLIC ELIMINATION */
    lmax = 0;
    if (mrkflg != 0 || mrgk == 0 || mrglnk[mrgk] != 0) 
      goto L350;
    xnzsub[k] = xnzsub[mrgk] + 1;
    knz = xlnz[mrgk + 1] - (xlnz[mrgk] + 1);
    goto L1400;


    /* LINK THROUGH EACH COLUMN I THAT AFFECTS L(*,K) */
L350:
    i = k;
    while ((i = mrglnk[i]) != 0) {
      inz = xlnz[i+1] - (xlnz[i]+1);
      jstrt = xnzsub[i] + 1;
      jstop = xnzsub[i] + inz;

      if (inz > lmax) { 
        lmax = inz;
        xnzsub[k] = jstrt;
      }

      /* MERGE STRUCTURE OF L(*,I) IN NZSUB INTO RCHLNK. */ 
      rchm = k;
      for (j = jstrt; j <= jstop; ++j) {
        nabor = nzsub[j];
        do {
          m = rchm;
          rchm = rchlnk[m];
        } while (rchm < nabor);

        if (rchm != nabor) {
          knz++;
          rchlnk[m] = nabor;
          rchlnk[nabor] = rchm;
          rchm = nabor;
        }
      }
    }

    /* CHECK IF SUBSCRIPTS DUPLICATE THOSE OF ANOTHER COLUMN */
    if (knz == lmax) 
      goto L1400;

    /* OR IF TAIL OF K-1ST COLUMN MATCHES HEAD OF KTH */
    if (nzbeg > nzend) 
      goto L1200;

    i = rchlnk[k];
    for (jstrt = nzbeg; jstrt <= nzend; ++jstrt) {
      if (nzsub[jstrt] < i) 
        continue;

      if (nzsub[jstrt] == i) 
        goto L1000;
      else 
        goto L1200;
    }
    goto L1200;

L1000:
    xnzsub[k] = jstrt;
    for (j = jstrt; j <= nzend; ++j) {
      if (nzsub[j] != i) 
        goto L1200;
      
      i = rchlnk[i];
      if (i > neqns) 
        goto L1400;
    }
    nzend = jstrt - 1;

    /* COPY THE STRUCTURE OF L(*,K) FROM RCHLNK TO THE DATA STRUCTURE (XNZSUB, NZSUB) */
L1200:
    nzbeg = nzend + 1;
    nzend += knz;

    if (nzend > *maxsub) {
      flag = 1; /* Out of memory */
      break;
    }

    i = k;
    for (j=nzbeg; j<=nzend; ++j) {
      i = rchlnk[i];
      nzsub[j] = i;
      marker[i] = k;
    }
    xnzsub[k] = nzbeg;
    marker[k] = k;

    /*
     * UPDATE THE VECTOR MRGLNK.  NOTE COLUMN L(*,K) JUST FOUND   
     * IS REQUIRED TO DETERMINE COLUMN L(*,J), WHERE              
     * L(J,K) IS THE FIRST NONZERO IN L(*,K) BELOW DIAGONAL.      
     */
L1400:
    if (knz > 1) { 
      kxsub = xnzsub[k];
      i = nzsub[kxsub];
      mrglnk[k] = mrglnk[i];
      mrglnk[i] = k;
    }

    xlnz[k + 1] = xlnz[k] + knz;
  }

  if (flag == 0) {
    *maxlnz = xlnz[neqns] - 1;
    *maxsub = xnzsub[neqns];
    xnzsub[neqns + 1] = xnzsub[neqns];
  }

  marker++;
  mrglnk++;
  rchlnk++;
  nzsub++;
  xnzsub++;
  xlnz++;
  invp++;
  perm++;
  adjncy++;
  xadj++;
  GKfree(&rchlnk, &mrglnk, &marker, LTERM);

  return flag;
  
} 

