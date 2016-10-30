/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mesh.c
 *
 * This file contains routines for constructing the dual graph of a mesh.
 * Assumes that each processor has at least one mesh element.
 *
 * Started 10/19/94
 * George
 *
 * $Id: mesh.c 10575 2011-07-14 14:46:42Z karypis $
 *
 */

#include <parmetislib.h>


/*************************************************************************
* This function converts a mesh into a dual graph
**************************************************************************/
int ParMETIS_V3_Mesh2Dual(idx_t *elmdist, idx_t *eptr, idx_t *eind, 
                 idx_t *numflag, idx_t *ncommon, idx_t **r_xadj, 
		 idx_t **r_adjncy, MPI_Comm *comm)
{
  idx_t i, j, jj, k, kk, m;
  idx_t npes, mype, pe, count, mask, pass;
  idx_t nelms, lnns, my_nns, node;
  idx_t firstelm, firstnode, lnode, nrecv, nsend;
  idx_t *scounts, *rcounts, *sdispl, *rdispl;
  idx_t *nodedist, *nmap, *auxarray;
  idx_t *gnptr, *gnind, *nptr, *nind, *myxadj=NULL, *myadjncy = NULL;
  idx_t *sbuffer, *rbuffer, *htable;
  ikv_t *nodelist, *recvbuffer;
  idx_t maxcount, *ind, *wgt;
  idx_t gmaxnode, gminnode;
  size_t curmem;

  gk_malloc_init();
  curmem = gk_GetCurMemoryUsed();
  
  /* Get basic comm info */
  gkMPI_Comm_size(*comm, &npes);
  gkMPI_Comm_rank(*comm, &mype);


  nelms = elmdist[mype+1]-elmdist[mype];

  if (*numflag > 0) 
    ChangeNumberingMesh(elmdist, eptr, eind, NULL, NULL, NULL, npes, mype, 1);

  mask = (1<<11)-1;

  /*****************************/
  /* Determine number of nodes */
  /*****************************/
  gminnode = GlobalSEMinComm(*comm, imin(eptr[nelms], eind));
  for (i=0; i<eptr[nelms]; i++)
    eind[i] -= gminnode;

  gmaxnode = GlobalSEMaxComm(*comm, imax(eptr[nelms], eind));


  /**************************/
  /* Check for input errors */
  /**************************/
  ASSERT(nelms > 0);

  /* construct node distribution array */
  nodedist = ismalloc(npes+1, 0, "nodedist");
  for (nodedist[0]=0, i=0,j=gmaxnode+1; i<npes; i++) {
    k = j/(npes-i);
    nodedist[i+1] = nodedist[i]+k;
    j -= k;
  }
  my_nns = nodedist[mype+1]-nodedist[mype];
  firstnode = nodedist[mype];

  nodelist = ikvmalloc(eptr[nelms], "nodelist");
  auxarray = imalloc(eptr[nelms], "auxarray");
  htable   = ismalloc(gk_max(my_nns, mask+1), -1, "htable");
  scounts  = imalloc(npes, "scounts");
  rcounts  = imalloc(npes, "rcounts");
  sdispl   = imalloc(npes+1, "sdispl");
  rdispl   = imalloc(npes+1, "rdispl");


  /*********************************************/
  /* first find a local numbering of the nodes */
  /*********************************************/
  for (i=0; i<nelms; i++) {
    for (j=eptr[i]; j<eptr[i+1]; j++) {
      nodelist[j].key = eind[j];
      nodelist[j].val = j;
      auxarray[j]     = i; /* remember the local element ID that uses this node */
    }
  }
  ikvsorti(eptr[nelms], nodelist);

  for (count=1, i=1; i<eptr[nelms]; i++) {
    if (nodelist[i].key > nodelist[i-1].key)
      count++;
  }

  lnns = count;
  nmap = imalloc(lnns, "nmap");

  /* renumber the nodes of the elements array */
  count = 1;
  nmap[0] = nodelist[0].key;
  eind[nodelist[0].val] = 0;
  nodelist[0].val = auxarray[nodelist[0].val];  /* Store the local element ID */
  for (i=1; i<eptr[nelms]; i++) {
    if (nodelist[i].key > nodelist[i-1].key) {
      nmap[count] = nodelist[i].key;
      count++;
    }
    eind[nodelist[i].val] = count-1;
    nodelist[i].val = auxarray[nodelist[i].val];  /* Store the local element ID */
  }
  gkMPI_Barrier(*comm);

  /**********************************************************/
  /* perform comms necessary to construct node-element list */
  /**********************************************************/
  iset(npes, 0, scounts);
  for (pe=i=0; i<eptr[nelms]; i++) {
    while (nodelist[i].key >= nodedist[pe+1])
      pe++;
    scounts[pe] += 2;
  }
  ASSERT(pe < npes);

  gkMPI_Alltoall((void *)scounts, 1, IDX_T, (void *)rcounts, 1, IDX_T, *comm);

  icopy(npes, scounts, sdispl);
  MAKECSR(i, npes, sdispl);

  icopy(npes, rcounts, rdispl);
  MAKECSR(i, npes, rdispl);

  ASSERT(sdispl[npes] == eptr[nelms]*2);

  nrecv = rdispl[npes]/2;
  recvbuffer = ikvmalloc(gk_max(1, nrecv), "recvbuffer");

  gkMPI_Alltoallv((void *)nodelist, scounts, sdispl, IDX_T, (void *)recvbuffer, 
      rcounts, rdispl, IDX_T, *comm);

  /**************************************/
  /* construct global node-element list */
  /**************************************/
  gnptr = ismalloc(my_nns+1, 0, "gnptr");

  for (i=0; i<npes; i++) {
    for (j=rdispl[i]/2; j<rdispl[i+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      ASSERT(lnode >= 0 && lnode < my_nns)

      gnptr[lnode]++;
    }
  }
  MAKECSR(i, my_nns, gnptr);

  gnind = imalloc(gk_max(1, gnptr[my_nns]), "gnind");
  for (pe=0; pe<npes; pe++) {
    firstelm = elmdist[pe];
    for (j=rdispl[pe]/2; j<rdispl[pe+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      gnind[gnptr[lnode]++] = recvbuffer[j].val+firstelm;
    }
  }
  SHIFTCSR(i, my_nns, gnptr);


  /*********************************************************/
  /* send the node-element info to the relevant processors */
  /*********************************************************/
  iset(npes, 0, scounts);

  /* use a hash table to ensure that each node is sent to a proc only once */
  for (pe=0; pe<npes; pe++) {
    for (j=rdispl[pe]/2; j<rdispl[pe+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      if (htable[lnode] == -1) {
        scounts[pe] += gnptr[lnode+1]-gnptr[lnode];
        htable[lnode] = 1;
      }
    }

    /* now reset the hash table */
    for (j=rdispl[pe]/2; j<rdispl[pe+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      htable[lnode] = -1;
    }
  }


  gkMPI_Alltoall((void *)scounts, 1, IDX_T, (void *)rcounts, 1, IDX_T, *comm);

  icopy(npes, scounts, sdispl);
  MAKECSR(i, npes, sdispl);

  /* create the send buffer */
  nsend = sdispl[npes];
  sbuffer = imalloc(gk_max(1, nsend), "sbuffer");

  count = 0;
  for (pe=0; pe<npes; pe++) {
    for (j=rdispl[pe]/2; j<rdispl[pe+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      if (htable[lnode] == -1) {
        for (k=gnptr[lnode]; k<gnptr[lnode+1]; k++) {
          if (k == gnptr[lnode])
            sbuffer[count++] = -1*(gnind[k]+1);
          else
            sbuffer[count++] = gnind[k];
        }
        htable[lnode] = 1;
      }
    }
    ASSERT(count == sdispl[pe+1]);

    /* now reset the hash table */
    for (j=rdispl[pe]/2; j<rdispl[pe+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      htable[lnode] = -1;
    }
  }

  icopy(npes, rcounts, rdispl);
  MAKECSR(i, npes, rdispl);

  nrecv   = rdispl[npes];
  rbuffer = imalloc(gk_max(1, nrecv), "rbuffer");

  gkMPI_Alltoallv((void *)sbuffer, scounts, sdispl, IDX_T, (void *)rbuffer, 
      rcounts, rdispl, IDX_T, *comm);

  k = -1;
  nptr = ismalloc(lnns+1, 0, "nptr");
  nind = rbuffer;
  for (pe=0; pe<npes; pe++) {
    for (j=rdispl[pe]; j<rdispl[pe+1]; j++) {
      if (nind[j] < 0) {
        k++;
        nind[j] = (-1*nind[j])-1;
      }
      nptr[k]++;
    }
  }
  MAKECSR(i, lnns, nptr);

  ASSERT(k+1 == lnns);
  ASSERT(nptr[lnns] == nrecv)

  myxadj = *r_xadj = (idx_t *)malloc(sizeof(idx_t)*(nelms+1));
  if (myxadj == NULL) 
    gk_errexit(SIGMEM, "Failed to allocate memory for the dual graph's xadj array.\n");
  iset(nelms+1, 0, myxadj);

  iset(mask+1, -1, htable);

  firstelm = elmdist[mype];

  /* Two passes -- in first pass, simply find out the memory requirements */
  maxcount = 200;
  ind = imalloc(maxcount, "ParMETIS_V3_Mesh2Dual: ind");
  wgt = imalloc(maxcount, "ParMETIS_V3_Mesh2Dual: wgt");

  for (pass=0; pass<2; pass++) {
    for (i=0; i<nelms; i++) {
      for (count=0, j=eptr[i]; j<eptr[i+1]; j++) {
        node = eind[j];

        for (k=nptr[node]; k<nptr[node+1]; k++) {
          if ((kk=nind[k]) == firstelm+i) 
	    continue;
	    
          m = htable[(kk&mask)];

          if (m == -1) {
            ind[count] = kk;
            wgt[count] = 1;
            htable[(kk&mask)] = count++;
          }
          else {
            if (ind[m] == kk) { 
              wgt[m]++;
            }
            else {
              for (jj=0; jj<count; jj++) {
                if (ind[jj] == kk) {
                  wgt[jj]++;
                  break;
	        }
              }
              if (jj == count) {
                ind[count]   = kk;
                wgt[count++] = 1;
              }
	    }
          }

          /* Adjust the memory. 
             This will be replaced by a idxrealloc() when GKlib will be incorporated */
          if (count == maxcount-1) {
            maxcount *= 2;
            ind = irealloc(ind, maxcount, "ParMETIS_V3_Mesh2Dual: ind");
            wgt = irealloc(wgt, maxcount, "ParMETIS_V3_Mesh2Dual: wgt");
          }
        }
      }

      for (j=0; j<count; j++) {
        htable[(ind[j]&mask)] = -1;
        if (wgt[j] >= *ncommon) {
          if (pass == 0) 
            myxadj[i]++;
          else 
            myadjncy[myxadj[i]++] = ind[j];
	}
      }
    }

    if (pass == 0) {
      MAKECSR(i, nelms, myxadj);
      myadjncy = *r_adjncy = (idx_t *)malloc(sizeof(idx_t)*myxadj[nelms]);
      if (myadjncy == NULL)
        gk_errexit(SIGMEM, "Failed to allocate memory for dual graph's adjncy array.\n");
    }
    else {
      SHIFTCSR(i, nelms, myxadj);
    }
  }

  /*****************************************/
  /* correctly renumber the elements array */
  /*****************************************/
  for (i=0; i<eptr[nelms]; i++)
    eind[i] = nmap[eind[i]] + gminnode;

  if (*numflag == 1) 
    ChangeNumberingMesh(elmdist, eptr, eind, myxadj, myadjncy, NULL, npes, mype, 0);

  /* do not free nodelist, recvbuffer, rbuffer */
  gk_free((void **)&nodedist, &nodelist, &auxarray, &htable, &scounts, &rcounts,
      &sdispl, &rdispl, &nmap, &recvbuffer, &gnptr, &gnind, &sbuffer, &rbuffer,
      &nptr, &ind, &wgt, LTERM);

  if (gk_GetCurMemoryUsed() - curmem > 0) {
    printf("ParMETIS appears to have a memory leak of %zdbytes. Report this.\n",
        (ssize_t)(gk_GetCurMemoryUsed() - curmem));
  }
  gk_malloc_cleanup(0);

  return METIS_OK;
}

