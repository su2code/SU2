/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * comm.c
 *
 * This function provides various high level communication functions 
 *
 * $Id: comm.c 10592 2011-07-16 21:17:53Z karypis $
 */

#include <parmetislib.h>


/*************************************************************************/
/*! This function performs the following functions:
    - determines the processors that contain adjacent vertices and setup
      the infrastructure for efficient communication.
    - localizes the numbering of the adjancency lists.
*/
/**************************************************************************/
void CommSetup(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, k, islocal, penum, gnvtxs, nvtxs, nlocal, firstvtx, lastvtx, 
        nsend, nrecv, nnbrs, nadj;
  idx_t npes=ctrl->npes, mype=ctrl->mype;
  idx_t *vtxdist, *xadj, *adjncy;
  idx_t *peind, *recvptr, *recvind, *sendptr, *sendind;
  idx_t *imap, *lperm;
  idx_t *pexadj, *peadjncy, *peadjloc, *startsind;
  ikv_t *recvrequests, *sendrequests, *adjpairs;

  if (graph->lperm != NULL)
    return; /* The communication structure has already been setup */

  STARTTIMER(ctrl, ctrl->SetupTmr);

  gnvtxs  = graph->gnvtxs;
  nvtxs   = graph->nvtxs;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  lperm   = graph->lperm = iincset(nvtxs, 0, imalloc(nvtxs, "CommSetup: graph->lperm"));

  vtxdist  = graph->vtxdist;
  firstvtx = vtxdist[mype];
  lastvtx  = vtxdist[mype+1];

  WCOREPUSH;
  /************************************************************* 
   * Determine what you need to receive 
   *************************************************************/
  /* first pass: determine nadj and interior/interface vertices */
  for (nlocal=0, nadj=0, i=0; i<nvtxs; i++) {
    islocal = 1;
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (k < firstvtx || k >= lastvtx) { /* remote vertex */
        nadj++;
        islocal = 0;
      }
    }
    if (islocal) {
      lperm[i]        = lperm[nlocal];
      lperm[nlocal++] = i;
    }
  }
  graph->nlocal = nlocal;

  adjpairs = ikvwspacemalloc(ctrl, nadj+1);

  /* second pass: rewrite locale entries and populate remote edges */
  for (nadj=0, i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (k >= firstvtx && k < lastvtx) { /* local vertex */
        adjncy[j] = k-firstvtx; 
      }
      else { /* remote vertex */
        adjpairs[nadj].key   = k;
        adjpairs[nadj++].val = j;
      }
    }
  }

  /* use a sort-based "unique" approach */
  ikvsorti(nadj, adjpairs);
  adjpairs[nadj].key = gnvtxs+1;  /* boundary condition */

  /* determine how many distinct vertices you need to receive */
  for (nrecv=0, i=0; i<nadj; i++) {
    if (adjpairs[i].key != adjpairs[i+1].key)
      nrecv++;
  }
  graph->nrecv = nrecv;


  /* allocate space for the to be received vertices part of the recvinfo */
  recvind = graph->recvind = imalloc(nrecv, "CommSetup: recvind");

  /* store distinct vertices into recvind array and re-write adjncy */
  for (nrecv=0, i=0; i<nadj; i++) {
    adjncy[adjpairs[i].val] = nvtxs+nrecv;
    if (adjpairs[i].key != adjpairs[i+1].key)
      recvind[nrecv++] = adjpairs[i].key;
  }
  PASSERT(ctrl, nrecv == graph->nrecv);


  /* determine the number of neighboring processors */
  for (i=0, nnbrs=0, penum=0; penum<npes; penum++) {
    for (j=i; j<nrecv; j++) {
      if (recvind[j] >= vtxdist[penum+1])
        break;
    }
    if (j > i) {
      nnbrs++;
      i = j;
    }
  }
  graph->nnbrs = nnbrs;

  /* Update the ctrl arrays that have to do with p2p communication */
  CommUpdateNnbrs(ctrl, nnbrs);

  /* allocate space for peind/recvptr part of the recvinfo */
  peind   = graph->peind   = imalloc(nnbrs, "CommSetup: peind");
  recvptr = graph->recvptr = imalloc(nnbrs+1, "CommSetup: recvptr");

  /* populate the peind/recvptr arrays */
  for (i=0, nnbrs=0, recvptr[0]=0, penum=0; penum<npes; penum++) {
    for (j=i; j<nrecv; j++) {
      if (recvind[j] >= vtxdist[penum+1])
        break;
    }
    if (j > i) {
      peind[nnbrs++] = penum;
      recvptr[nnbrs] = j;
      i = j;
    }
  }
  PASSERT(ctrl, nnbrs == graph->nnbrs);

  WCOREPOP;

  /* PrintVector(ctrl, nnbrs+1, 0, recvptr, "recvptr"); */


  WCOREPUSH;
  /************************************************************* 
   * Determine what you need to send 
   *************************************************************/
  /* GKTODO - This can be replaced via a sparse communication */
  /* Tell the other processors what they need to send you */
  recvrequests = ikvwspacemalloc(ctrl, npes);
  sendrequests = ikvwspacemalloc(ctrl, npes);
  memset(recvrequests, 0, sizeof(ikv_t)*npes);
  for (i=0; i<nnbrs; i++) {
    recvrequests[peind[i]].key = recvptr[i+1]-recvptr[i];
    recvrequests[peind[i]].val = nvtxs+recvptr[i];
  }
  gkMPI_Alltoall((void *)recvrequests, 2, IDX_T, (void *)sendrequests, 
      2, IDX_T, ctrl->comm);

  /* PrintPairs(ctrl, npes, recvrequests, "recvrequests"); */
  /* PrintPairs(ctrl, npes, sendrequests, "sendrequests"); */

  startsind = iwspacemalloc(ctrl, nnbrs);
  sendptr   = graph->sendptr = imalloc(nnbrs+1, "CommSetup: sendptr");

  for (j=0, i=0; i<npes; i++) {
    if (sendrequests[i].key > 0) {
      sendptr[j]   = sendrequests[i].key;
      startsind[j] = sendrequests[i].val;
      j++;
    }
  }
  PASSERT(ctrl, j == nnbrs);

  MAKECSR(i, nnbrs, sendptr);

  nsend   = graph->nsend   = sendptr[nnbrs];
  sendind = graph->sendind = imalloc(nsend, "CommSetup: sendind");


  /* Issue the receives for sendind */
  for (i=0; i<nnbrs; i++) {
    gkMPI_Irecv((void *)(sendind+sendptr[i]), sendptr[i+1]-sendptr[i], IDX_T, 
              peind[i], 1, ctrl->comm, ctrl->rreq+i);
  }

  /* Issue the sends. My recvind[penum] becomes penum's sendind[mype] */
  for (i=0; i<nnbrs; i++) {
    gkMPI_Isend((void *)(recvind+recvptr[i]), recvptr[i+1]-recvptr[i], IDX_T,
              peind[i], 1, ctrl->comm, ctrl->sreq+i);
  }

  gkMPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
  gkMPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


  /* Create the peadjncy data structure for sparse boundary exchanges */
  pexadj   = graph->pexadj   = ismalloc(nvtxs+1, 0, "CommSetup: pexadj");
  peadjncy = graph->peadjncy = imalloc(nsend, "CommSetup: peadjncy");
  peadjloc = graph->peadjloc = imalloc(nsend, "CommSetup: peadjloc");

  for (i=0; i<nsend; i++) {
    PASSERTP(ctrl, sendind[i] >= firstvtx && sendind[i] < lastvtx, 
            (ctrl, "%"PRIDX" %"PRIDX" %"PRIDX"\n", sendind[i], firstvtx, lastvtx));
    pexadj[sendind[i]-firstvtx]++;
  }
  MAKECSR(i, nvtxs, pexadj);

  for (i=0; i<nnbrs; i++) {
    for (j=sendptr[i]; j<sendptr[i+1]; j++) {
      k = pexadj[sendind[j]-firstvtx]++;
      peadjncy[k] = i;  /* peind[i] is the actual PE number */
      peadjloc[k] = startsind[i]++;
    }
  }
  PASSERT(ctrl, pexadj[nvtxs] == nsend);

  SHIFTCSR(i, nvtxs, pexadj);

  WCOREPOP;

  /* Create the inverse map from ladjncy to adjncy */
  imap = graph->imap = imalloc(nvtxs+nrecv, "CommSetup: imap");
  for (i=0; i<nvtxs; i++)
    imap[i] = firstvtx+i;
  for (i=0; i<nrecv; i++)
    imap[nvtxs+i] = recvind[i];

  STOPTIMER(ctrl, ctrl->SetupTmr);

#ifdef DEBUG_SETUPINFO
  rprintf(ctrl, "[%5"PRIDX" %5"PRIDX"] \tl:[%5"PRIDX" %5"PRIDX"] \ts:[%5"PRIDX", %5"PRIDX"] \tr:[%5"PRIDX", %5"PRIDX"]\n", 
            GlobalSEMin(ctrl, nvtxs), GlobalSEMax(ctrl, nvtxs),
            GlobalSEMin(ctrl, nlocal), GlobalSEMax(ctrl, nlocal),
            GlobalSEMin(ctrl, nsend), GlobalSEMax(ctrl, nsend),
            GlobalSEMin(ctrl, nrecv), GlobalSEMax(ctrl, nrecv));

  PrintSetUpInfo(ctrl, graph);
#endif

}


/*************************************************************************/
/*! This function updates the sreq/rreq/statuses arrays in ctrl based on
    the new number of neighbors.
*/
/*************************************************************************/
void CommUpdateNnbrs(ctrl_t *ctrl, idx_t nnbrs)
{
  if (ctrl->ncommpes >= nnbrs)
    return;

  ctrl->ncommpes = nnbrs;
  ctrl->sreq = (MPI_Request *)gk_realloc(ctrl->sreq, sizeof(MPI_Request)*nnbrs, "sreq");
  ctrl->rreq = (MPI_Request *)gk_realloc(ctrl->rreq, sizeof(MPI_Request)*nnbrs, "rreq");
  ctrl->statuses = (MPI_Status *)gk_realloc(ctrl->statuses, sizeof(MPI_Status)*nnbrs, "statuses");

}


/*************************************************************************/
/*! This function performs the gather/scatter for the boundary vertices 
*/
/*************************************************************************/
void CommInterfaceData(ctrl_t *ctrl, graph_t *graph, idx_t *data, 
         idx_t *recvvector)
{
  idx_t i, k, nnbrs, firstvtx;
  idx_t *peind, *sendptr, *sendind, *sendvector, *recvptr, *recvind;

  WCOREPUSH;

  firstvtx = graph->vtxdist[ctrl->mype];
  nnbrs    = graph->nnbrs;
  peind    = graph->peind;
  sendptr  = graph->sendptr;
  sendind  = graph->sendind;
  recvptr  = graph->recvptr;
  recvind  = graph->recvind;

  /* Issue the receives first */
  for (i=0; i<nnbrs; i++) {
    gkMPI_Irecv((void *)(recvvector+recvptr[i]), recvptr[i+1]-recvptr[i], IDX_T, 
              peind[i], 1, ctrl->comm, ctrl->rreq+i);
  }

  /* Issue the sends next */
  k = sendptr[nnbrs];
  sendvector = iwspacemalloc(ctrl, k);
  for (i=0; i<k; i++) 
    sendvector[i] = data[sendind[i]-firstvtx];

  for (i=0; i<nnbrs; i++) {
    gkMPI_Isend((void *)(sendvector+sendptr[i]), sendptr[i+1]-sendptr[i], IDX_T, 
              peind[i], 1, ctrl->comm, ctrl->sreq+i); 
  }

  /* OK, now get into the loop waiting for the operations to finish */
  gkMPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses); 
  gkMPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses); 

  WCOREPOP;
}



/*************************************************************************
* This function performs the gather/scatter for the boundary vertices
**************************************************************************/
void CommChangedInterfaceData(ctrl_t *ctrl, graph_t *graph, idx_t nchanged, 
         idx_t *changed, idx_t *data, ikv_t *sendpairs, ikv_t *recvpairs)
{
  idx_t i, j, k, nnbrs, firstvtx, nrecv, penum, nreceived;
  idx_t *peind, *sendptr, *recvptr, *recvind, *pexadj, *peadjncy, 
        *peadjloc, *psendptr;
  ikv_t *pairs;


  firstvtx = graph->vtxdist[ctrl->mype];
  nnbrs    = graph->nnbrs;
  nrecv    = graph->nrecv;
  peind    = graph->peind;
  sendptr  = graph->sendptr;
  recvptr  = graph->recvptr;
  recvind  = graph->recvind;
  pexadj   = graph->pexadj;
  peadjncy = graph->peadjncy;
  peadjloc = graph->peadjloc;

  /* Issue the receives first */
  for (i=0; i<nnbrs; i++) {
    gkMPI_Irecv((void *)(recvpairs+recvptr[i]), 2*(recvptr[i+1]-recvptr[i]), IDX_T, 
              peind[i], 1, ctrl->comm, ctrl->rreq+i);
  }

  if (nchanged != 0) {
    WCOREPUSH;

    psendptr = icopy(nnbrs, sendptr, iwspacemalloc(ctrl, nnbrs));

    /* Copy the changed values into the sendvector */
    for (i=0; i<nchanged; i++) {
      j = changed[i];
      for (k=pexadj[j]; k<pexadj[j+1]; k++) {
        penum = peadjncy[k];
        sendpairs[psendptr[penum]].key = peadjloc[k];
        sendpairs[psendptr[penum]].val = data[j];
        psendptr[penum]++;
      }
    }

    for (i=0; i<nnbrs; i++) {
      gkMPI_Isend((void *)(sendpairs+sendptr[i]), 2*(psendptr[i]-sendptr[i]), IDX_T, 
                peind[i], 1, ctrl->comm, ctrl->sreq+i);
    }

    WCOREPOP;
  }
  else {
    for (i=0; i<nnbrs; i++) 
      gkMPI_Isend((void *)(sendpairs), 0, IDX_T, peind[i], 1, ctrl->comm, ctrl->sreq+i);
  }


  /* OK, now get into the loop waiting for the operations to finish */
  for (i=0; i<nnbrs; i++) {
    gkMPI_Wait(ctrl->rreq+i, &(ctrl->status));
    gkMPI_Get_count(&ctrl->status, IDX_T, &nreceived);
    if (nreceived != 0) {
      nreceived = nreceived/2;
      pairs     = recvpairs+graph->recvptr[i];
      for (k=0; k<nreceived; k++) 
        data[pairs[k].key] = pairs[k].val;
    }
  }

  gkMPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);
}



/*************************************************************************
* This function computes the max of a single element
**************************************************************************/
idx_t GlobalSEMax(ctrl_t *ctrl, idx_t value)
{
  idx_t max;

  gkMPI_Allreduce((void *)&value, (void *)&max, 1, IDX_T, MPI_MAX, ctrl->comm);

  return max;
}

/*************************************************************************
* This function computes the max of a single element
**************************************************************************/
idx_t GlobalSEMaxComm(MPI_Comm comm, idx_t value)
{
  idx_t max;

  gkMPI_Allreduce((void *)&value, (void *)&max, 1, IDX_T, MPI_MAX, comm);

  return max;
}

/*************************************************************************
* This function computes the max of a single element
**************************************************************************/
idx_t GlobalSEMin(ctrl_t *ctrl, idx_t value)
{
  idx_t min;

  gkMPI_Allreduce((void *)&value, (void *)&min, 1, IDX_T, MPI_MIN, ctrl->comm);

  return min;
}

/*************************************************************************
* This function computes the max of a single element
**************************************************************************/
idx_t GlobalSEMinComm(MPI_Comm comm, idx_t value)
{
  idx_t min;

  gkMPI_Allreduce((void *)&value, (void *)&min, 1, IDX_T, MPI_MIN, comm);

  return min;
}


/*************************************************************************
* This function computes the max of a single element
**************************************************************************/
idx_t GlobalSESum(ctrl_t *ctrl, idx_t value)
{
  idx_t sum;

  gkMPI_Allreduce((void *)&value, (void *)&sum, 1, IDX_T, MPI_SUM, ctrl->comm);

  return sum;
}

/*************************************************************************
* This function computes the max of a single element
**************************************************************************/
idx_t GlobalSESumComm(MPI_Comm comm, idx_t value)
{
  idx_t min;

  gkMPI_Allreduce((void *)&value, (void *)&min, 1, IDX_T, MPI_SUM, comm);

  return min;
}



/*************************************************************************
* This function computes the max of a single element
**************************************************************************/
real_t GlobalSEMaxFloat(ctrl_t *ctrl, real_t value)
{
  real_t max;

  gkMPI_Allreduce((void *)&value, (void *)&max, 1, REAL_T, MPI_MAX, ctrl->comm);

  return max;
}



/*************************************************************************
* This function computes the max of a single element
**************************************************************************/
real_t GlobalSEMinFloat(ctrl_t *ctrl, real_t value)
{
  real_t min;

  gkMPI_Allreduce((void *)&value, (void *)&min, 1, REAL_T, MPI_MIN, ctrl->comm);

  return min;
}

/*************************************************************************
* This function computes the max of a single element
**************************************************************************/
real_t GlobalSESumFloat(ctrl_t *ctrl, real_t value)
{
  real_t sum;

  gkMPI_Allreduce((void *)&value, (void *)&sum, 1, REAL_T, MPI_SUM, ctrl->comm);

  return sum;
}

