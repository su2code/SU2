/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * node_refine.c
 *
 * This file contains code that performs the k-way refinement
 *
 * Started 3/1/96
 * George
 *
 * $Id: node_refine.c 10391 2011-06-23 19:00:08Z karypis $
 */

#include <parmetislib.h>


/************************************************************************************/
/*! 
  This function allocates the memory required for the nodeND refinement code.
  The only refinement-related information that is has is the \c graph->where vector
  and allocates the memory for the remaining of the refinement-related data-structures.

*/
/************************************************************************************/
void AllocateNodePartitionParams(ctrl_t *ctrl, graph_t *graph)
{
  idx_t nparts, nvtxs;
  idx_t *vwgt;
  NRInfoType *rinfo, *myrinfo;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayInitTmr));

  nvtxs  = graph->nvtxs;
  nparts = ctrl->nparts;

  graph->nrinfo  = (NRInfoType *)gk_malloc(sizeof(NRInfoType)*nvtxs, 
                                     "AllocateNodePartitionParams: rinfo");
  graph->lpwgts  = imalloc(2*nparts, "AllocateNodePartitionParams: lpwgts");
  graph->gpwgts  = imalloc(2*nparts, "AllocateNodePartitionParams: gpwgts");
  graph->sepind  = imalloc(nvtxs, "AllocateNodePartitionParams: sepind");

  /* Allocate additional memory for graph->vwgt in order to store the weights
     of the remote vertices */
  vwgt        = graph->vwgt;
  graph->vwgt = imalloc(nvtxs+graph->nrecv, "AllocateNodePartitionParams: graph->vwgt");
  icopy(nvtxs, vwgt, graph->vwgt);
  gk_free((void **)&vwgt, LTERM);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayInitTmr));
}


/************************************************************************************/
/*! 
  This function computes the initial node refinment information for the parallel
  nodeND code. It requires that the required data-structures have already been
  allocated via a call to AllocateNodePartitionParams.

*/
/************************************************************************************/
void ComputeNodePartitionParams(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, nparts, nvtxs, nsep;
  idx_t *xadj, *adjncy, *adjwgt, *vtxdist, *vwgt, *lpwgts, *gpwgts, *sepind;
  idx_t *where;
  NRInfoType *rinfo, *myrinfo;
  idx_t me, other, otherwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayInitTmr));

  nvtxs  = graph->nvtxs;
  nparts = ctrl->nparts;

  vtxdist = graph->vtxdist;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  vwgt    = graph->vwgt;

  where   = graph->where;
  rinfo   = graph->nrinfo;
  lpwgts  = graph->lpwgts;
  gpwgts  = graph->gpwgts;
  sepind  = graph->sepind;

  /* Reset refinement data structures */
  iset(2*nparts, 0, lpwgts);

  /* Send/Receive the where and vwgt information of interface vertices. */
  CommInterfaceData(ctrl, graph, where, where+nvtxs); 
  CommInterfaceData(ctrl, graph, vwgt, vwgt+nvtxs); 


  /*------------------------------------------------------------
  / Compute now the degrees
  /------------------------------------------------------------*/
  for (nsep=i=0; i<nvtxs; i++) {
    me = where[i];
    PASSERT(ctrl, me >= 0 && me < 2*nparts);
    lpwgts[me] += vwgt[i];

    if (me >= nparts) {  /* If it is a separator vertex */
      sepind[nsep++] = i;
      lpwgts[2*nparts-1] += vwgt[i];  /* Keep track of total separator weight */

      myrinfo = rinfo+i;
      myrinfo->edegrees[0] = myrinfo->edegrees[1] = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (me != other)
          myrinfo->edegrees[other%2] += vwgt[adjncy[j]];
      }
    }
  }
  graph->nsep = nsep;

  /* Finally, sum-up the partition weights */
  gkMPI_Allreduce((void *)lpwgts, (void *)gpwgts, 2*nparts, IDX_T, MPI_SUM, ctrl->comm);
  graph->mincut = gpwgts[2*nparts-1];

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayInitTmr));
}


/************************************************************************************/
/*! 
  This function updates the node refinment information after a where[] change 

*/
/************************************************************************************/
void UpdateNodePartitionParams(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, nparts, nvtxs, nsep;
  idx_t *xadj, *adjncy, *adjwgt, *vtxdist, *vwgt, *lpwgts, *gpwgts, *sepind;
  idx_t *where;
  NRInfoType *rinfo, *myrinfo;
  idx_t me, other, otherwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayInitTmr));

  nvtxs  = graph->nvtxs;
  nparts = ctrl->nparts;

  vtxdist = graph->vtxdist;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  vwgt    = graph->vwgt;

  where   = graph->where;
  rinfo   = graph->nrinfo;
  lpwgts  = graph->lpwgts;
  gpwgts  = graph->gpwgts;
  sepind  = graph->sepind;

  /* Reset refinement data structures */
  iset(2*nparts, 0, lpwgts);

  /* Send/Receive the where and vwgt information of interface vertices. */
  CommInterfaceData(ctrl, graph, where, where+nvtxs); 


  /*------------------------------------------------------------
  / Compute now the degrees
  /------------------------------------------------------------*/
  for (nsep=i=0; i<nvtxs; i++) {
    me = where[i];
    PASSERT(ctrl, me >= 0 && me < 2*nparts);
    lpwgts[me] += vwgt[i];

    if (me >= nparts) {  /* If it is a separator vertex */
      sepind[nsep++] = i;
      lpwgts[2*nparts-1] += vwgt[i];  /* Keep track of total separator weight */

      myrinfo = rinfo+i;
      myrinfo->edegrees[0] = myrinfo->edegrees[1] = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (me != other)
          myrinfo->edegrees[other%2] += vwgt[adjncy[j]];
      }
    }
  }
  graph->nsep = nsep;

  /* Finally, sum-up the partition weights */
  gkMPI_Allreduce((void *)lpwgts, (void *)gpwgts, 2*nparts, IDX_T, MPI_SUM, ctrl->comm);
  graph->mincut = gpwgts[2*nparts-1];

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayInitTmr));
}



/************************************************************************************/
/*! 
  This function performs k-way node-based refinement. The refinement is done
  concurrently for all the different partitions. It works because each of the
  partitions is disconnected from each other due to the removal of the previous level
  separators.

  This version uses a priority queue to order the nodes and incorporates gain updates 
  from the local information within each inner iteration.

  The refinement exits when there is no improvement in two successive itereations
  in order to account for the fact that a 0 => 1 iteration may have no gain but a
  1 => 0 iteration may have a gain.
*/
/************************************************************************************/
void KWayNodeRefine_Greedy(ctrl_t *ctrl, graph_t *graph, idx_t npasses, real_t ubfrac)
{
  idx_t i, ii, iii, j, jj, k, pass, nvtxs, nrecv, firstvtx, lastvtx, otherlastvtx, 
      side, c, cc, nmoves, nlupd, nsupd, nnbrs, nchanged, nsep, nzerogainiterations;
  idx_t npes = ctrl->npes, mype = ctrl->mype, nparts = ctrl->nparts;
  idx_t *xadj, *adjncy, *adjwgt, *vtxdist, *vwgt;
  idx_t *where, *lpwgts, *gpwgts, *sepind;
  idx_t *peind, *recvptr, *sendptr;
  idx_t *update, *supdate, *rupdate, *pe_updates, *marker, *changed;
  idx_t *badmaxpwgt;
  ikv_t *swchanges, *rwchanges;
  idx_t *nupds_pe;
  NRInfoType *rinfo, *myrinfo;
  idx_t from, to, me, other, oldcut;
  rpq_t *queue;
  idx_t *inqueue;
  idx_t *rxadj, *radjncy;
  char title[1024];

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayTmr));

  WCOREPUSH;

  nvtxs = graph->nvtxs;
  nrecv = graph->nrecv;

  vtxdist = graph->vtxdist;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  vwgt    = graph->vwgt;

  firstvtx = vtxdist[mype];
  lastvtx  = vtxdist[mype+1];

  where  = graph->where;
  rinfo  = graph->nrinfo;
  lpwgts = graph->lpwgts;
  gpwgts = graph->gpwgts;

  nsep   = graph->nsep;
  sepind = graph->sepind;

  nnbrs   = graph->nnbrs;
  peind   = graph->peind;
  recvptr = graph->recvptr;
  sendptr = graph->sendptr;

  badmaxpwgt = iwspacemalloc(ctrl, nparts);
  nupds_pe   = iwspacemalloc(ctrl, npes);
  changed    = iwspacemalloc(ctrl, nvtxs);
  update     = iwspacemalloc(ctrl, nvtxs);
  marker     = iset(nvtxs+nrecv, 0, iwspacemalloc(ctrl, nvtxs+nrecv));
  inqueue    = iwspacemalloc(ctrl, nvtxs+nrecv);
  rwchanges  = ikvwspacemalloc(ctrl, graph->nrecv);
  swchanges  = ikvwspacemalloc(ctrl, graph->nsend);
  supdate    = iwspacemalloc(ctrl, graph->nrecv);
  rupdate    = iwspacemalloc(ctrl, graph->nsend);

  queue = rpqCreate(nvtxs);

  for (i=0; i<nparts; i+=2) {
    //badmaxpwgt[i] = badmaxpwgt[i+1] = ubfrac*(gpwgts[i]+gpwgts[i+1]+gpwgts[nparts+i])/2;
    badmaxpwgt[i] = badmaxpwgt[i+1] = ubfrac*(gpwgts[i]+gpwgts[i+1])/2;
  }

  /* construct the local adjacency list of the interface nodes */
  rxadj = iset(nrecv+1, 0, iwspacemalloc(ctrl, nrecv+1));
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if ((k = adjncy[j]-nvtxs) >= 0)
        rxadj[k]++;
    }
  }
  MAKECSR(i, nrecv, rxadj);

  radjncy = iwspacemalloc(ctrl, rxadj[nrecv]);
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if ((k = adjncy[j]-nvtxs) >= 0)
        radjncy[rxadj[k]++] = i;
    }
  }
  SHIFTCSR(i, nrecv, rxadj);


  IFSET(ctrl->dbglvl, DBG_REFINEINFO, 
      PrintNodeBalanceInfo(ctrl, nparts, gpwgts, badmaxpwgt, "K-way sep-refinement"));

  c = GlobalSESum(ctrl, RandomInRange(2))%2;
  nzerogainiterations = 0;
  for (pass=0; pass<npasses; pass++) {
    oldcut = graph->mincut;

    for (side=0; side<2; side++) {
      cc = (c+1)%2;

#ifndef NDEBUG
      /* This is for debugging purposes */
      for (ii=0, i=0; i<nvtxs; i++) {
        if (where[i] >= nparts)
          ii++;
      }
      PASSERT(ctrl, ii == nsep);
#endif

      /* Put the separator nodes in queue */
      rpqReset(queue);
      iset(nvtxs+nrecv, 0, inqueue);
      for (ii=0; ii<nsep; ii++) {
        i = sepind[ii];
        PASSERT(ctrl, inqueue[i] == 0);
        PASSERT(ctrl, where[i] >= nparts);
        rpqInsert(queue, i, vwgt[i] - rinfo[i].edegrees[cc]);
        inqueue[i] = 1;
      }

      nlupd = nsupd = nmoves = nchanged = nsep = 0;
      while ((i = rpqGetTop(queue)) != -1) {
        inqueue[i] = 0; 

        from = where[i];
        PASSERT(ctrl, from >= nparts);

        /* It is a one-sided move so it will go to the other partition. 
           Look at the comments in InitMultisection to understand the meaning 
           of from%nparts */
        to    = from%nparts+c;   /* where to move the separator node */
        other = from%nparts+cc;  /* the other partition involved in the 3-way view */

        /* Go through the loop to see if gain is possible for the separator vertex */
        if (gpwgts[to]+vwgt[i] <= badmaxpwgt[to] && vwgt[i] - rinfo[i].edegrees[cc] >= 0) {
          /* Update the where information of the vertex you moved */
          where[i] = to;

          lpwgts[from]       -= vwgt[i];
          lpwgts[2*nparts-1] -= vwgt[i];
          lpwgts[to]         += vwgt[i];
          gpwgts[to]         += vwgt[i];

          /* Update and record future updates */
          for (j=xadj[i]; j<xadj[i+1]; j++) {
            ii = adjncy[j];

            /* If vertex ii is being pulled into the separator for the first time, 
               then update the edegrees[] of the nodes currently in the queue that
               vertex ii is connected to */
            if (marker[ii] == 0 && where[ii] == other) {
              if (ii < nvtxs) { /* local vertices */
                for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
                  iii = adjncy[jj];
                  if (inqueue[iii] == 1) {
                    PASSERT(ctrl, rinfo[iii].edegrees[cc] >= vwgt[ii]);
                    rinfo[iii].edegrees[cc] -= vwgt[ii];
                    rpqUpdate(queue, iii, vwgt[iii]-rinfo[iii].edegrees[cc]);
                  }
                }
              }
              else { /* remote vertices */
                for (jj=rxadj[ii-nvtxs]; jj<rxadj[ii-nvtxs+1]; jj++) {
                  iii = radjncy[jj];
                  if (inqueue[iii] == 1) {
                    PASSERT(ctrl, rinfo[iii].edegrees[cc] >= vwgt[ii]);
                    rinfo[iii].edegrees[cc] -= vwgt[ii];
                    rpqUpdate(queue, iii, vwgt[iii]-rinfo[iii].edegrees[cc]);
                  }
                }
              }
            }


            /* Put the vertices adjacent to i that belong to either the separator or
               the cc partition into the update array */
            if (marker[ii] == 0 && where[ii] != to) {
              marker[ii] = 1;
              if (ii<nvtxs)
                update[nlupd++] = ii;
              else
                supdate[nsupd++] = ii;
            }
          }
          nmoves++;
          if (graph->pexadj[i+1]-graph->pexadj[i] > 0)
            changed[nchanged++] = i;
        }
        else {
          /* This node will remain in the separator for the next iteration */
          sepind[nsep++] = i;
        }
      }

      /* myprintf(ctrl, "nmoves: %"PRIDX", nlupd: %"PRIDX", nsupd: %"PRIDX"\n", nmoves, nlupd, nsupd); */

      IFSET(ctrl->dbglvl, DBG_RMOVEINFO, rprintf(ctrl, 
            "\t[%"PRIDX" %"PRIDX"], [%"PRIDX" %"PRIDX" %"PRIDX"]\n", 
            pass, c, GlobalSESum(ctrl, nmoves), GlobalSESum(ctrl, nsupd), 
            GlobalSESum(ctrl, nlupd)));


      /*-----------------------------------------------------------------------
      / Time to communicate with processors to send the vertices whose degrees 
      / need to be updated.
      /-----------------------------------------------------------------------*/
      /* Issue the receives first */
      for (i=0; i<nnbrs; i++) {
        gkMPI_Irecv((void *)(rupdate+sendptr[i]), sendptr[i+1]-sendptr[i], IDX_T,
            peind[i], 1, ctrl->comm, ctrl->rreq+i);
      }

      /* Issue the sends next. This needs some preporcessing */
      for (i=0; i<nsupd; i++) {
        marker[supdate[i]] = 0;
        supdate[i] = graph->imap[supdate[i]];
      }
      isorti(nsupd, supdate);

      for (j=i=0; i<nnbrs; i++) {
        otherlastvtx = vtxdist[peind[i]+1];
        for (k=j; k<nsupd && supdate[k] < otherlastvtx; k++); 
        gkMPI_Isend((void *)(supdate+j), k-j, IDX_T, peind[i], 1, ctrl->comm, 
            ctrl->sreq+i);
        j = k;
      }

      /* OK, now get into the loop waiting for the send/recv operations to finish */
      gkMPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
      for (i=0; i<nnbrs; i++) 
        gkMPI_Get_count(ctrl->statuses+i, IDX_T, nupds_pe+i);
      gkMPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


      /*-------------------------------------------------------------
      / Place the received to-be updated vertices into update[] 
      /-------------------------------------------------------------*/
      for (i=0; i<nnbrs; i++) {
        pe_updates = rupdate+sendptr[i];
        for (j=0; j<nupds_pe[i]; j++) {
          k = pe_updates[j];
          if (marker[k-firstvtx] == 0) {
            marker[k-firstvtx] = 1;
            update[nlupd++] = k-firstvtx;
          }
        }
      }


      /*-------------------------------------------------------------
      / Update the where information of the vertices that are pulled
      / into the separator.
      /-------------------------------------------------------------*/
      for (ii=0; ii<nlupd; ii++) {
        i = update[ii];
        me = where[i];
        if (me < nparts && me%2 == cc) { /* This vertex is pulled into the separator */
          lpwgts[me] -= vwgt[i];
          where[i] = nparts+me-(me%2); 
          sepind[nsep++] = i;  /* Put the vertex into the sepind array */
          if (graph->pexadj[i+1]-graph->pexadj[i] > 0)
            changed[nchanged++] = i;

          lpwgts[where[i]]   += vwgt[i];
          lpwgts[2*nparts-1] += vwgt[i];
          /* myprintf(ctrl, "Vertex %"PRIDX" moves into the separator from %"PRIDX" to %"PRIDX"\n", 
                 i+firstvtx, me, where[i]); */
        }
      }

      /* Tell everybody interested what the new where[] info is for the 
         interface vertices */
      CommChangedInterfaceData(ctrl, graph, nchanged, changed, where, 
          swchanges, rwchanges); 


      /*-------------------------------------------------------------
      / Update the rinfo of the vertices in the update[] array
      /-------------------------------------------------------------*/
      for (ii=0; ii<nlupd; ii++) {
        i = update[ii];
        PASSERT(ctrl, marker[i] == 1);

        marker[i] = 0;

        me = where[i];
        if (me >= nparts) {  /* If it is a separator vertex */
          /* myprintf(ctrl, "Updating %"PRIDX" %"PRIDX"\n", i+firstvtx, me); */

          myrinfo = rinfo+i;
          myrinfo->edegrees[0] = myrinfo->edegrees[1] = 0;

          for (j=xadj[i]; j<xadj[i+1]; j++) {
            other = where[adjncy[j]];
            if (me != other)
              myrinfo->edegrees[other%2] += vwgt[adjncy[j]];
          }
        }
      }

      /* Finally, sum-up the partition weights */
      gkMPI_Allreduce((void *)lpwgts, (void *)gpwgts, 2*nparts, IDX_T, MPI_SUM, 
          ctrl->comm);
      graph->mincut = gpwgts[2*nparts-1];

      sprintf(title, "\tTotalSep [%"PRIDX"]", c);
      IFSET(ctrl->dbglvl, DBG_REFINEINFO, 
          PrintNodeBalanceInfo(ctrl, nparts, gpwgts, badmaxpwgt, title));

      /* break out if there is no improvement in two successive inner iterations that
         can span successive outer iterations */
      if (graph->mincut == oldcut) {
        if (++nzerogainiterations == 2)
          break;
      }
      else {
        nzerogainiterations = 0;
      }

      c = cc;
    }

    /* break out if there is no improvement in two successive inner iterations that
       can span successive outer iterations */
    if (graph->mincut == oldcut && nzerogainiterations == 2) 
      break;
  }

  rpqDestroy(queue);

  WCOREPOP;

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayTmr));
}


/************************************************************************************/
/*! 
  This function performs k-way node-based refinement. The key difference between
  this and the previous routine is that the well-interior nodes are refined using 
  a serial node-based refinement algortihm.
*/
/************************************************************************************/
void KWayNodeRefine2Phase(ctrl_t *ctrl, graph_t *graph, idx_t npasses, real_t ubfrac)
{
  idx_t i, oldcut;

  oldcut = graph->mincut+1;
  for (i=0; i<npasses; i++) {
    KWayNodeRefine_Greedy(ctrl, graph, npasses, ubfrac);
    if (graph->mincut == oldcut)
      break;
    oldcut = graph->mincut;

    KWayNodeRefineInterior(ctrl, graph, 2, ubfrac);
    UpdateNodePartitionParams(ctrl, graph);
    if (graph->mincut == oldcut)
      break;
    oldcut = graph->mincut;
  }
}


/************************************************************************************/
/*! This function performs k-way node-based refinement of the interior nodes of the 
    graph assigned to each processor using a serial node-refinement algorithm. */
/************************************************************************************/
void KWayNodeRefineInterior(ctrl_t *ctrl, graph_t *graph, idx_t npasses, real_t ubfrac)
{
  idx_t i, j, k, ii, gnnz, gid, qsize;
  idx_t npes = ctrl->npes, mype = ctrl->mype, nparts = ctrl->nparts;
  idx_t nvtxs, *xadj, *adjncy, *vwgt, *where, *pexadj;
  idx_t gnvtxs, *gxadj, *gadjncy, *gvwgt, *gwhere, *ghmarker;
  idx_t *gmap, *gimap;
  idx_t *pptr, *pind;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->AuxTmr1));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayTmr));

  WCOREPUSH;

  nvtxs   = graph->nvtxs;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  vwgt    = graph->vwgt;
  where   = graph->where;
  pexadj  = graph->pexadj;

  gxadj    = iwspacemalloc(ctrl, nvtxs+1);
  gvwgt    = iwspacemalloc(ctrl, nvtxs);
  gadjncy  = iwspacemalloc(ctrl, xadj[nvtxs]);
  gwhere   = iwspacemalloc(ctrl, nvtxs);
  ghmarker = iwspacemalloc(ctrl, nvtxs);
  gmap     = iset(nvtxs, -1, iwspacemalloc(ctrl, nvtxs));
  gimap    = iwspacemalloc(ctrl, nvtxs);
  pptr     = iset(2*nparts+1, 0, iwspacemalloc(ctrl, 2*nparts+1));
  pind     = iwspacemalloc(ctrl, nvtxs);


  /* Set pptr/pind to contain the vertices in each one of the partitions */
  for (i=0; i<nvtxs; i++)
    pptr[where[i]]++;
  MAKECSR(i, 2*nparts, pptr);
  for (i=0; i<nvtxs; i++) 
    pind[pptr[where[i]]++] = i;
  SHIFTCSR(i, 2*nparts, pptr);


  /* Extract each individual graph and refine it */
  for (gid=0; gid<nparts; gid+=2) {
    if (graph->lpwgts[nparts+gid] == 0)
      continue;

    /* a quick test to determine if there are sufficient non-interface separator nodes */
    for (qsize=0, ii=pptr[nparts+gid]; ii<pptr[nparts+gid+1]; ii++) {
      if (pexadj[pind[ii]] == pexadj[pind[ii]+1]) {
        if (++qsize >= 2)
          break;
      }
    }
    if (ii == pptr[nparts+gid+1]) 
      break;  /* no need to proceed. not enough non-interface separator nodes */
    
    /* compute the gmap/gimap and other node info for the extracted graph */
    for (gnvtxs=0, ii=pptr[nparts+gid]; ii<pptr[nparts+gid+1]; ii++, gnvtxs++) {
      i                = pind[ii];
      gmap[i]          = gnvtxs;
      gimap[gnvtxs]    = i;
      gvwgt[gnvtxs]    = vwgt[i];
      gwhere[gnvtxs]   = 2;
      ghmarker[gnvtxs] = (pexadj[i+1]-pexadj[i] > 0 ? gwhere[gnvtxs] : -1);
    }

    for (ii=pptr[gid]; ii<pptr[gid+2]; ii++, gnvtxs++) {
      i                = pind[ii];
      gmap[i]          = gnvtxs;
      gimap[gnvtxs]    = i;
      gvwgt[gnvtxs]    = vwgt[i];
      gwhere[gnvtxs]   = where[i] - gid;
      ghmarker[gnvtxs] = (pexadj[i+1]-pexadj[i] > 0 ? gwhere[gnvtxs] : -1);
      PASSERT(ctrl, gwhere[gnvtxs] >= 0 && gwhere[gnvtxs] <= 1);
    }

    gxadj[0]=0; gnvtxs=0; gnnz=0;
    /* go over the separator nodes */
    for (ii=pptr[nparts+gid]; ii<pptr[nparts+gid+1]; ii++) {
      i = pind[ii];
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (adjncy[j] < nvtxs) 
          gadjncy[gnnz++] = gmap[adjncy[j]];
      }
      gxadj[++gnvtxs] = gnnz;
    }

    /* go over the interior nodes */
    for (ii=pptr[gid]; ii<pptr[gid+2]; ii++) {
      i = pind[ii];
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (adjncy[j] < nvtxs) 
          gadjncy[gnnz++] = gmap[adjncy[j]];
      }
      gxadj[++gnvtxs] = gnnz;
    }

    if (gnnz == 0)
      continue;

    /* The 1.03 is here by choice as it is better to refine using 
       tight constraints */
    METIS_NodeRefine(gnvtxs, gxadj, gvwgt, gadjncy, gwhere, ghmarker, 1.03);

    for (i=0; i<gnvtxs; i++) {
      if (gwhere[i] == 2)
        where[gimap[i]] = nparts + gid;
      else
        where[gimap[i]] = gwhere[i] + gid;
    }
  }

  WCOREPOP;

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayTmr));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->AuxTmr1));
}


/************************************************************************************/
/*! This function prints balance information for the parallel k-section refinement 
    algorithm. */
/************************************************************************************/
void PrintNodeBalanceInfo(ctrl_t *ctrl, idx_t nparts, idx_t *gpwgts, idx_t *badmaxpwgt, 
         char *title)
{
  idx_t i;

  if (ctrl->mype == 0) {
    printf("%s: %"PRIDX", ", title, gpwgts[2*nparts-1]);

    for (i=0; i<nparts; i+=2) 
      printf(" [%5"PRIDX" %5"PRIDX" %5"PRIDX" %5"PRIDX"]", gpwgts[i], gpwgts[i+1], gpwgts[nparts+i], badmaxpwgt[i]); 
    printf("\n");
  }
  gkMPI_Barrier(ctrl->comm);
}


