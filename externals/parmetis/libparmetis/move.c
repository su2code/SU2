/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mmove.c
 *
 * This file contains functions that move the graph given a partition
 *
 * Started 11/22/96
 * George
 *
 * $Id: move.c 10657 2011-08-03 14:34:35Z karypis $
 *
 */

#include <parmetislib.h>

/*************************************************************************/
/*! This function moves the graph, and returns a new graph.
    This routine can be called with or without performing refinement.
    In the latter case it allocates and computes lpwgts itself.
*/
/*************************************************************************/
graph_t *MoveGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t h, i, ii, j, jj, nvtxs, ncon, npes, nsnbrs, nrnbrs;
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *mvtxdist;
  idx_t *where, *newlabel, *lpwgts, *gpwgts;
  idx_t *sgraph, *rgraph;
  ikv_t *sinfo, *rinfo;
  graph_t *mgraph;

  WCOREPUSH;

  /* this routine only works when nparts <= npes */
  PASSERT(ctrl, ctrl->nparts <= ctrl->npes);

  npes = ctrl->npes;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  mvtxdist = imalloc(npes+1, "MoveGraph: mvtxdist");

  /* Let's do a prefix scan to determine the labeling of the nodes given */
  lpwgts = iwspacemalloc(ctrl, npes+1);
  gpwgts = iwspacemalloc(ctrl, npes+1);
  sinfo  = ikvwspacemalloc(ctrl, npes);
  rinfo  = ikvwspacemalloc(ctrl, npes);

  for (i=0; i<npes; i++)
    sinfo[i].key = sinfo[i].val = 0;

  for (i=0; i<nvtxs; i++) {
    sinfo[where[i]].key++;
    sinfo[where[i]].val += xadj[i+1]-xadj[i];
  }
  for (i=0; i<npes; i++)
    lpwgts[i] = sinfo[i].key;

  gkMPI_Scan((void *)lpwgts, (void *)gpwgts, npes, IDX_T, MPI_SUM, ctrl->comm);
  gkMPI_Allreduce((void *)lpwgts, (void *)mvtxdist, npes, IDX_T, MPI_SUM, ctrl->comm);
  MAKECSR(i, npes, mvtxdist);

  /* gpwgts[i] will store the label of the first vertex for each domain 
     in each processor */
  for (i=0; i<npes; i++) 
    /* We were interested in an exclusive scan */
    gpwgts[i] = mvtxdist[i] + gpwgts[i] - lpwgts[i];

  newlabel = iwspacemalloc(ctrl, nvtxs+graph->nrecv);
  for (i=0; i<nvtxs; i++) 
    newlabel[i] = gpwgts[where[i]]++;

  /* OK, now send the newlabel info to processors storing adjacent interface nodes */
  CommInterfaceData(ctrl, graph, newlabel, newlabel+nvtxs);

  /* Now lets tell everybody what and from where he will get it. */
  gkMPI_Alltoall((void *)sinfo, 2, IDX_T, (void *)rinfo, 2, IDX_T, ctrl->comm);

  /* Use lpwgts and gpwgts as pointers to where data will be received and send */
  lpwgts[0] = 0;  /* Send part */
  gpwgts[0] = 0;  /* Received part */
  for (nsnbrs=nrnbrs=0, i=0; i<npes; i++) {
    lpwgts[i+1] = lpwgts[i] + (1+ncon)*sinfo[i].key + 2*sinfo[i].val;
    gpwgts[i+1] = gpwgts[i] + (1+ncon)*rinfo[i].key + 2*rinfo[i].val;
    if (rinfo[i].key > 0)
      nrnbrs++;
    if (sinfo[i].key > 0)
      nsnbrs++;
  }

  /* Update the max # of sreq/rreq/statuses */
  CommUpdateNnbrs(ctrl, gk_max(nsnbrs, nrnbrs));

  rgraph = iwspacemalloc(ctrl, gpwgts[npes]);
  WCOREPUSH;  /* for freeing the send part early */
  sgraph = iwspacemalloc(ctrl, lpwgts[npes]);

  /* Issue the receives first */
  for (j=0, i=0; i<npes; i++) {
    if (rinfo[i].key > 0) 
      gkMPI_Irecv((void *)(rgraph+gpwgts[i]), gpwgts[i+1]-gpwgts[i], IDX_T, 
          i, 1, ctrl->comm, ctrl->rreq+j++);
    else 
      PASSERT(ctrl, gpwgts[i+1]-gpwgts[i] == 0);
  }

  /* Assemble the graph to be sent and send it */
  for (i=0; i<nvtxs; i++) {
    PASSERT(ctrl, where[i] >= 0 && where[i] < npes);
    ii = lpwgts[where[i]];
    sgraph[ii++] = xadj[i+1]-xadj[i];
    for (h=0; h<ncon; h++)
      sgraph[ii++] = vwgt[i*ncon+h];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      sgraph[ii++] = newlabel[adjncy[j]];
      sgraph[ii++] = adjwgt[j];
    }
    lpwgts[where[i]] = ii;
  }

  SHIFTCSR(i, npes, lpwgts);

  for (j=0, i=0; i<npes; i++) {
    if (sinfo[i].key > 0)
      gkMPI_Isend((void *)(sgraph+lpwgts[i]), lpwgts[i+1]-lpwgts[i], IDX_T, 
          i, 1, ctrl->comm, ctrl->sreq+j++);
    else 
      PASSERT(ctrl, lpwgts[i+1]-lpwgts[i] == 0);
  }

  /* Wait for the send/recv to finish */
  gkMPI_Waitall(nrnbrs, ctrl->rreq, ctrl->statuses);
  gkMPI_Waitall(nsnbrs, ctrl->sreq, ctrl->statuses);

  WCOREPOP;  /* frees sgraph */

  /* OK, now go and put the graph into graph_t Format */
  mgraph = CreateGraph();
  
  mgraph->vtxdist = mvtxdist;
  mgraph->gnvtxs  = graph->gnvtxs;
  mgraph->ncon    = ncon;
  mgraph->level   = 0;
  mgraph->nvtxs   = mgraph->nedges = 0;
  for (i=0; i<npes; i++) {
    mgraph->nvtxs  += rinfo[i].key;
    mgraph->nedges += rinfo[i].val;
  }
  nvtxs  = mgraph->nvtxs;
  xadj   = mgraph->xadj   = imalloc(nvtxs+1, "MMG: mgraph->xadj");
  vwgt   = mgraph->vwgt   = imalloc(nvtxs*ncon, "MMG: mgraph->vwgt");
  adjncy = mgraph->adjncy = imalloc(mgraph->nedges, "MMG: mgraph->adjncy");
  adjwgt = mgraph->adjwgt = imalloc(mgraph->nedges, "MMG: mgraph->adjwgt");

  for (jj=ii=i=0; i<nvtxs; i++) {
    xadj[i] = rgraph[ii++];
    for (h=0; h<ncon; h++)
      vwgt[i*ncon+h] = rgraph[ii++]; 
    for (j=0; j<xadj[i]; j++, jj++) {
      adjncy[jj] = rgraph[ii++];
      adjwgt[jj] = rgraph[ii++];
    }
  }
  MAKECSR(i, nvtxs, xadj);

  PASSERT(ctrl, jj == mgraph->nedges);
  PASSERT(ctrl, ii == gpwgts[npes]);
  PASSERTP(ctrl, jj == mgraph->nedges, (ctrl, "%"PRIDX" %"PRIDX"\n", jj, mgraph->nedges));
  PASSERTP(ctrl, ii == gpwgts[npes], (ctrl, "%"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"\n", 
      ii, gpwgts[npes], jj, mgraph->nedges, nvtxs));

#ifdef DEBUG
  IFSET(ctrl->dbglvl, DBG_INFO, rprintf(ctrl, "Checking moved graph...\n"));
  CheckMGraph(ctrl, mgraph);
  IFSET(ctrl->dbglvl, DBG_INFO, rprintf(ctrl, "Moved graph is consistent.\n"));
#endif

  WCOREPOP;

  return mgraph;
}



/*************************************************************************
* This function is used to transfer information from the moved graph
* back to the original graph. The information is transfered from array
* minfo to array info. The routine assumes that graph->where is left intact
* and it is used to get the inverse mapping information.
* The routine assumes that graph->where corresponds to a npes-way partition.
**************************************************************************/
void ProjectInfoBack(ctrl_t *ctrl, graph_t *graph, idx_t *info, idx_t *minfo)
{
  idx_t i, nvtxs, nparts, nrecvs, nsends;
  idx_t *where, *auxinfo, *sinfo, *rinfo;

  WCOREPUSH;

  nparts = ctrl->npes;

  nvtxs = graph->nvtxs;
  where = graph->where;

  sinfo = iwspacemalloc(ctrl, nparts+1);
  rinfo = iwspacemalloc(ctrl, nparts+1);

  /* Find out in rinfo how many entries are received per partition */
  iset(nparts, 0, rinfo);
  for (i=0; i<nvtxs; i++)
    rinfo[where[i]]++;

  /* The rinfo are transposed and become the sinfo for the back-projection */
  gkMPI_Alltoall((void *)rinfo, 1, IDX_T, (void *)sinfo, 1, IDX_T, ctrl->comm);

  MAKECSR(i, nparts, sinfo);
  MAKECSR(i, nparts, rinfo);

  /* allocate memory for auxinfo */
  auxinfo = iwspacemalloc(ctrl, rinfo[nparts]);

  /*-----------------------------------------------------------------
   * Now, go and send back the minfo
   -----------------------------------------------------------------*/
  for (nrecvs=0, i=0; i<nparts; i++) {
    if (rinfo[i+1]-rinfo[i] > 0)
      gkMPI_Irecv((void *)(auxinfo+rinfo[i]), rinfo[i+1]-rinfo[i], IDX_T, 
          i, 1, ctrl->comm, ctrl->rreq+nrecvs++);
  }

  for (nsends=0, i=0; i<nparts; i++) {
    if (sinfo[i+1]-sinfo[i] > 0) 
      gkMPI_Isend((void *)(minfo+sinfo[i]), sinfo[i+1]-sinfo[i], IDX_T, 
          i, 1, ctrl->comm, ctrl->sreq+nsends++);
  }
  PASSERT(ctrl, nrecvs <= ctrl->ncommpes);
  PASSERT(ctrl, nsends <= ctrl->ncommpes);

  /* Wait for the send/recv to finish */
  gkMPI_Waitall(nrecvs, ctrl->rreq, ctrl->statuses);
  gkMPI_Waitall(nsends, ctrl->sreq, ctrl->statuses);

  /* Scatter the info received in auxinfo back to info. */
  for (i=0; i<nvtxs; i++)
    info[i] = auxinfo[rinfo[where[i]]++];

  WCOREPOP;
}



/*************************************************************************
* This function is used to convert a partition vector to a permutation
* vector.
**************************************************************************/
void FindVtxPerm(ctrl_t *ctrl, graph_t *graph, idx_t *perm)
{
  idx_t i, nvtxs, nparts;
  idx_t *xadj, *adjncy, *adjwgt, *mvtxdist;
  idx_t *where, *lpwgts, *gpwgts;

  WCOREPUSH;

  nparts = ctrl->nparts;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  mvtxdist = iwspacemalloc(ctrl, nparts+1);
  lpwgts   = iwspacemalloc(ctrl, nparts+1);
  gpwgts   = iwspacemalloc(ctrl, nparts+1);


  /* Here we care about the count and not total weight (diff since graph may 
     be weighted */
  iset(nparts, 0, lpwgts);
  for (i=0; i<nvtxs; i++)
    lpwgts[where[i]]++;

  /* Let's do a prefix scan to determine the labeling of the nodes given */
  gkMPI_Scan((void *)lpwgts, (void *)gpwgts, nparts, IDX_T, MPI_SUM, ctrl->comm);
  gkMPI_Allreduce((void *)lpwgts, (void *)mvtxdist, nparts, IDX_T, MPI_SUM, ctrl->comm);

  MAKECSR(i, nparts, mvtxdist);

  for (i=0; i<nparts; i++)
    gpwgts[i] = mvtxdist[i] + gpwgts[i] - lpwgts[i];  /* We were interested in an exclusive Scan */

  for (i=0; i<nvtxs; i++)
    perm[i] = gpwgts[where[i]]++;

  WCOREPOP;
}


/*************************************************************************
* This function quickly performs a check on the consistency of moved graph.
**************************************************************************/
void CheckMGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, jj, k, nvtxs, firstvtx, lastvtx;
  idx_t *xadj, *adjncy, *vtxdist;

  nvtxs   = graph->nvtxs;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  vtxdist = graph->vtxdist;

  firstvtx = vtxdist[ctrl->mype];
  lastvtx  = vtxdist[ctrl->mype+1];

  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (firstvtx+i == adjncy[j])
        myprintf(ctrl, "(%"PRIDX" %"PRIDX") diagonal entry\n", i, i);

      if (adjncy[j] >= firstvtx && adjncy[j] < lastvtx) {
        k = adjncy[j]-firstvtx;
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          if (adjncy[jj] == firstvtx+i)
            break;
        }
        if (jj == xadj[k+1])
          myprintf(ctrl, "(%"PRIDX" %"PRIDX") but not (%"PRIDX" %"PRIDX") [%"PRIDX" %"PRIDX"] [%"PRIDX" %"PRIDX"]\n", 
              i, k, k, i, firstvtx+i, firstvtx+k, 
              xadj[i+1]-xadj[i], xadj[k+1]-xadj[k]);
      }
    }
  }
}



