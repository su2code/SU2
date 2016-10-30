/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * initpart.c
 *
 * This file contains code that performs log(p) parallel multilevel
 * recursive bissection
 *
 * Started 3/4/96
 * George
 *
 * $Id: initpart.c 10542 2011-07-11 16:56:22Z karypis $
 */

#include <parmetislib.h>


#define DEBUG_IPART_



/*************************************************************************
* This function is the entry point of the initial partition algorithm
* that does recursive bissection.
* This algorithm assembles the graph to all the processors and preceeds
* by parallelizing the recursive bisection step.
**************************************************************************/
void InitPartition(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, ncon, mype, npes, gnvtxs, ngroups;
  idx_t *xadj, *adjncy, *adjwgt, *vwgt;
  idx_t *part, *gwhere0, *gwhere1;
  idx_t *tmpwhere, *tmpvwgt, *tmpxadj, *tmpadjncy, *tmpadjwgt;
  graph_t *agraph;
  idx_t lnparts, fpart, fpe, lnpes; 
  idx_t twoparts=2, moptions[METIS_NOPTIONS], edgecut, max_cut;
  real_t *tpwgts, *tpwgts2, *lbvec, lbsum, min_lbsum, wsum;
  MPI_Comm ipcomm;
  struct {
    double sum;
    int rank;
  } lpesum, gpesum;

  WCOREPUSH;

  ncon = graph->ncon;

  ngroups = gk_max(gk_min(RIP_SPLIT_FACTOR, ctrl->npes), 1);

  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));

  lbvec = rwspacemalloc(ctrl, ncon);

  /* assemble the graph to all the processors */
  agraph = AssembleAdaptiveGraph(ctrl, graph);
  gnvtxs = agraph->nvtxs;

  /* make a copy of the graph's structure for later */
  xadj   = icopy(gnvtxs+1, agraph->xadj, iwspacemalloc(ctrl, gnvtxs+1));
  vwgt   = icopy(gnvtxs*ncon, agraph->vwgt, iwspacemalloc(ctrl, gnvtxs*ncon));
  adjncy = icopy(agraph->nedges, agraph->adjncy, iwspacemalloc(ctrl, agraph->nedges));
  adjwgt = icopy(agraph->nedges, agraph->adjwgt, iwspacemalloc(ctrl, agraph->nedges));
  part   = iwspacemalloc(ctrl, gnvtxs);

  /* create different processor groups */
  gkMPI_Comm_split(ctrl->gcomm, ctrl->mype % ngroups, 0, &ipcomm);
  gkMPI_Comm_rank(ipcomm, &mype);
  gkMPI_Comm_size(ipcomm, &npes);


  /* Go into the recursive bisection */
  METIS_SetDefaultOptions(moptions);
  moptions[METIS_OPTION_SEED] = ctrl->sync + (ctrl->mype % ngroups) + 1;

  tpwgts  = ctrl->tpwgts;
  tpwgts2 = rwspacemalloc(ctrl, 2*ncon);

  lnparts = ctrl->nparts;
  fpart = fpe = 0;
  lnpes = npes;
  while (lnpes > 1 && lnparts > 1) {
    /* determine the weights of the two partitions as a function of the 
       weight of the target partition weights */
    for (j=(lnparts>>1), i=0; i<ncon; i++) {
      tpwgts2[i]      = rsum(j, tpwgts+fpart*ncon+i, ncon);
      tpwgts2[ncon+i] = rsum(lnparts-j, tpwgts+(fpart+j)*ncon+i, ncon);
      wsum            = 1.0/(tpwgts2[i] + tpwgts2[ncon+i]);
      tpwgts2[i]      *= wsum;
      tpwgts2[ncon+i] *= wsum;
    }

    METIS_PartGraphRecursive(&agraph->nvtxs, &ncon, agraph->xadj, agraph->adjncy, 
          agraph->vwgt, NULL, agraph->adjwgt, &twoparts, tpwgts2, NULL, moptions, 
          &edgecut, part);

    /* pick one of the branches */
    if (mype < fpe+lnpes/2) {
      KeepPart(ctrl, agraph, part, 0);
      lnpes   = lnpes/2;
      lnparts = lnparts/2;
    }
    else {
      KeepPart(ctrl, agraph, part, 1);
      fpart   = fpart + lnparts/2;
      fpe     = fpe + lnpes/2;
      lnpes   = lnpes - lnpes/2;
      lnparts = lnparts - lnparts/2;
    }
  }

  gwhere0 = iset(gnvtxs, 0, iwspacemalloc(ctrl, gnvtxs));
  gwhere1 = iwspacemalloc(ctrl, gnvtxs);

  if (lnparts == 1) { /* Case npes is greater than or equal to nparts */
    /* Only the first process will assign labels (for the reduction to work) */
    if (mype == fpe) {
      for (i=0; i<agraph->nvtxs; i++) 
        gwhere0[agraph->label[i]] = fpart;
    }
  }
  else { /* Case in which npes is smaller than nparts */
    /* create the normalized tpwgts for the lnparts from ctrl->tpwgts */
    tpwgts = rwspacemalloc(ctrl, lnparts*ncon);
    for (j=0; j<ncon; j++) {
      for (wsum=0.0, i=0; i<lnparts; i++) {
        tpwgts[i*ncon+j] = ctrl->tpwgts[(fpart+i)*ncon+j];
        wsum += tpwgts[i*ncon+j];
      }
      for (wsum=1.0/wsum, i=0; i<lnparts; i++) 
        tpwgts[i*ncon+j] *= wsum;
    }

    METIS_PartGraphKway(&agraph->nvtxs, &ncon, agraph->xadj, agraph->adjncy, 
          agraph->vwgt, NULL, agraph->adjwgt, &lnparts, tpwgts, NULL, moptions, 
          &edgecut, part);

    for (i=0; i<agraph->nvtxs; i++) 
      gwhere0[agraph->label[i]] = fpart + part[i];
  }

  gkMPI_Allreduce((void *)gwhere0, (void *)gwhere1, gnvtxs, IDX_T, MPI_SUM, ipcomm);

  if (ngroups > 1) {
    tmpxadj   = agraph->xadj;
    tmpadjncy = agraph->adjncy;
    tmpadjwgt = agraph->adjwgt;
    tmpvwgt   = agraph->vwgt;
    tmpwhere  = agraph->where;

    agraph->xadj   = xadj;
    agraph->adjncy = adjncy;
    agraph->adjwgt = adjwgt;
    agraph->vwgt   = vwgt;
    agraph->where  = gwhere1;
    agraph->vwgt   = vwgt;
    agraph->nvtxs  = gnvtxs;

    edgecut = ComputeSerialEdgeCut(agraph);
    ComputeSerialBalance(ctrl, agraph, gwhere1, lbvec);
    lbsum = rsum(ncon, lbvec, 1);

    gkMPI_Allreduce((void *)&edgecut, (void *)&max_cut,   1, IDX_T,  MPI_MAX, ctrl->gcomm);
    gkMPI_Allreduce((void *)&lbsum,   (void *)&min_lbsum, 1, REAL_T, MPI_MIN, ctrl->gcomm);

    lpesum.sum = lbsum;
    if (min_lbsum < UNBALANCE_FRACTION*ncon) {
      if (lbsum < UNBALANCE_FRACTION*ncon)
        lpesum.sum = edgecut;
      else
        lpesum.sum = max_cut;
    } 
    lpesum.rank = ctrl->mype;
    
    gkMPI_Allreduce((void *)&lpesum, (void *)&gpesum, 1, MPI_DOUBLE_INT,
        MPI_MINLOC, ctrl->gcomm);
    gkMPI_Bcast((void *)gwhere1, gnvtxs, IDX_T, gpesum.rank, ctrl->gcomm);

    agraph->xadj   = tmpxadj;
    agraph->adjncy = tmpadjncy;
    agraph->adjwgt = tmpadjwgt;
    agraph->vwgt   = tmpvwgt;
    agraph->where  = tmpwhere;
  }

  icopy(graph->nvtxs, gwhere1+graph->vtxdist[ctrl->mype], graph->where);

  FreeGraph(agraph);
  gkMPI_Comm_free(&ipcomm);

  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));

  WCOREPOP;
}


/*************************************************************************
* This function keeps one parts
**************************************************************************/
void KeepPart(ctrl_t *ctrl, graph_t *graph, idx_t *part, idx_t mypart)
{
  idx_t h, i, j, k;
  idx_t nvtxs, ncon, mynvtxs, mynedges;
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *label;
  idx_t *rename;

  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  label  = graph->label;

  rename = iwspacemalloc(ctrl, nvtxs);
 
  for (mynvtxs=0, i=0; i<nvtxs; i++) {
    if (part[i] == mypart)
      rename[i] = mynvtxs++;
  }

  for (mynvtxs=0, mynedges=0, j=xadj[0], i=0; i<nvtxs; i++) {
    if (part[i] == mypart) {
      for (; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (part[k] == mypart) {
          adjncy[mynedges] = rename[k];
          adjwgt[mynedges++] = adjwgt[j];
        }
      }
      j = xadj[i+1];  /* Save xadj[i+1] for later use */

      for (h=0; h<ncon; h++)
        vwgt[mynvtxs*ncon+h] = vwgt[i*ncon+h];

      label[mynvtxs] = label[i];
      xadj[++mynvtxs] = mynedges;
    }
    else {
      j = xadj[i+1];  /* Save xadj[i+1] for later use */
    }
  }

  graph->nvtxs  = mynvtxs;
  graph->nedges = mynedges;

  WCOREPOP;
}


