/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * initbalance.c
 *
 * This file contains code that computes an initial partitioning
 *
 * Started 3/4/96
 * George
 *
 * $Id: initbalance.c 10592 2011-07-16 21:17:53Z karypis $
 */

#include <parmetislib.h>


/*************************************************************************
* This function is the entry point of the initial balancing algorithm.
* This algorithm assembles the graph to all the processors and preceeds
* with the balancing step.
**************************************************************************/
void Balance_Partition(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, nvtxs, nedges, ncon;
  idx_t mype, npes, srnpes, srmype; 
  idx_t *vtxdist, *xadj, *adjncy, *adjwgt, *vwgt, *vsize;
  idx_t *part, *lwhere, *home;
  idx_t lnparts, fpart, fpe, lnpes, ngroups;
  idx_t *rcounts, *rdispls;
  idx_t twoparts=2, moptions[METIS_NOPTIONS], edgecut, max_cut;
  idx_t sr_pe, gd_pe, sr, gd, who_wins;
  real_t my_cut, my_totalv, my_cost = -1.0, my_balance = -1.0, wsum;
  real_t rating, max_rating, your_cost = -1.0, your_balance = -1.0;
  real_t lbsum, min_lbsum, *lbvec, *tpwgts, *tpwgts2, buffer[2];
  graph_t *agraph, cgraph;
  ctrl_t *myctrl;
  MPI_Status status;
  MPI_Comm ipcomm, srcomm;
  struct {
    double cost;
    int rank;
  } lpecost, gpecost;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));
  WCOREPUSH;

  vtxdist = graph->vtxdist;
  agraph  = AssembleAdaptiveGraph(ctrl, graph);
  nvtxs   = cgraph.nvtxs  = agraph->nvtxs;
  nedges  = cgraph.nedges = agraph->nedges;
  ncon    = cgraph.ncon   = agraph->ncon;
  xadj    = cgraph.xadj   = icopy(nvtxs+1, agraph->xadj, iwspacemalloc(ctrl, nvtxs+1));
  vwgt    = cgraph.vwgt   = icopy(nvtxs*ncon, agraph->vwgt, iwspacemalloc(ctrl, nvtxs*ncon));
  vsize   = cgraph.vsize  = icopy(nvtxs, agraph->vsize, iwspacemalloc(ctrl, nvtxs));
  adjncy  = cgraph.adjncy = icopy(nedges, agraph->adjncy, iwspacemalloc(ctrl, nedges));
  adjwgt  = cgraph.adjwgt = icopy(nedges, agraph->adjwgt, iwspacemalloc(ctrl, nedges));
  part    = cgraph.where  = agraph->where = iwspacemalloc(ctrl, nvtxs);

  lwhere = iwspacemalloc(ctrl, nvtxs);
  home   = iwspacemalloc(ctrl, nvtxs);
  lbvec  = rwspacemalloc(ctrl, graph->ncon);


  /****************************************/
  /****************************************/
  if (ctrl->ps_relation == PARMETIS_PSR_UNCOUPLED) {
    WCOREPUSH;
    rcounts = iwspacemalloc(ctrl, ctrl->npes);
    rdispls = iwspacemalloc(ctrl, ctrl->npes+1);

    for (i=0; i<ctrl->npes; i++) 
      rdispls[i] = rcounts[i] = vtxdist[i+1]-vtxdist[i];
    MAKECSR(i, ctrl->npes, rdispls);

    gkMPI_Allgatherv((void *)graph->home, graph->nvtxs, IDX_T,
        (void *)part, rcounts, rdispls, IDX_T, ctrl->comm);

    for (i=0; i<agraph->nvtxs; i++)
      home[i] = part[i];

    WCOREPOP;  /* local frees */
  }
  else {
    for (i=0; i<ctrl->npes; i++) {
      for (j=vtxdist[i]; j<vtxdist[i+1]; j++)
        part[j] = home[j] = i;
    }
  }

  /* Ensure that the initial partitioning is legal */
  for (i=0; i<agraph->nvtxs; i++) {
    if (part[i] >= ctrl->nparts)
      part[i] = home[i] = part[i] % ctrl->nparts;
    if (part[i] < 0)
      part[i] = home[i] = (-1*part[i]) % ctrl->nparts;
  }
  /****************************************/
  /****************************************/

  IFSET(ctrl->dbglvl, DBG_REFINEINFO, 
      ComputeSerialBalance(ctrl, agraph, agraph->where, lbvec));
  IFSET(ctrl->dbglvl, DBG_REFINEINFO, 
      rprintf(ctrl, "input cut: %"PRIDX", balance: ", ComputeSerialEdgeCut(agraph)));
  for (i=0; i<agraph->ncon; i++)
    IFSET(ctrl->dbglvl, DBG_REFINEINFO, rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]));
  IFSET(ctrl->dbglvl, DBG_REFINEINFO, rprintf(ctrl, "\n"));

  /****************************************/
  /* Split the processors into two groups */
  /****************************************/
  sr = (ctrl->mype % 2 == 0) ? 1 : 0;
  gd = (ctrl->mype % 2 == 1) ? 1 : 0;

  if (graph->ncon > MAX_NCON_FOR_DIFFUSION || ctrl->npes == 1) {
    sr = 1;
    gd = 0;
  }

  sr_pe = 0;
  gd_pe = 1;

  gkMPI_Comm_split(ctrl->gcomm, sr, 0, &ipcomm);
  gkMPI_Comm_rank(ipcomm, &mype);
  gkMPI_Comm_size(ipcomm, &npes);

  if (sr == 1) { /* Half of the processors do scratch-remap */
    ngroups = gk_max(gk_min(RIP_SPLIT_FACTOR, npes), 1);
    gkMPI_Comm_split(ipcomm, mype % ngroups, 0, &srcomm);
    gkMPI_Comm_rank(srcomm, &srmype);
    gkMPI_Comm_size(srcomm, &srnpes);

    METIS_SetDefaultOptions(moptions);
    moptions[METIS_OPTION_SEED] = ctrl->sync + (mype % ngroups) + 1;

    tpwgts  = ctrl->tpwgts;
    tpwgts2 = rwspacemalloc(ctrl, 2*ncon);

    iset(nvtxs, 0, lwhere);
    lnparts = ctrl->nparts;
    fpart = fpe = 0;
    lnpes = srnpes;
    while (lnpes > 1 && lnparts > 1) {
      PASSERT(ctrl, agraph->nvtxs > 1);
      /* determine the weights of the two partitions as a function of the 
         weight of the target partition weights */
      for (j=(lnparts>>1), i=0; i<ncon; i++) {
        tpwgts2[i]      = rsum(j, tpwgts+fpart*ncon+i, ncon);
        tpwgts2[ncon+i] = rsum(lnparts-j, tpwgts+(fpart+j)*ncon+i, ncon);
        wsum            = 1.0/(tpwgts2[i] + tpwgts2[ncon+i]);
        tpwgts2[i]      *= wsum;
        tpwgts2[ncon+i] *= wsum;
      }

      METIS_PartGraphRecursive(&agraph->nvtxs, &ncon, agraph->xadj, 
            agraph->adjncy, agraph->vwgt, NULL, agraph->adjwgt, 
            &twoparts, tpwgts2, NULL, moptions, &edgecut, part);

      /* pick one of the branches */
      if (srmype < fpe+lnpes/2) {
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

    if (lnparts == 1) { /* Case in which srnpes is greater or equal to nparts */
      /* Only the first process will assign labels (for the reduction to work) */
      if (srmype == fpe) {
        for (i=0; i<agraph->nvtxs; i++) 
          lwhere[agraph->label[i]] = fpart;
      }
    }
    else { /* Case in which srnpes is smaller than nparts */
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
        lwhere[agraph->label[i]] = fpart + part[i];
    }

    gkMPI_Allreduce((void *)lwhere, (void *)part, nvtxs, IDX_T, MPI_SUM, srcomm);

    edgecut = ComputeSerialEdgeCut(&cgraph);
    ComputeSerialBalance(ctrl, &cgraph, part, lbvec);
    lbsum = rsum(ncon, lbvec, 1);
    gkMPI_Allreduce((void *)&edgecut, (void *)&max_cut, 1, IDX_T, MPI_MAX, ipcomm);
    gkMPI_Allreduce((void *)&lbsum, (void *)&min_lbsum, 1, REAL_T, MPI_MIN, ipcomm);
    lpecost.rank = ctrl->mype;
    lpecost.cost = lbsum;
    if (min_lbsum < UNBALANCE_FRACTION * (real_t)(ncon)) {
      if (lbsum < UNBALANCE_FRACTION * (real_t)(ncon))
        lpecost.cost = (double)edgecut;
      else
        lpecost.cost = (double)max_cut + lbsum;
    }
    gkMPI_Allreduce((void *)&lpecost, (void *)&gpecost, 1, MPI_DOUBLE_INT,
        MPI_MINLOC, ipcomm);

    if (ctrl->mype == gpecost.rank && ctrl->mype != sr_pe) 
      gkMPI_Send((void *)part, nvtxs, IDX_T, sr_pe, 1, ctrl->comm);

    if (ctrl->mype != gpecost.rank && ctrl->mype == sr_pe) 
      gkMPI_Recv((void *)part, nvtxs, IDX_T, gpecost.rank, 1, ctrl->comm, &status);

    if (ctrl->mype == sr_pe) {
      icopy(nvtxs, part, lwhere);
      SerialRemap(ctrl, &cgraph, ctrl->nparts, home, lwhere, part, ctrl->tpwgts);
    }

    gkMPI_Comm_free(&srcomm);
  }
  else { /* The other half do global diffusion */
    /* setup a ctrl for the diffusion */
    myctrl = (ctrl_t *)gk_malloc(sizeof(ctrl_t), "myctrl");
    memset(myctrl, 0, sizeof(ctrl_t));
    myctrl->mype          = mype;
    myctrl->npes          = npes;
    myctrl->comm          = ipcomm;
    myctrl->sync          = ctrl->sync;
    myctrl->seed          = ctrl->seed;
    myctrl->nparts        = ctrl->nparts;
    myctrl->ncon          = ctrl->ncon;
    myctrl->ipc_factor    = ctrl->ipc_factor;
    myctrl->redist_factor = ctrl->redist_base;
    myctrl->partType      = ADAPTIVE_PARTITION;
    myctrl->ps_relation   = PARMETIS_PSR_UNCOUPLED;
    myctrl->tpwgts        = rmalloc(myctrl->nparts*myctrl->ncon, "myctrl->tpwgts");
    myctrl->ubvec         = rmalloc(myctrl->ncon, "myctrl->ubvec");
    myctrl->invtvwgts     = rmalloc(myctrl->ncon, "myctrl->invtvwgts");

    rcopy(myctrl->nparts*myctrl->ncon, ctrl->tpwgts, myctrl->tpwgts);
    rcopy(myctrl->ncon, ctrl->ubvec, myctrl->ubvec);
    rcopy(myctrl->ncon, ctrl->invtvwgts, myctrl->invtvwgts);

    AllocateWSpace(myctrl, 10*agraph->nvtxs);
    AllocateRefinementWorkSpace(myctrl, agraph->nvtxs);

    /* This stmt is required to balance out the sr gkMPI_Comm_split */
    gkMPI_Comm_split(ipcomm, MPI_UNDEFINED, 0, &srcomm);

    if (ncon == 1) {
      rating = WavefrontDiffusion(myctrl, agraph, home);
      ComputeSerialBalance(ctrl, &cgraph, part, lbvec);
      lbsum = rsum(ncon, lbvec, 1);

      /* Determine which PE computed the best partitioning */
      gkMPI_Allreduce((void *)&rating, (void *)&max_rating, 1, REAL_T, MPI_MAX, ipcomm);
      gkMPI_Allreduce((void *)&lbsum, (void *)&min_lbsum, 1, REAL_T, MPI_MIN, ipcomm);

      lpecost.rank = ctrl->mype;
      lpecost.cost = lbsum;
      if (min_lbsum < UNBALANCE_FRACTION * (real_t)(ncon)) {
        if (lbsum < UNBALANCE_FRACTION * (real_t)(ncon))
          lpecost.cost = rating;
        else
          lpecost.cost = max_rating + lbsum;
      }

      gkMPI_Allreduce((void *)&lpecost, (void *)&gpecost, 1, MPI_DOUBLE_INT, 
          MPI_MINLOC, ipcomm);

      /* Now send this to the coordinating processor */
      if (ctrl->mype == gpecost.rank && ctrl->mype != gd_pe)
        gkMPI_Send((void *)part, nvtxs, IDX_T, gd_pe, 1, ctrl->comm);

      if (ctrl->mype != gpecost.rank && ctrl->mype == gd_pe)
        gkMPI_Recv((void *)part, nvtxs, IDX_T, gpecost.rank, 1, ctrl->comm, &status);

      if (ctrl->mype == gd_pe) {
        icopy(nvtxs, part, lwhere);
        SerialRemap(ctrl, &cgraph, ctrl->nparts, home, lwhere, part, ctrl->tpwgts);
      }
    }
    else {
      Mc_Diffusion(myctrl, agraph, graph->vtxdist, agraph->where, home, 
          N_MOC_GD_PASSES);
    }

    FreeCtrl(&myctrl);
  }


  if (graph->ncon <= MAX_NCON_FOR_DIFFUSION) {
    if (ctrl->mype == sr_pe  || ctrl->mype == gd_pe) {
      /********************************************************************/
      /* The coordinators from each group decide on the best partitioning */
      /********************************************************************/
      my_cut = (real_t) ComputeSerialEdgeCut(&cgraph);
      my_totalv = (real_t) Mc_ComputeSerialTotalV(&cgraph, home);
      ComputeSerialBalance(ctrl, &cgraph, part, lbvec);
      my_balance = rsum(cgraph.ncon, lbvec, 1);
      my_balance /= (real_t) cgraph.ncon;
      my_cost = ctrl->ipc_factor * my_cut + REDIST_WGT * ctrl->redist_base * my_totalv;

      IFSET(ctrl->dbglvl, DBG_REFINEINFO, 
          printf("%s initial cut: %.1"PRREAL", totalv: %.1"PRREAL", balance: %.3"PRREAL"\n",
                 (ctrl->mype == sr_pe ? "scratch-remap" : "diffusion"), 
                 my_cut, my_totalv, my_balance));

      if (ctrl->mype == gd_pe) {
        buffer[0] = my_cost;
        buffer[1] = my_balance;
        gkMPI_Send((void *)buffer, 2, REAL_T, sr_pe, 1, ctrl->comm);
      }
      else {
        gkMPI_Recv((void *)buffer, 2, REAL_T, gd_pe, 1, ctrl->comm, &status);
        your_cost    = buffer[0];
        your_balance = buffer[1];
      }
    }

    if (ctrl->mype == sr_pe) {
      who_wins = gd_pe;
      if ((my_balance < 1.1 && your_balance > 1.1) ||
          (my_balance < 1.1 && your_balance < 1.1 && my_cost < your_cost) ||
          (my_balance > 1.1 && your_balance > 1.1 && my_balance < your_balance)) {
        who_wins = sr_pe;
      }
    }

    gkMPI_Bcast((void *)&who_wins, 1, IDX_T, sr_pe, ctrl->comm);
  }
  else {
    who_wins = sr_pe;
  }

  gkMPI_Bcast((void *)part, nvtxs, IDX_T, who_wins, ctrl->comm);
  icopy(graph->nvtxs, part+vtxdist[ctrl->mype], graph->where);

  gkMPI_Comm_free(&ipcomm);

  agraph->where = NULL;
  FreeGraph(agraph);
          
  WCOREPOP;
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));

}


/*************************************************************************/
/*! This function assembles the graph into a single processor. It should
    work for static, adaptive, single-, and multi-contraint */
/*************************************************************************/
graph_t *AssembleAdaptiveGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, k, l, gnvtxs, nvtxs, ncon, gnedges, nedges, gsize;
  idx_t *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *vtxdist, *imap;
  idx_t *axadj, *aadjncy, *aadjwgt, *avwgt, *avsize = NULL, *alabel;
  idx_t *mygraph, *ggraph;
  idx_t *rcounts, *rdispls, mysize;
  real_t *anvwgt;
  graph_t *agraph;

  WCOREPUSH;

  gnvtxs  = graph->gnvtxs;
  nvtxs   = graph->nvtxs;
  ncon    = graph->ncon;
  nedges  = graph->xadj[nvtxs];
  xadj    = graph->xadj;
  vwgt    = graph->vwgt;
  vsize   = graph->vsize;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  vtxdist = graph->vtxdist;
  imap    = graph->imap;

  /*************************************************************/
  /* Determine the # of idx_t to receive from each processor */
  /*************************************************************/
  rcounts = iwspacemalloc(ctrl, ctrl->npes);
  switch (ctrl->partType) {
    case STATIC_PARTITION:
      mysize = (1+ncon)*nvtxs + 2*nedges;
      break;
    case ADAPTIVE_PARTITION:
    case REFINE_PARTITION:
      mysize = (2+ncon)*nvtxs + 2*nedges;
      break;
    default:
      printf("WARNING: bad value for ctrl->partType %"PRIDX"\n", ctrl->partType);
      break;
  }
  gkMPI_Allgather((void *)(&mysize), 1, IDX_T, (void *)rcounts, 1, IDX_T, ctrl->comm);

  rdispls = iwspacemalloc(ctrl, ctrl->npes+1);
  for (rdispls[0]=0, i=1; i<ctrl->npes+1; i++)
    rdispls[i] = rdispls[i-1] + rcounts[i-1];

  /* allocate memory for the recv buffer of the assembled graph */
  gsize   = rdispls[ctrl->npes];
  ggraph  = iwspacemalloc(ctrl, gsize);

  /* Construct the one-array storage format of the assembled graph */
  WCOREPUSH;  /* for freeing mygraph */
  mygraph = iwspacemalloc(ctrl, mysize);

  for (k=i=0; i<nvtxs; i++) {
    mygraph[k++] = xadj[i+1]-xadj[i];
    for (j=0; j<ncon; j++)
      mygraph[k++] = vwgt[i*ncon+j];
    if (ctrl->partType == ADAPTIVE_PARTITION || ctrl->partType == REFINE_PARTITION)
      mygraph[k++] = vsize[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      mygraph[k++] = imap[adjncy[j]];
      mygraph[k++] = adjwgt[j];
    }
  }
  PASSERT(ctrl, mysize == k);

  /**************************************/
  /* Assemble and send the entire graph */
  /**************************************/
  gkMPI_Allgatherv((void *)mygraph, mysize, IDX_T, (void *)ggraph, 
      rcounts, rdispls, IDX_T, ctrl->comm);

  WCOREPOP;  /* free mygraph */


  agraph = CreateGraph();
  agraph->nvtxs = gnvtxs;
  agraph->ncon  = ncon;

  switch (ctrl->partType) {
    case STATIC_PARTITION:
      agraph->nedges = gnedges = (gsize-(1+ncon)*gnvtxs)/2;
      break;
    case ADAPTIVE_PARTITION:
    case REFINE_PARTITION:
      agraph->nedges = gnedges = (gsize-(2+ncon)*gnvtxs)/2;
      break;
    default:
      printf("WARNING: bad value for ctrl->partType %"PRIDX"\n", ctrl->partType);
      agraph->nedges = gnedges = -1;
      break;
  }


  /*******************************************/
  /* Allocate memory for the assembled graph */
  /*******************************************/
  axadj   = agraph->xadj   = imalloc(gnvtxs+1, "AssembleGraph: axadj");
  avwgt   = agraph->vwgt   = imalloc(gnvtxs*ncon, "AssembleGraph: avwgt");
  anvwgt  = agraph->nvwgt  = rmalloc(gnvtxs*ncon, "AssembleGraph: anvwgt");
  aadjncy = agraph->adjncy = imalloc(gnedges, "AssembleGraph: adjncy");
  aadjwgt = agraph->adjwgt = imalloc(gnedges, "AssembleGraph: adjwgt");
  alabel  = agraph->label  = imalloc(gnvtxs, "AssembleGraph: alabel");
  if (ctrl->partType == ADAPTIVE_PARTITION || ctrl->partType == REFINE_PARTITION)
    avsize = agraph->vsize = imalloc(gnvtxs, "AssembleGraph: avsize");

  for (k=j=i=0; i<gnvtxs; i++) {
    axadj[i] = ggraph[k++];
    for (l=0; l<ncon; l++)
      avwgt[i*ncon+l] = ggraph[k++];
    if (ctrl->partType == ADAPTIVE_PARTITION || ctrl->partType == REFINE_PARTITION)
      avsize[i] = ggraph[k++];
    for (l=0; l<axadj[i]; l++) {
      aadjncy[j] = ggraph[k++];
      aadjwgt[j] = ggraph[k++];
      j++;
    }
  }

  /*********************************/
  /* Now fix up the received graph */
  /*********************************/
  MAKECSR(i, gnvtxs, axadj);

  for (i=0; i<gnvtxs; i++) {
    for (j=0; j<ncon; j++)
      anvwgt[i*ncon+j] = ctrl->invtvwgts[j]*agraph->vwgt[i*ncon+j];
  }

  iincset(gnvtxs, 0, alabel);

  WCOREPOP;

  return agraph;
}


