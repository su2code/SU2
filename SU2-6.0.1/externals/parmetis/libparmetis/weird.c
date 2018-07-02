/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * weird.c
 *
 * This file contain various graph setting up routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: weird.c 10592 2011-07-16 21:17:53Z karypis $
 *
 */

#include <parmetislib.h>

#define RIFN(x) \
  if ((x) == NULL) {\
    printf("PARMETIS ERROR " #x " is NULL.\n");\
    return 0;\
  }

#define RIFNP(x) \
  if ((*x) <= 0) {\
    printf("PARMETIS ERROR " #x " is <= 0.\n");\
    return 0;\
  }


/*************************************************************************/
/*! This function checks the validity of the inputs for PartKway
  */
/*************************************************************************/
int CheckInputsPartKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
        idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts,
        real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part,
        MPI_Comm *comm)
{
  idx_t i, j, mype;
  real_t sum;

  /* Check that the supplied information is actually non-NULL */
  if (comm == NULL) {
    printf("PARMETIS ERROR: comm is NULL. Aborting\n");
    abort();
  }
  gkMPI_Comm_rank(*comm, &mype);

  RIFN(vtxdist);
  RIFN(xadj);
  RIFN(adjncy);
  RIFN(wgtflag);
  RIFN(numflag);
  RIFN(ncon);
  RIFN(nparts);
  RIFN(tpwgts);
  RIFN(ubvec);
  RIFN(options);
  RIFN(edgecut);
  RIFN(part);

  if (*wgtflag == 2 || *wgtflag == 3) {
    RIFN(vwgt);
    for (j=0; j<*ncon; j++) {
      if (GlobalSESumComm(*comm, isum(vtxdist[mype+1]-vtxdist[mype], vwgt+j, *ncon)) == 0) {
        printf("PARMETIS ERROR: sum weight for constraint %"PRIDX" is zero.\n", j);
        return 0;
      }
    }
  }
  if (*wgtflag == 1 || *wgtflag == 3) 
    RIFN(adjwgt);


  /* Check that the supplied information is actually valid/reasonable */
  if (vtxdist[mype+1]-vtxdist[mype] < 1) {
    printf("PARMETIS ERROR: Poor initial vertex distribution. "
           "Processor %"PRIDX" has no vertices assigned to it!\n", mype);
    return 0;
  }

  RIFNP(ncon);
  RIFNP(nparts);


  for (j=0; j<*ncon; j++) {
    sum = rsum(*nparts, tpwgts+j, *ncon);
    if (sum < 0.999 || sum > 1.001) {
      printf("PARMETIS ERROR: The sum of tpwgts for constraint #%"PRIDX" is not 1.0\n", j);
      return 0;
    }
  }
  for (j=0; j<*ncon; j++) {
    for (i=0; i<*nparts; i++) {
      if (tpwgts[i*(*ncon)+j] < 0.0 || tpwgts[i] > 1.001) {
        printf("PARMETIS ERROR: The tpwgts for constraint #%"PRIDX" and partition #%"PRIDX" is out of bounds.\n", j, i);
        return 0;
      }
    }
  }


  for (j=0; j<*ncon; j++) {
    if (ubvec[j] <= 1.0) {
      printf("PARMETIS ERROR: The ubvec for constraint #%"PRIDX" must be > 1.0\n", j);
      return 0;
    }
  }

  return 1;
}


/*************************************************************************/
/*! This function checks the validity of the inputs for PartGeomKway
  */
/*************************************************************************/
int CheckInputsPartGeomKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
        idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ndims, real_t *xyz, 
        idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
        idx_t *edgecut, idx_t *part, MPI_Comm *comm)
{
  idx_t i, j, mype;
  real_t sum;

  /* Check that the supplied information is actually non-NULL */
  if (comm == NULL) {
    printf("PARMETIS ERROR: comm is NULL. Aborting\n");
    abort();
  }
  gkMPI_Comm_rank(*comm, &mype);

  RIFN(vtxdist);
  RIFN(xadj);
  RIFN(adjncy);
  RIFN(xyz);
  RIFN(ndims);
  RIFN(wgtflag);
  RIFN(numflag);
  RIFN(ncon);
  RIFN(nparts);
  RIFN(tpwgts);
  RIFN(ubvec);
  RIFN(options);
  RIFN(edgecut);
  RIFN(part);

  if (*wgtflag == 2 || *wgtflag == 3) {
    RIFN(vwgt);
    for (j=0; j<*ncon; j++) {
      if (GlobalSESumComm(*comm, isum(vtxdist[mype+1]-vtxdist[mype], vwgt+j, *ncon)) == 0) {
        printf("PARMETIS ERROR: sum weight for constraint %"PRIDX" is zero.\n", j);
        return 0;
      }
    }
  }
  if (*wgtflag == 1 || *wgtflag == 3)
    RIFN(adjwgt);


  /* Check that the supplied information is actually valid/reasonable */
  if (vtxdist[mype+1]-vtxdist[mype] < 1) {
    printf("PARMETIS ERROR: Poor initial vertex distribution. "
           "Processor %"PRIDX" has no vertices assigned to it!\n", mype);
    return 0;
  }

  RIFNP(ncon);
  RIFNP(nparts);
  RIFNP(ndims);

  if (*ndims > 3) {
    printf("PARMETIS ERROR: The ndims should be <= 3.\n");
    return 0;
  }

  for (j=0; j<*ncon; j++) {
    sum = rsum(*nparts, tpwgts+j, *ncon);
    if (sum < 0.999 || sum > 1.001) {
      printf("PARMETIS ERROR: The sum of tpwgts for constraint #%"PRIDX" is not 1.0\n", j);
      return 0;
    }
  }
  for (j=0; j<*ncon; j++) {
    for (i=0; i<*nparts; i++) {
      if (tpwgts[i*(*ncon)+j] < 0.0 || tpwgts[i] > 1.001) {
        printf("PARMETIS ERROR: The tpwgts for constraint #%"PRIDX" and partition #%"PRIDX" is out of bounds.\n", j, i);
        return 0;
      }
    }
  }


  for (j=0; j<*ncon; j++) {
    if (ubvec[j] <= 1.0) {
      printf("PARMETIS ERROR: The ubvec for constraint #%"PRIDX" must be > 1.0\n", j);
      return 0;
    }
  }

  return 1;
}


/*************************************************************************/
/*! This function checks the validity of the inputs for PartGeom
  */
/*************************************************************************/
int CheckInputsPartGeom(idx_t *vtxdist, idx_t *ndims, real_t *xyz, 
        idx_t *part, MPI_Comm *comm)
{
  idx_t mype;

  /* Check that the supplied information is actually non-NULL */
  if (comm == NULL) {
    printf("PARMETIS ERROR: comm is NULL. Aborting\n");
    abort();
  }

  RIFN(vtxdist);
  RIFN(xyz);
  RIFN(ndims);
  RIFN(part);

  /* Check that the supplied information is actually valid/reasonable */
  gkMPI_Comm_rank(*comm, &mype);
  if (vtxdist[mype+1]-vtxdist[mype] < 1) {
    printf("PARMETIS ERROR: Poor initial vertex distribution. "
           "Processor %"PRIDX" has no vertices assigned to it!\n", mype);
    return 0;
  }

  RIFNP(ndims);

  if (*ndims > 3) {
    printf("PARMETIS ERROR: The ndims should be <= 3.\n");
    return 0;
  }

  return 1;
}


/*************************************************************************/
/*! This function checks the validity of the inputs for AdaptiveRepart
  */
/*************************************************************************/
int CheckInputsAdaptiveRepart(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, 
        idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *wgtflag, 
        idx_t *numflag, idx_t *ncon, idx_t *nparts, real_t *tpwgts, 
        real_t *ubvec, real_t *ipc2redist, idx_t *options, idx_t *edgecut, 
        idx_t *part, MPI_Comm *comm)
{
  idx_t i, j, mype;
  real_t sum;

  /* Check that the supplied information is actually non-NULL */
  if (comm == NULL) {
    printf("PARMETIS ERROR: comm is NULL. Aborting\n");
    abort();
  }
  gkMPI_Comm_rank(*comm, &mype);

  RIFN(vtxdist);
  RIFN(xadj);
  RIFN(adjncy);
  /*RIFN(vsize);*/
  RIFN(wgtflag);
  RIFN(numflag);
  RIFN(ncon);
  RIFN(nparts);
  RIFN(tpwgts);
  RIFN(ubvec);
  RIFN(options);
  RIFN(edgecut);
  RIFN(part);

  if (*wgtflag == 2 || *wgtflag == 3) {
    RIFN(vwgt);
    for (j=0; j<*ncon; j++) {
      if (GlobalSESumComm(*comm, isum(vtxdist[mype+1]-vtxdist[mype], vwgt+j, *ncon)) == 0) {
        printf("PARMETIS ERROR: sum weight for constraint %"PRIDX" is zero.\n", j);
        return 0;
      }
    }
  }
  if (*wgtflag == 1 || *wgtflag == 3)
    RIFN(adjwgt);


  /* Check that the supplied information is actually valid/reasonable */
  if (vtxdist[mype+1]-vtxdist[mype] < 1) {
    printf("PARMETIS ERROR: Poor initial vertex distribution. "
           "Processor %"PRIDX" has no vertices assigned to it!\n", mype);
    return 0;
  }

  RIFNP(ncon);
  RIFNP(nparts);


  for (j=0; j<*ncon; j++) {
    sum = rsum(*nparts, tpwgts+j, *ncon);
    if (sum < 0.999 || sum > 1.001) {
      printf("PARMETIS ERROR: The sum of tpwgts for constraint #%"PRIDX" is not 1.0\n", j);
      return 0;
    }
  }
  for (j=0; j<*ncon; j++) {
    for (i=0; i<*nparts; i++) {
      if (tpwgts[i*(*ncon)+j] < 0.0 || tpwgts[i] > 1.001) {
        printf("PARMETIS ERROR: The tpwgts for constraint #%"PRIDX" and partition #%"PRIDX" is out of bounds.\n", j, i);
        return 0;
      }
    }
  }


  for (j=0; j<*ncon; j++) {
    if (ubvec[j] <= 1.0) {
      printf("PARMETIS ERROR: The ubvec for constraint #%"PRIDX" must be > 1.0\n", j);
      return 0;
    }
  }

  if (*ipc2redist < 0.0001 || *ipc2redist > 1000000.0) {
    printf("PARMETIS ERROR: The ipc2redist value should be between [0.0001, 1000000.0]\n");
    return 0;
  }

  return 1;
}


/*************************************************************************/
/*! This function checks the validity of the inputs for NodeND
  */
/*************************************************************************/
int CheckInputsNodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, 
        idx_t *numflag, idx_t *options, idx_t *order, idx_t *sizes,
        MPI_Comm *comm)
{
  idx_t mype;

  /* Check that the supplied information is actually non-NULL */
  if (comm == NULL) {
    printf("PARMETIS ERROR: comm is NULL. Aborting\n");
    abort();
  }

  RIFN(vtxdist);
  RIFN(xadj);
  RIFN(adjncy);
  RIFN(numflag);
  RIFN(options);
  RIFN(order);
  RIFN(sizes);

  /* Check that the supplied information is actually valid/reasonable */
  gkMPI_Comm_rank(*comm, &mype);
  if (vtxdist[mype+1]-vtxdist[mype] < 1) {
    printf("PARMETIS ERROR: Poor initial vertex distribution. "
           "Processor %"PRIDX" has no vertices assigned to it!\n", mype);
    return 0;
  }

  return 1;
}



/*************************************************************************/
/*! This function checks the validity of the inputs for PartMeshKway
  */
/*************************************************************************/
int CheckInputsPartMeshKway(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt,
        idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *ncommon, idx_t *nparts,
        real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part,
        MPI_Comm *comm)
{
  idx_t i, j, mype;
  real_t sum;

  /* Check that the supplied information is actually non-NULL */
  if (comm == NULL) {
    printf("PARMETIS ERROR: comm is NULL. Aborting\n");
    abort();
  }

  RIFN(elmdist);
  RIFN(eptr);
  RIFN(eind);
  RIFN(wgtflag);
  RIFN(numflag);
  RIFN(ncon);
  RIFN(nparts);
  RIFN(tpwgts);
  RIFN(ubvec);
  RIFN(options);
  RIFN(edgecut);
  RIFN(part);

  if (*wgtflag == 2 || *wgtflag == 3)
    RIFN(elmwgt);


  /* Check that the supplied information is actually valid/reasonable */
  gkMPI_Comm_rank(*comm, &mype);
  if (elmdist[mype+1]-elmdist[mype] < 1) {
    printf("PARMETIS ERROR: Poor initial element distribution. "
           "Processor %"PRIDX" has no elements assigned to it!\n", mype);
    return 0;
  }

  RIFNP(ncon);
  RIFNP(nparts);


  for (j=0; j<*ncon; j++) {
    sum = rsum(*nparts, tpwgts+j, *ncon);
    if (sum < 0.999 || sum > 1.001) {
      printf("PARMETIS ERROR: The sum of tpwgts for constraint #%"PRIDX" is not 1.0\n", j);
      return 0;
    }
  }
  for (j=0; j<*ncon; j++) {
    for (i=0; i<*nparts; i++) {
      if (tpwgts[i*(*ncon)+j] < 0.0 || tpwgts[i] > 1.001) {
        printf("PARMETIS ERROR: The tpwgts for constraint #%"PRIDX" and partition #%"PRIDX" is out of bounds.\n", j, i);
        return 0;
      }
    }
  }


  for (j=0; j<*ncon; j++) {
    if (ubvec[j] <= 1.0) {
      printf("PARMETIS ERROR: The ubvec for constraint #%"PRIDX" must be > 1.0\n", j);
      return 0;
    }
  }

  return 1;
}



/*************************************************************************/
/*! This function computes a partitioning of a small graph
  */
/*************************************************************************/
void PartitionSmallGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, h, ncon, nparts, npes, mype;
  idx_t moptions[METIS_NOPTIONS];
  idx_t me;
  idx_t *mypart;
  int lpecut[2], gpecut[2];
  graph_t *agraph;
  idx_t *sendcounts, *displs;
  real_t *gnpwgts, *lnpwgts;

  ncon   = graph->ncon;
  nparts = ctrl->nparts;
  npes   = ctrl->npes;
  mype   = ctrl->mype;

  WCOREPUSH;

  CommSetup(ctrl, graph);
  graph->where = imalloc(graph->nvtxs+graph->nrecv, "PartitionSmallGraph: where");
  agraph       = AssembleAdaptiveGraph(ctrl, graph);
  mypart       = iwspacemalloc(ctrl, agraph->nvtxs);

  METIS_SetDefaultOptions(moptions);
  moptions[METIS_OPTION_SEED] = ctrl->sync + mype;

  METIS_PartGraphKway(&agraph->nvtxs, &ncon, agraph->xadj, agraph->adjncy, 
        agraph->vwgt, NULL, agraph->adjwgt, &nparts, ctrl->tpwgts, NULL,
	moptions, &graph->mincut, mypart);

  lpecut[0] = graph->mincut;
  lpecut[1] = mype;
  gkMPI_Allreduce(lpecut, gpecut, 1, MPI_2INT, MPI_MINLOC, ctrl->comm);
  graph->mincut = gpecut[0];

  if (lpecut[1] == gpecut[1] && gpecut[1] != 0)
    gkMPI_Send((void *)mypart, agraph->nvtxs, IDX_T, 0, 1, ctrl->comm);
  if (lpecut[1] == 0 && gpecut[1] != 0)
    gkMPI_Recv((void *)mypart, agraph->nvtxs, IDX_T, gpecut[1], 1, ctrl->comm, &ctrl->status);

  sendcounts = iwspacemalloc(ctrl, npes);
  displs     = iwspacemalloc(ctrl, npes);

  for (i=0; i<npes; i++) {
    sendcounts[i] = graph->vtxdist[i+1]-graph->vtxdist[i];
    displs[i] = graph->vtxdist[i];
  }

  gkMPI_Scatterv((void *)mypart, sendcounts, displs, IDX_T,
      (void *)graph->where, graph->nvtxs, IDX_T, 0, ctrl->comm);

  lnpwgts = graph->lnpwgts = rmalloc(nparts*ncon, "lnpwgts");
  gnpwgts = graph->gnpwgts = rmalloc(nparts*ncon, "gnpwgts");
  rset(nparts*ncon, 0, lnpwgts);
  for (i=0; i<graph->nvtxs; i++) {
    me = graph->where[i];
    for (h=0; h<ncon; h++)
      lnpwgts[me*ncon+h] += graph->nvwgt[i*ncon+h];
  }
  gkMPI_Allreduce((void *)lnpwgts, (void *)gnpwgts, nparts*ncon, REAL_T, MPI_SUM, ctrl->comm);

  FreeGraph(agraph);

  WCOREPOP;

  return;
}

