/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * gkmetis.c
 *
 * This is the entry point of parallel geometry based partitioning
 * routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: gkmetis.c 10663 2011-08-04 03:54:49Z karypis $
 *
 */

#include <parmetislib.h>



/***********************************************************************************
* This function is the entry point of the parallel kmetis algorithm that uses
* coordinates to compute an initial graph distribution.
************************************************************************************/
int ParMETIS_V3_PartGeomKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy,
        idx_t *vwgt, idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ndims, 
	real_t *xyz, idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, 
	idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm)
{
  idx_t h, i, j, npes, mype, status, nvtxs, seed, dbglvl;
  idx_t cut, gcut, maxnvtxs;
  idx_t moptions[METIS_NOPTIONS];
  ctrl_t *ctrl;
  graph_t *graph, *mgraph;
  real_t balance;
  size_t curmem;

  /* Check the input parameters and return if an error */
  status = CheckInputsPartGeomKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag,
                numflag, ndims, xyz, ncon, nparts, tpwgts, ubvec, options, 
                edgecut, part, comm);
  if (GlobalSEMinComm(*comm, status) == 0)
    return METIS_ERROR;

  status = METIS_OK;
  gk_malloc_init();
  curmem = gk_GetCurMemoryUsed();

  /* Setup the ctrl */
  ctrl = SetupCtrl(PARMETIS_OP_GKMETIS, options, *ncon, *nparts, tpwgts, ubvec, *comm);
  npes = ctrl->npes;
  mype = ctrl->mype;

  /* Take care the nparts == 1 case */
  if (*nparts == 1) {
    iset(vtxdist[mype+1]-vtxdist[mype], (*numflag == 0 ? 0 : 1), part);
    *edgecut = 0;
    goto DONE;
  }


  /* Take care of npes == 1 case */
  if (npes == 1) {
    nvtxs = vtxdist[1] - vtxdist[0];  /* subtraction is required when numflag==1 */

    METIS_SetDefaultOptions(moptions);
    moptions[METIS_OPTION_NUMBERING] = *numflag;

    status = METIS_PartGraphKway(&nvtxs, ncon, xadj, adjncy, vwgt, NULL, adjwgt, 
                 nparts, tpwgts, ubvec, moptions, edgecut, part);

    goto DONE;
  }


  /* Setup the graph */
  if (*numflag > 0)
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  graph = SetupGraph(ctrl, *ncon, vtxdist, xadj, vwgt, NULL, adjncy, adjwgt, *wgtflag);
  gk_free((void **)&graph->nvwgt, LTERM); 


  /* Allocate the workspace */
  AllocateWSpace(ctrl, 10*graph->nvtxs);


  /* Compute the initial npes-way partitioning geometric partitioning */
  STARTTIMER(ctrl, ctrl->TotalTmr);

  Coordinate_Partition(ctrl, graph, *ndims, xyz, 1);

  STOPTIMER(ctrl, ctrl->TotalTmr);


  /* Move the graph according to the partitioning */
  STARTTIMER(ctrl, ctrl->MoveTmr);

  ctrl->nparts = npes;
  mgraph = MoveGraph(ctrl, graph);
  ctrl->nparts = *nparts;

  SetupGraph_nvwgts(ctrl, mgraph); /* compute nvwgts for the moved graph */

  if (ctrl->dbglvl&DBG_INFO) {
    CommInterfaceData(ctrl, graph, graph->where, graph->where+graph->nvtxs);
    for (cut=0, i=0; i<graph->nvtxs; i++) {
      for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++) {
        if (graph->where[i] != graph->where[graph->adjncy[j]])
          cut += graph->adjwgt[j];
      }
    }
    gcut     = GlobalSESum(ctrl, cut)/2;
    maxnvtxs = GlobalSEMax(ctrl, mgraph->nvtxs);
    balance  = (real_t)(maxnvtxs)/((real_t)(graph->gnvtxs)/(real_t)(npes));
    rprintf(ctrl, "XYZ Cut: %6"PRIDX" \tBalance: %6.3"PRREAL" [%"PRIDX" %"PRIDX" %"PRIDX"]\n",
       gcut, balance, maxnvtxs, graph->gnvtxs, npes);
  }

  STOPTIMER(ctrl, ctrl->MoveTmr);


  /* Compute the partition of the moved graph */
  STARTTIMER(ctrl, ctrl->TotalTmr);

  ctrl->CoarsenTo = gk_min(vtxdist[npes]+1, 25*(*ncon)*gk_max(npes, *nparts));

  if (vtxdist[npes] < SMALLGRAPH 
      || vtxdist[npes] < npes*20 
      || GlobalSESum(ctrl, mgraph->nedges) == 0) { /* serially */
    IFSET(ctrl->dbglvl, DBG_INFO, 
        rprintf(ctrl, "Partitioning a graph of size %"PRIDX" serially\n", vtxdist[npes]));
    PartitionSmallGraph(ctrl, mgraph);
  }
  else { /* in parallel */
    Global_Partition(ctrl, mgraph);
  }

  ParallelReMapGraph(ctrl, mgraph);

  /* Invert the ordering back to the original graph */
  ctrl->nparts = npes;
  ProjectInfoBack(ctrl, graph, part, mgraph->where);
  ctrl->nparts = *nparts;

  *edgecut = mgraph->mincut;

  STOPTIMER(ctrl, ctrl->TotalTmr);


  /* Print some stats */
  IFSET(ctrl->dbglvl, DBG_TIME, PrintTimingInfo(ctrl));
  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->gcomm));
  IFSET(ctrl->dbglvl, DBG_INFO, PrintPostPartInfo(ctrl, mgraph, 0));

  FreeGraph(mgraph);
  FreeInitialGraphAndRemap(graph);

  if (*numflag > 0)
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

DONE:
  FreeCtrl(&ctrl);
  if (gk_GetCurMemoryUsed() - curmem > 0) {
    printf("ParMETIS appears to have a memory leak of %zdbytes. Report this.\n",
        (ssize_t)(gk_GetCurMemoryUsed() - curmem));
  }
  gk_malloc_cleanup(0);

  return (int)status;
}



/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
int ParMETIS_V3_PartGeom(idx_t *vtxdist, idx_t *ndims, real_t *xyz, idx_t *part, 
         MPI_Comm *comm)
{
  idx_t i, nvtxs, firstvtx, npes, mype, status;
  idx_t *xadj, *adjncy;
  ctrl_t *ctrl=NULL;
  graph_t *graph=NULL;
  size_t curmem;


  /* Check the input parameters and return if an error */
  status = CheckInputsPartGeom(vtxdist, ndims, xyz, part, comm);
  if (GlobalSEMinComm(*comm, status) == 0)
    return METIS_ERROR;

  status = METIS_OK;
  gk_malloc_init();
  curmem = gk_GetCurMemoryUsed();

  /* Setup the ctrl */
  ctrl = SetupCtrl(PARMETIS_OP_GMETIS, NULL, 1, 1, NULL, NULL, *comm);
  /*ctrl->dbglvl=15;*/
  npes = ctrl->npes;
  mype = ctrl->mype;


  /* Trivial case when npes == 1 */
  if (npes == 1) {
    iset(vtxdist[mype+1]-vtxdist[mype], 0, part);
    goto DONE;
  }


  /* Setup a fake graph to allow the rest of the code to work unchanged */
  nvtxs    = vtxdist[mype+1]-vtxdist[mype];
  firstvtx = vtxdist[mype];
  xadj     = imalloc(nvtxs+1, "ParMETIS_PartGeom: xadj");
  adjncy   = imalloc(nvtxs, "ParMETIS_PartGeom: adjncy");
  for (i=0; i<nvtxs; i++) {
    xadj[i] = i;
    adjncy[i] = firstvtx + (i+1)%nvtxs;
  }
  xadj[nvtxs] = nvtxs;

  graph = SetupGraph(ctrl, 1, vtxdist, xadj, NULL, NULL, adjncy, NULL, 0);


  /* Allocate workspace memory */
  AllocateWSpace(ctrl, 5*graph->nvtxs);


  /* Compute the initial geometric partitioning */
  STARTTIMER(ctrl, ctrl->TotalTmr);

  Coordinate_Partition(ctrl, graph, *ndims, xyz, 0);
  icopy(graph->nvtxs, graph->where, part);

  STOPTIMER(ctrl, ctrl->TotalTmr);
  IFSET(ctrl->dbglvl, DBG_TIME, PrintTimingInfo(ctrl));


  gk_free((void **)&xadj, (void **)&adjncy, LTERM);
  FreeInitialGraphAndRemap(graph);


DONE:
  FreeCtrl(&ctrl);
  if (gk_GetCurMemoryUsed() - curmem > 0) {
    printf("ParMETIS appears to have a memory leak of %zdbytes. Report this.\n",
        (ssize_t)(gk_GetCurMemoryUsed() - curmem));
  }
  gk_malloc_cleanup(0);

  return (int)status;
}




