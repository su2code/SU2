/*
 * pspases.c
 *
 * This file contains ordering routines that are to be used with the
 * parallel Cholesky factorization code PSPASES
 *
 * Started 10/14/97
 * George
 *
 * $Id: pspases.c 10535 2011-07-11 04:29:44Z karypis $
 *
 */

#include <parmetislib.h>


/***********************************************************************************
* This function is the entry point of the serial ordering algorithm.
************************************************************************************/
int ParMETIS_SerialNodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, 
        idx_t *numflag, idx_t *options, idx_t *order, idx_t *sizes, 
        MPI_Comm *comm)
{
  idx_t i, npes, mype;
  ctrl_t *ctrl=NULL;
  graph_t *agraph=NULL;
  idx_t *perm=NULL, *iperm=NULL;
  idx_t *sendcount, *displs;

  /* Setup the ctrl */
  ctrl = SetupCtrl(PARMETIS_OP_OMETIS, options, 1, 1, NULL, NULL, *comm);
  npes = ctrl->npes;
  mype = ctrl->mype;

  if (!ispow2(npes)) {
    if (mype == 0)
      printf("Error: The number of processors must be a power of 2!\n");
    FreeCtrl(&ctrl);
    return METIS_ERROR;
  }


  if (*numflag > 0) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 1);

  STARTTIMER(ctrl, ctrl->TotalTmr);
  STARTTIMER(ctrl, ctrl->MoveTmr);

  agraph = AssembleEntireGraph(ctrl, vtxdist, xadj, adjncy);

  STOPTIMER(ctrl, ctrl->MoveTmr);

  if (mype == 0) {
    perm  = imalloc(agraph->nvtxs, "PAROMETISS: perm");
    iperm = imalloc(agraph->nvtxs, "PAROMETISS: iperm");

    METIS_NodeNDP(agraph->nvtxs, agraph->xadj, agraph->adjncy, 
        agraph->vwgt, npes, NULL, perm, iperm, sizes);
  }

  STARTTIMER(ctrl, ctrl->MoveTmr);

  /* Broadcast the sizes array */
  gkMPI_Bcast((void *)sizes, 2*npes, IDX_T, 0, ctrl->gcomm);

  /* Scatter the iperm */
  sendcount = imalloc(npes, "PAROMETISS: sendcount");
  displs = imalloc(npes, "PAROMETISS: displs");
  for (i=0; i<npes; i++) {
    sendcount[i] = vtxdist[i+1]-vtxdist[i];
    displs[i] = vtxdist[i];
  }

  gkMPI_Scatterv((void *)iperm, sendcount, displs, IDX_T, (void *)order, 
      vtxdist[mype+1]-vtxdist[mype], IDX_T, 0, ctrl->gcomm);

  STOPTIMER(ctrl, ctrl->MoveTmr);
  STOPTIMER(ctrl, ctrl->TotalTmr);
  IFSET(ctrl->dbglvl, DBG_TIME, PrintTimingInfo(ctrl));
  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->gcomm));

  gk_free((void **)&agraph->xadj, &agraph->adjncy, &perm, &iperm, 
      &sendcount, &displs, &agraph, LTERM);

  if (*numflag > 0) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 0);

  goto DONE;

DONE:
  FreeCtrl(&ctrl);
  return METIS_OK;
}



/*************************************************************************
* This function assembles the graph into a single processor
**************************************************************************/
graph_t *AssembleEntireGraph(ctrl_t *ctrl, idx_t *vtxdist, idx_t *xadj, idx_t *adjncy)
{
  idx_t i, gnvtxs, nvtxs, gnedges, nedges;
  idx_t npes = ctrl->npes, mype = ctrl->mype;
  idx_t *axadj, *aadjncy;
  idx_t *recvcounts, *displs;
  graph_t *agraph;

  gnvtxs = vtxdist[npes];
  nvtxs = vtxdist[mype+1]-vtxdist[mype];
  nedges = xadj[nvtxs];

  recvcounts = imalloc(npes, "AssembleGraph: recvcounts");
  displs = imalloc(npes+1, "AssembleGraph: displs");

  /* Gather all the xadj arrays first */
  for (i=0; i<nvtxs; i++)
    xadj[i] = xadj[i+1]-xadj[i];

  axadj = imalloc(gnvtxs+1, "AssembleEntireGraph: axadj");

  for (i=0; i<npes; i++) {
    recvcounts[i] = vtxdist[i+1]-vtxdist[i];
    displs[i] = vtxdist[i];
  }

  /* Assemble the xadj and then the adjncy */
  gkMPI_Gatherv((void *)xadj, nvtxs, IDX_T, axadj, recvcounts, displs, 
      IDX_T, 0, ctrl->comm);

  MAKECSR(i, nvtxs, xadj);
  MAKECSR(i, gnvtxs, axadj);

  /* Gather all the adjncy arrays next */
  /* Determine the # of edges stored at each processor */
  gkMPI_Allgather((void *)(&nedges), 1, IDX_T, (void *)recvcounts, 1, IDX_T, ctrl->comm);
  
  displs[0] = 0;
  for (i=1; i<npes+1; i++) 
    displs[i] = displs[i-1] + recvcounts[i-1];
  gnedges = displs[npes];

  aadjncy = imalloc(gnedges, "AssembleEntireGraph: aadjncy");

  /* Assemble the xadj and then the adjncy */
  gkMPI_Gatherv((void *)adjncy, nedges, IDX_T, aadjncy, recvcounts, displs, IDX_T, 0, ctrl->comm);

  /* myprintf(ctrl, "Gnvtxs: %"PRIDX", Gnedges: %"PRIDX"\n", gnvtxs, gnedges); */

  agraph = CreateGraph();
  agraph->nvtxs = gnvtxs;
  agraph->nedges = gnedges;
  agraph->xadj = axadj;
  agraph->adjncy = aadjncy; 

  return agraph;
}
