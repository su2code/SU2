/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * graph.c
 *
 * This file contains routines that deal with setting up and freeing the
 * graph structure.
 *
 * Started 2/24/96
 * George
 *
 * $Id: memory.c 10412 2011-06-25 23:12:57Z karypis $
 *
 */

#include <parmetislib.h>



/*************************************************************************/
/*! This function creates the graph from the user's inputs 
*/
/*************************************************************************/
graph_t *SetupGraph(ctrl_t *ctrl, idx_t ncon, idx_t *vtxdist, idx_t *xadj, 
             idx_t *vwgt, idx_t *vsize, idx_t *adjncy, idx_t *adjwgt, 
             idx_t wgtflag)
{
  idx_t i, j;
  graph_t *graph;

  graph          = CreateGraph();
  graph->level   = 0;
  graph->gnvtxs  = vtxdist[ctrl->npes];
  graph->nvtxs   = vtxdist[ctrl->mype+1]-vtxdist[ctrl->mype];
  graph->ncon    = ncon;
  graph->nedges  = xadj[graph->nvtxs];
  graph->xadj    = xadj;
  graph->vwgt    = vwgt;
  graph->vsize   = vsize;
  graph->adjncy  = adjncy;
  graph->adjwgt  = adjwgt;
  graph->vtxdist = vtxdist;


  /* allocate memory for weight arrays if not provided */
  if ((wgtflag&2) == 0 || vwgt == NULL) 
    graph->vwgt = ismalloc(graph->nvtxs*ncon, 1, "SetupGraph: vwgt");
  else 
    graph->free_vwgt = 0;

  if ((wgtflag&1) == 0 || adjwgt == NULL) 
    graph->adjwgt = ismalloc(graph->nedges, 1, "SetupGraph: adjwgt");
  else 
    graph->free_adjwgt = 0;


  /* allocate memory for special arrays that apply only to some optypes */
  if (ctrl->optype == PARMETIS_OP_AMETIS || ctrl->optype == PARMETIS_OP_RMETIS) { 
    if (vsize == NULL) 
      graph->vsize = ismalloc(graph->nvtxs, 1, "vsize");
    else
      graph->free_vsize = 0;

    graph->home = ismalloc(graph->nvtxs, 1, "home");

    /* determine edge_size_ratio */
    ctrl->edge_size_ratio = 
        (.1+(real_t)GlobalSESum(ctrl, isum(graph->nedges, graph->adjwgt, 1))) /
        (.1+(real_t)GlobalSESum(ctrl, isum(graph->nvtxs,  graph->vsize,  1)));
  }


  /* compute invtvwgts */
  SetupCtrl_invtvwgts(ctrl, graph); 


  /* compute nvwgts */
  SetupGraph_nvwgts(ctrl, graph); 

  return graph;
}


/*************************************************************************/
/*! This function creates the nvwgt of the graph 
*/
/*************************************************************************/
void SetupGraph_nvwgts(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, nvtxs, ncon;
  idx_t *vwgt;
  real_t *nvwgt, *invtvwgts;

  nvtxs = graph->nvtxs;
  ncon  = graph->ncon;
  vwgt  = graph->vwgt;

  invtvwgts = ctrl->invtvwgts;

  /* compute nvwgts */
  nvwgt = graph->nvwgt = rmalloc(nvtxs*ncon, "SetupGraph_nvwgts: graph->nvwgt");
  for (i=0; i<nvtxs; i++) {
    for (j=0; j<ncon; j++)
      nvwgt[i*ncon+j] = invtvwgts[j]*vwgt[i*ncon+j];
  }
}



/*************************************************************************/
/*! This function creates a Coarsegraph_t data structure and initializes
    the various fields.
*/
/*************************************************************************/
graph_t *CreateGraph(void)
{
  graph_t *graph;

  graph = (graph_t *)gk_malloc(sizeof(graph_t), "CreateCoarseGraph: graph");

  InitGraph(graph);

  return graph;
}


/*************************************************************************/
/*! This function creates a Coarsegraph_t data structure and initializes
    the various fields
*/
/*************************************************************************/
void InitGraph(graph_t *graph) 
{
  memset(graph, 0, sizeof(graph_t));

  graph->gnvtxs = graph->nvtxs = graph->nedges = graph->nsep = -1;
  graph->nnbrs = graph->nrecv = graph->nsend = graph->nlocal = -1;
  graph->xadj = graph->vwgt = graph->vsize = graph->adjncy = graph->adjwgt = NULL;
  graph->nvwgt = NULL;
  graph->vtxdist = NULL;
  graph->match = graph->cmap = NULL;
  graph->label = NULL;

  graph->peind = NULL;
  graph->sendptr = graph->sendind = graph->recvptr = graph->recvind = NULL;
  graph->imap = NULL;
  graph->pexadj = graph->peadjncy = graph->peadjloc = NULL;
  graph->lperm = NULL;

  graph->slens = graph->rlens = NULL;
  graph->rcand = NULL;

  graph->where = graph->home = graph->lpwgts = graph->gpwgts = NULL;
  graph->lnpwgts = graph->gnpwgts = NULL;
  graph->ckrinfo = NULL;

  graph->nrinfo  = NULL;
  graph->sepind  = NULL;

  graph->coarser = graph->finer = NULL;

  graph->free_vwgt = graph->free_adjwgt = graph->free_vsize = 1;
}


/*************************************************************************/
/*! This function deallocates any memory stored in a graph 
*/
/*************************************************************************/
void FreeGraph(graph_t *graph) 
{

  /* Graph structure fields */
  gk_free((void **)&graph->xadj, 
         (void **)&graph->vwgt,
         (void **)&graph->nvwgt,
         (void **)&graph->vsize,
         (void **)&graph->adjncy,
         (void **)&graph->adjwgt,
         (void **)&graph->vtxdist, 
         (void **)&graph->home, 
         LTERM);

  FreeNonGraphFields(graph); 

  gk_free((void **)&graph, LTERM);
}


/*************************************************************************/
/*! This function deallocates the non-graph structure fields of a graph
    data structure */
/*************************************************************************/
void FreeNonGraphFields(graph_t *graph) 
{

  gk_free(
      /* Coarsening fields */
      (void **)&graph->match, 
      (void **)&graph->cmap, 

      /* Initial partitioning fields */
      (void **)&graph->label, 

      /* Communication/Setup fields */
      (void **)&graph->peind, 
      (void **)&graph->sendptr, 
      (void **)&graph->sendind, 
      (void **)&graph->recvptr, 
      (void **)&graph->recvind, 
      (void **)&graph->imap,
      (void **)&graph->pexadj,
      (void **)&graph->peadjncy,
      (void **)&graph->peadjloc,
      (void **)&graph->lperm, 

      /* Projection fields */
      (void **)&graph->rlens,
      (void **)&graph->slens,
      (void **)&graph->rcand,

      /* Refinement fields */
      (void **)&graph->where, 
      (void **)&graph->lpwgts, 
      (void **)&graph->gpwgts, 
      (void **)&graph->lnpwgts, 
      (void **)&graph->gnpwgts, 
      (void **)&graph->ckrinfo, 
      (void **)&graph->nrinfo, 
      (void **)&graph->sepind,

      LTERM);
}


/*************************************************************************/
/*! This function deallocates the non-graph and non-setup structure fields 
    of a graph data structure */
/*************************************************************************/
void FreeNonGraphNonSetupFields(graph_t *graph) 
{

  gk_free(
      /* Coarsening fields */
      (void **)&graph->match, 
      (void **)&graph->cmap, 

      /* Initial partitioning fields */
      (void **)&graph->label, 

      /* Projection fields */
      (void **)&graph->rlens,
      (void **)&graph->slens,
      (void **)&graph->rcand,

      /* Refinement fields */
      (void **)&graph->where, 
      (void **)&graph->lpwgts, 
      (void **)&graph->gpwgts, 
      (void **)&graph->lnpwgts, 
      (void **)&graph->gnpwgts, 
      (void **)&graph->ckrinfo, 
      (void **)&graph->nrinfo, 
      (void **)&graph->sepind,

      LTERM);
}


/*************************************************************************/
/*! This function frees any memory allocated for storing the initial graph
    and performs the local to global (i.e., original numbering of the
    adjacency list)
*/
/*************************************************************************/
void FreeInitialGraphAndRemap(graph_t *graph) 
{
  idx_t i, nedges;
  idx_t *adjncy, *imap;

  nedges = graph->nedges;
  adjncy = graph->adjncy;
  imap   = graph->imap;

  if (imap != NULL) {
    for (i=0; i<nedges; i++)
      adjncy[i] = imap[adjncy[i]];  /* Apply local to global transformation */
  }

  /* Free fields that are not related to the structure of the graph */
  FreeNonGraphFields(graph); 

  /* Free some derived graph-structure fields */
  gk_free((void **)&graph->nvwgt, &graph->home, &graph->lnpwgts, 
      &graph->gnpwgts, LTERM);

  if (graph->free_vwgt)
    gk_free((void **)&graph->vwgt, LTERM);
  if (graph->free_adjwgt)
    gk_free((void **)&graph->adjwgt, LTERM);
  if (graph->free_vsize)
    gk_free((void **)&graph->vsize, LTERM);

  gk_free((void **)&graph, LTERM);
}
