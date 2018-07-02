/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * initmsection.c
 *
 * This file contains code that performs the k-way multisection
 *
 * Started 6/3/97
 * George
 *
 * $Id: initmsection.c 10361 2011-06-21 19:16:22Z karypis $
 */

#include <parmetislib.h>


#define DEBUG_IPART_


/************************************************************************************/
/*! 
 The entry point of the algorithm that finds the separators of the coarsest graph.
 This algorithm first assembles the graph to all the processors, it then splits the
 processors into groups depending on the number of partitions for which we want to
 compute the separator. The processors of each group compute the separator of their
 corresponding graph and the smallest separator is selected.

 The parameter nparts on calling this routine indicates the number of desired
 partitions after the multisection (excluding the nodes that end up on the
 separator). The initial bisection is achieved when nparts==2 and upon entry
 graph->where[] = 0 for all vertices. Similarly, if nparts==4, it indicates that we
 have a graph that is already partitioned into two parts (as specified in
 graph->where) and we need to find the separator of each one of these parts.

 The final partitioning is encoded in the graph->where vector as follows. If nparts
 is the number of partitions, the left, right, and separator subpartitions of the
 original partition i will be labeled 2*i, 2*i+1, and nparts+2*i, respectively. Note
 that in the above expressions, i goes from [0...nparts/2]. As a result of this
 encoding, the left (i.e., 0th) partition of a node \c i on the separator will be
 given by where[i]%nparts. 

*/
/************************************************************************************/
void InitMultisection(ctrl_t *ctrl, graph_t *graph) 
{ 
  idx_t i, myrank, mypart, options[METIS_NOPTIONS]; 
  idx_t *vtxdist, *gwhere = NULL, *part, *label; 
  graph_t *agraph; 
  idx_t *sendcounts, *displs; 
  MPI_Comm newcomm, labelcomm;
  struct {
    double cut;
    int rank;
  } lpecut, gpecut;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));
  WCOREPUSH;

  /* Assemble the graph and do the necessary pre-processing */
  agraph = AssembleMultisectedGraph(ctrl, graph);
  part   = agraph->where;
  agraph->where = NULL;

  /* Split the processors into groups so that each one can do a bisection */
  mypart = ctrl->mype%(ctrl->nparts/2);
  gkMPI_Comm_split(ctrl->comm, mypart, 0, &newcomm);
  gkMPI_Comm_rank(newcomm, &myrank);

  /* Each processor keeps the graph that it only needs and bisects it */
  KeepPart(ctrl, agraph, part, mypart);
  label = agraph->label;  /* Save this because ipart may need it */
  agraph->label = NULL;

  /* Bisect the graph and construct the separator */
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_SEED]    = (ctrl->mype+8)*101;
  options[METIS_OPTION_NSEPS]   = 5;
  options[METIS_OPTION_UFACTOR] = (idx_t)(1000.0*(ctrl->ubfrac - 1.0)); 

  WCOREPUSH;  /* for freeing agraph->where and gwhere */
  agraph->where = iwspacemalloc(ctrl, agraph->nvtxs);
  METIS_ComputeVertexSeparator(&agraph->nvtxs, agraph->xadj, agraph->adjncy, 
        agraph->vwgt, options, &agraph->mincut, agraph->where);

  for (i=0; i<agraph->nvtxs; i++) {
    PASSERT(ctrl, agraph->where[i]>=0 && agraph->where[i]<=2);
    if (agraph->where[i] == 2)
      agraph->where[i] = ctrl->nparts+2*mypart;
    else
      agraph->where[i] += 2*mypart;
  }

  /* Determine which PE got the minimum cut */
  lpecut.cut  = agraph->mincut;
  lpecut.rank = myrank;
  gkMPI_Allreduce(&lpecut, &gpecut, 1, MPI_DOUBLE_INT, MPI_MINLOC, newcomm);

  /* myprintf(ctrl, "Nvtxs: %"PRIDX", Mincut: %"PRIDX", GMincut: %"PRIDX", %"PRIDX"\n", agraph->nvtxs, 
         agraph->mincut, (idx_t)gpecut.cut, (idx_t)gpecut.rank); */

  /* Send the best where to the root processor of this partition */
  if (myrank != 0 && myrank == gpecut.rank) 
    gkMPI_Send((void *)agraph->where, agraph->nvtxs, IDX_T, 0, 1, newcomm);
  if (myrank == 0 && myrank != gpecut.rank)
    gkMPI_Recv((void *)agraph->where, agraph->nvtxs, IDX_T, gpecut.rank, 1, 
        newcomm, &ctrl->status);

  /* Create a communicator that stores all the i-th processors of the newcomm */
  gkMPI_Comm_split(ctrl->comm, myrank, 0, &labelcomm);

  /* Map the separator back to agraph. This is inefficient! */
  if (myrank == 0) {
    gwhere = iset(graph->gnvtxs, 0, iwspacemalloc(ctrl, graph->gnvtxs));
    for (i=0; i<agraph->nvtxs; i++)
      gwhere[label[i]] = agraph->where[i];

    gkMPI_Reduce((void *)gwhere, (void *)part, graph->gnvtxs, IDX_T, 
        MPI_SUM, 0, labelcomm);
  }

  WCOREPOP; /* free agraph->where & gwhere */

  agraph->where = part;

  /* The minimum PE performs the Scatter */
  vtxdist = graph->vtxdist;
  PASSERT(ctrl, graph->where != NULL);
  gk_free((void **)&graph->where, LTERM);  /* Remove the propagated down where info */
  graph->where = imalloc(graph->nvtxs+graph->nrecv, "InitPartition: where");

  sendcounts = iwspacemalloc(ctrl, ctrl->npes);
  displs     = iwspacemalloc(ctrl, ctrl->npes);

  for (i=0; i<ctrl->npes; i++) {
    sendcounts[i] = vtxdist[i+1]-vtxdist[i];
    displs[i]     = vtxdist[i];
  }

  gkMPI_Scatterv((void *)agraph->where, sendcounts, displs, IDX_T, 
               (void *)graph->where, graph->nvtxs, IDX_T, 0, ctrl->comm);

  agraph->label = label;
  FreeGraph(agraph);

  gkMPI_Comm_free(&newcomm);
  gkMPI_Comm_free(&labelcomm);

  WCOREPOP;

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));
}



/*************************************************************************/
/*! This function assembles the graph into a single processor */
/*************************************************************************/
graph_t *AssembleMultisectedGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, k, l, gnvtxs, nvtxs, gnedges, nedges, gsize;
  idx_t *xadj, *vwgt, *where, *adjncy, *adjwgt, *vtxdist, *imap;
  idx_t *axadj, *aadjncy, *aadjwgt, *avwgt, *awhere, *alabel;
  idx_t *mygraph, *ggraph;
  idx_t *recvcounts, *displs, mysize;
  graph_t *agraph;

  WCOREPUSH;

  gnvtxs  = graph->gnvtxs;
  nvtxs   = graph->nvtxs;
  nedges  = graph->xadj[nvtxs];
  xadj    = graph->xadj;
  vwgt    = graph->vwgt;
  where   = graph->where;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  vtxdist = graph->vtxdist;
  imap    = graph->imap;

  /* Determine the # of idx_t to receive from each processor */
  recvcounts = iwspacemalloc(ctrl, ctrl->npes);
  mysize = 3*nvtxs + 2*nedges;
  gkMPI_Allgather((void *)(&mysize), 1, IDX_T, (void *)recvcounts, 1, IDX_T, ctrl->comm);
  
  displs = iwspacemalloc(ctrl, ctrl->npes+1);
  for (displs[0]=0, i=1; i<ctrl->npes+1; i++) 
    displs[i] = displs[i-1] + recvcounts[i-1];

  /* allocate memory for the recv buffer of the assembled graph */
  gsize  = displs[ctrl->npes];
  ggraph = iwspacemalloc(ctrl, gsize);

  /* Construct the one-array storage format of the assembled graph */
  WCOREPUSH;  /* for freeing mygraph */
  mygraph = iwspacemalloc(ctrl, mysize);

  for (k=i=0; i<nvtxs; i++) {
    mygraph[k++] = xadj[i+1]-xadj[i];
    mygraph[k++] = vwgt[i];
    mygraph[k++] = where[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      mygraph[k++] = imap[adjncy[j]];
      mygraph[k++] = adjwgt[j];
    }
  }
  PASSERT(ctrl, mysize == k);

  /* Assemble the entire graph */
  gkMPI_Allgatherv((void *)mygraph, mysize, IDX_T, (void *)ggraph, 
      recvcounts, displs, IDX_T, ctrl->comm);

  WCOREPOP;  /* free mygraph */


  agraph = CreateGraph();
  agraph->nvtxs  = gnvtxs;
  agraph->ncon   = 1;
  agraph->nedges = gnedges = (gsize-3*gnvtxs)/2;

  /* Allocate memory for the assembled graph */
  axadj   = agraph->xadj   = imalloc(gnvtxs+1, "AssembleGraph: axadj");
  avwgt   = agraph->vwgt   = imalloc(gnvtxs, "AssembleGraph: avwgt");
  awhere  = agraph->where  = imalloc(gnvtxs, "AssembleGraph: awhere");
  aadjncy = agraph->adjncy = imalloc(gnedges, "AssembleGraph: adjncy");
  aadjwgt = agraph->adjwgt = imalloc(gnedges, "AssembleGraph: adjwgt");
  alabel  = agraph->label  = imalloc(gnvtxs, "AssembleGraph: alabel");

  for (k=j=i=0; i<gnvtxs; i++) {
    axadj[i] = ggraph[k++];
    avwgt[i] = ggraph[k++];
    awhere[i] = ggraph[k++];
    for (l=0; l<axadj[i]; l++) {
      aadjncy[j] = ggraph[k++];
      aadjwgt[j] = ggraph[k++];
      j++;
    }
  }

  /* Now fix up the received graph */
  MAKECSR(i, gnvtxs, axadj);

  iincset(gnvtxs, 0, alabel);

  WCOREPOP;

  return agraph;
}

