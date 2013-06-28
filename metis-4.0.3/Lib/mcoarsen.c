/*
 * mcoarsen.c
 *
 * This file contains the driving routines for the coarsening process 
 *
 * Started 7/23/97
 * George
 *
 * $Id: mcoarsen.c,v 1.2 2003/07/31 16:23:29 karypis Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function takes a graph and creates a sequence of coarser graphs
**************************************************************************/
GraphType *MCCoarsen2Way(CtrlType *ctrl, GraphType *graph)
{
  int i, clevel;
  GraphType *cgraph;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->CoarsenTmr));

  cgraph = graph;

  clevel = 0;
  do {
    if (ctrl->dbglvl&DBG_COARSEN) {
      printf("%6d %7d %10d [%d] [%6.4f", cgraph->nvtxs, cgraph->nedges, 
              idxsum(cgraph->nvtxs, cgraph->adjwgtsum), ctrl->CoarsenTo, ctrl->nmaxvwgt);
      for (i=0; i<graph->ncon; i++)
        printf(" %5.3f", ssum_strd(cgraph->nvtxs, cgraph->nvwgt+i, cgraph->ncon));
      printf("]\n");
    }

    switch (ctrl->CType) {
      case MATCH_RM:
        MCMatch_RM(ctrl, cgraph);
        break;
      case MATCH_HEM:
        if (clevel < 1 || cgraph->nedges == 0)
          MCMatch_RM(ctrl, cgraph);
        else
          MCMatch_HEM(ctrl, cgraph);
        break;
      case MATCH_SHEM:
        if (clevel < 1 || cgraph->nedges == 0)
          MCMatch_RM(ctrl, cgraph);
        else
          MCMatch_SHEM(ctrl, cgraph);
        break;
      case MATCH_SHEMKWAY:
        if (clevel < 1 || cgraph->nedges == 0)
          MCMatch_RM(ctrl, cgraph);
        else
          MCMatch_SHEM(ctrl, cgraph);
        break;
      case MATCH_SHEBM_ONENORM:
        if (clevel < 1 || cgraph->nedges == 0)
          MCMatch_RM(ctrl, cgraph);
        else
          MCMatch_SHEBM(ctrl, cgraph, 1);
        break;
      case MATCH_SHEBM_INFNORM:
        if (clevel < 1 || cgraph->nedges == 0)
          MCMatch_RM(ctrl, cgraph);
        else
          MCMatch_SHEBM(ctrl, cgraph, -1);
        break;
      case MATCH_SBHEM_ONENORM:
        if (clevel < 1 || cgraph->nedges == 0)
          MCMatch_RM(ctrl, cgraph);
        else
          MCMatch_SBHEM(ctrl, cgraph, 1);
        break;
      case MATCH_SBHEM_INFNORM:
        if (clevel < 1 || cgraph->nedges == 0)
          MCMatch_RM(ctrl, cgraph);
        else
          MCMatch_SBHEM(ctrl, cgraph, -1);
        break;
      default:
        errexit("Unknown CType: %d\n", ctrl->CType);
    }

    cgraph = cgraph->coarser;
    clevel++;

  } while (cgraph->nvtxs > ctrl->CoarsenTo && cgraph->nvtxs < COARSEN_FRACTION2*cgraph->finer->nvtxs && cgraph->nedges > cgraph->nvtxs/2); 

  if (ctrl->dbglvl&DBG_COARSEN) {
    printf("%6d %7d %10d [%d] [%6.4f", cgraph->nvtxs, cgraph->nedges, 
            idxsum(cgraph->nvtxs, cgraph->adjwgtsum), ctrl->CoarsenTo, ctrl->nmaxvwgt);
    for (i=0; i<graph->ncon; i++)
      printf(" %5.3f", ssum_strd(cgraph->nvtxs, cgraph->nvwgt+i, cgraph->ncon));
    printf("]\n");
  }


  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->CoarsenTmr));

  return cgraph;
}

