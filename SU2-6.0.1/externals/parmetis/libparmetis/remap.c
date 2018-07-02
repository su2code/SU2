/*
 * premap.c
 *
 * This file contains code that computes the assignment of processors to
 * partition numbers so that it will minimize the redistribution cost
 *
 * Started 4/16/98
 * George
 *
 * $Id: remap.c 10361 2011-06-21 19:16:22Z karypis $
 *
 */

#include <parmetislib.h>

/*************************************************************************
* This function remaps that graph so that it will minimize the 
* redistribution cost
**************************************************************************/
void ParallelReMapGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, nvtxs, nparts;
  idx_t *where, *vsize, *map, *lpwgts;

  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->RemapTmr));

  if (ctrl->npes != ctrl->nparts) {
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->RemapTmr));
    return;
  }

  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  where  = graph->where;
  vsize  = graph->vsize;
  nparts = ctrl->nparts;

  map    = iwspacemalloc(ctrl, nparts);
  lpwgts = iset(nparts, 0, iwspacemalloc(ctrl, nparts));

  for (i=0; i<nvtxs; i++)
    lpwgts[where[i]] += (vsize == NULL) ? 1 : vsize[i];

  ParallelTotalVReMap(ctrl, lpwgts, map, NREMAP_PASSES, graph->ncon);

  for (i=0; i<nvtxs; i++)
    where[i] = map[where[i]];

  WCOREPOP;

  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->RemapTmr));
}


/*************************************************************************
* This function computes the assignment using the the objective the 
* minimization of the total volume of data that needs to move
**************************************************************************/
void ParallelTotalVReMap(ctrl_t *ctrl, idx_t *lpwgts, idx_t *map, idx_t npasses, idx_t ncon)
{
  idx_t i, ii, j, k, nparts, mype;
  idx_t pass, maxipwgt, nmapped, oldwgt, newwgt, done;
  idx_t *rowmap, *mylpwgts;
  ikv_t *recv, send;
  idx_t nsaved, gnsaved;

  WCOREPUSH;

  mype   = ctrl->mype;
  nparts = ctrl->nparts;

  rowmap   = iset(nparts, -1, iwspacemalloc(ctrl, nparts));
  mylpwgts = icopy(nparts, lpwgts, iwspacemalloc(ctrl, nparts));
  recv     = ikvwspacemalloc(ctrl, nparts);

  iset(nparts, -1, map);

  done = nmapped = 0;
  for (pass=0; pass<npasses; pass++) {
    maxipwgt = iargmax(nparts, mylpwgts);

    if (mylpwgts[maxipwgt] > 0 && !done) {
      send.key = -mylpwgts[maxipwgt];
      send.val = mype*nparts+maxipwgt;
    }
    else {
      send.key = 0;
      send.val = -1;
    }

    /* each processor sends its selection */
    gkMPI_Allgather((void *)&send, 2, IDX_T, (void *)recv, 2, IDX_T, ctrl->comm); 

    ikvsorti(nparts, recv);
    if (recv[0].key == 0)
      break;

    /* now make as many assignments as possible */
    for (ii=0; ii<nparts; ii++) {
      i = recv[ii].val;

      if (i == -1)
        continue;

      j = i%nparts;
      k = i/nparts;
      if (map[j] == -1 && rowmap[k] == -1 && SimilarTpwgts(ctrl->tpwgts, ncon, j, k)) {
        map[j] = k;
        rowmap[k] = j;
        nmapped++;
        mylpwgts[j] = 0;
        if (mype == k)
          done = 1;
      }

      if (nmapped == nparts)
        break;
    }

    if (nmapped == nparts)
      break;
  }

  /* Map unmapped partitions */
  if (nmapped < nparts) {
    for (i=j=0; j<nparts && nmapped<nparts; j++) {
      if (map[j] == -1) {
        for (; i<nparts; i++) {
          if (rowmap[i] == -1 && SimilarTpwgts(ctrl->tpwgts, ncon, i, j)) {
            map[j] = i;
            rowmap[i] = j;
            nmapped++;
            break;
          }
        }
      }
    }
  }

  /* check to see if remapping fails (due to dis-similar tpwgts) */
  /* if remapping fails, revert to original mapping */
  if (nmapped < nparts) {
    for (i=0; i<nparts; i++)
      map[i] = i; 
    IFSET(ctrl->dbglvl, DBG_REMAP, rprintf(ctrl, "Savings from parallel remapping: %0\n")); 
  }
  else {
    /* check for a savings */
    oldwgt  = lpwgts[mype];
    newwgt  = lpwgts[rowmap[mype]];
    nsaved  = newwgt - oldwgt;
    gnsaved = GlobalSESum(ctrl, nsaved);

    /* undo everything if we don't see a savings */
    if (gnsaved <= 0) {
      for (i=0; i<nparts; i++)
        map[i] = i;
    }
    IFSET(ctrl->dbglvl, DBG_REMAP, rprintf(ctrl, 
          "Savings from parallel remapping: %"PRIDX"\n",gk_max(0,gnsaved))); 
  }

  WCOREPOP;
}


/*************************************************************************
* This function computes the assignment using the the objective the
* minimization of the total volume of data that needs to move
**************************************************************************/
idx_t SimilarTpwgts(real_t *tpwgts, idx_t ncon, idx_t s1, idx_t s2)
{
  idx_t i;

  for (i=0; i<ncon; i++)
    if (fabs(tpwgts[s1*ncon+i]-tpwgts[s2*ncon+i]) > SMALLFLOAT)
      break;

  if (i == ncon)
    return 1;

  return 0;
}

