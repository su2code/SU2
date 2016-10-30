/*
 * serial.c
 *
 * This file contains code that implements k-way refinement
 *
 * Started 7/28/97
 * George
 *
 * $Id: serial.c 13927 2013-03-27 22:42:41Z karypis $
 *
 */

#include <parmetislib.h>


/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void Mc_ComputeSerialPartitionParams(ctrl_t *ctrl, graph_t *graph, idx_t nparts)
{
  idx_t i, j, k;
  idx_t nvtxs, nedges, ncon, mincut, me, other;
  idx_t *xadj, *adjncy, *adjwgt, *where;
  ckrinfo_t *myrinfo;
  cnbr_t *mynbrs;
  real_t *nvwgt, *npwgts;
  idx_t mype;

  gkMPI_Comm_rank(MPI_COMM_WORLD, &mype);

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  nvwgt  = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  npwgts = rset(ncon*nparts, 0.0, graph->gnpwgts);

  PASSERT(ctrl, graph->ckrinfo != NULL);
  PASSERT(ctrl, ctrl->cnbrpool != NULL);

  memset(graph->ckrinfo, 0, sizeof(ckrinfo_t)*nvtxs);
  cnbrpoolReset(ctrl);

  /*------------------------------------------------------------
  / Compute now the id/ed degrees
  /------------------------------------------------------------*/
  nedges = mincut = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    raxpy(ncon, 1.0, nvwgt+i*ncon, 1, npwgts+me*ncon, 1);

    myrinfo = graph->ckrinfo+i;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[adjncy[j]]) {
        myrinfo->id += adjwgt[j];
      }
      else {
        myrinfo->ed += adjwgt[j];
      }
    }

    mincut += myrinfo->ed;

    /* Time to compute the particular external degrees */
    if (myrinfo->ed > 0) {
      myrinfo->inbr = cnbrpoolGetNext(ctrl, xadj[i+1]-xadj[i]+1);
      mynbrs        = ctrl->cnbrpool + myrinfo->inbr;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (me != other) {
          for (k=0; k<myrinfo->nnbrs; k++) {
            if (mynbrs[k].pid == other) {
              mynbrs[k].ed += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->nnbrs) {
            mynbrs[k].pid = other;
            mynbrs[k].ed  = adjwgt[j];
            myrinfo->nnbrs++;
          }
        }
      }
    }
    else {
      myrinfo->inbr = -1;
    }
  }

  graph->mincut = mincut/2;

  return;
}


/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void Mc_SerialKWayAdaptRefine(ctrl_t *ctrl, graph_t *graph, idx_t nparts, 
          idx_t *home, real_t *orgubvec, idx_t npasses)
{
  idx_t i, ii, iii, j, k;
  idx_t nvtxs, ncon, pass, nmoves;
  idx_t from, me, myhome, to, oldcut, gain, tmp;
  idx_t *xadj, *adjncy, *adjwgt;
  idx_t *where;
  real_t *npwgts, *nvwgt, *minwgt, *maxwgt, *ubvec;
  idx_t gain_is_greater, gain_is_same, fit_in_to, fit_in_from, going_home;
  idx_t zero_gain, better_balance_ft, better_balance_tt;
  ikv_t *cand;
  idx_t mype;
  ckrinfo_t *myrinfo;
  cnbr_t *mynbrs;

  WCOREPUSH;

  gkMPI_Comm_rank(MPI_COMM_WORLD, &mype);

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;
  npwgts = graph->gnpwgts;
  
  /* Setup the weight intervals of the various subdomains */
  cand   = ikvwspacemalloc(ctrl, nvtxs);
  minwgt = rwspacemalloc(ctrl, nparts*ncon);
  maxwgt = rwspacemalloc(ctrl, nparts*ncon);
  ubvec  = rwspacemalloc(ctrl, ncon);

  ComputeHKWayLoadImbalance(ncon, nparts, npwgts, ubvec);
  for (i=0; i<ncon; i++)
    ubvec[i] =gk_max(ubvec[i], orgubvec[i]);

  for (i=0; i<nparts; i++) {
    for (j=0; j<ncon; j++) {
      maxwgt[i*ncon+j] = ubvec[j]/(real_t)nparts;
      minwgt[i*ncon+j] = ubvec[j]*(real_t)nparts;
    }
  }

  for (pass=0; pass<npasses; pass++) {
    oldcut = graph->mincut;

    for (i=0; i<nvtxs; i++) {
      cand[i].key = graph->ckrinfo[i].ed-graph->ckrinfo[i].id;
      cand[i].val = i;
    }
    ikvsortd(nvtxs, cand);

    nmoves = 0;
    for (iii=0; iii<nvtxs; iii++) {
      i = cand[iii].val;

      myrinfo = graph->ckrinfo+i;

      if (myrinfo->ed >= myrinfo->id) {
        from   = where[i];
        myhome = home[i];
        nvwgt  = graph->nvwgt+i*ncon;

        if (myrinfo->id > 0 &&
            AreAllHVwgtsBelow(ncon, 1.0, npwgts+from*ncon, -1.0, nvwgt, minwgt+from*ncon)) 
          continue;

        mynbrs = ctrl->cnbrpool + myrinfo->inbr;

        for (k=myrinfo->nnbrs-1; k>=0; k--) {
          to = mynbrs[k].pid;
          gain = mynbrs[k].ed - myrinfo->id; 
          if (gain >= 0 && 
             (AreAllHVwgtsBelow(ncon, 1.0, npwgts+to*ncon, 1.0, nvwgt, maxwgt+to*ncon) ||
             IsHBalanceBetterFT(ncon,npwgts+from*ncon,npwgts+to*ncon,nvwgt,ubvec))) {
            break;
          }
        }

        /* break out if you did not find a candidate */
        if (k < 0)
          continue;

        for (j=k-1; j>=0; j--) {
          to = mynbrs[j].pid;
          going_home        = (myhome == to);
          gain_is_same      = (mynbrs[j].ed == mynbrs[k].ed);
          gain_is_greater   = (mynbrs[j].ed > mynbrs[k].ed);
          fit_in_to         = AreAllHVwgtsBelow(ncon, 1.0, npwgts+to*ncon, 1.0, nvwgt, 
                                  maxwgt+to*ncon);
          better_balance_ft = IsHBalanceBetterFT(ncon, npwgts+from*ncon,
                                  npwgts+to*ncon, nvwgt, ubvec);
          better_balance_tt = IsHBalanceBetterTT(ncon, npwgts+mynbrs[k].pid*ncon,
                                  npwgts+to*ncon,nvwgt,ubvec);

          if (
               (gain_is_greater &&
                 (fit_in_to ||
                  better_balance_ft)
               )
            ||
               (gain_is_same &&
                 (
                   (fit_in_to &&
                    going_home)
                ||
                    better_balance_tt
                 )
               )
             ) {
            k = j;
          }
        }

        to = mynbrs[k].pid;
        going_home = (myhome == to);
        zero_gain  = (mynbrs[k].ed == myrinfo->id);

        fit_in_from       = AreAllHVwgtsBelow(ncon, 1.0, npwgts+from*ncon, 0.0, 
                                npwgts+from*ncon, maxwgt+from*ncon);
        better_balance_ft = IsHBalanceBetterFT(ncon, npwgts+from*ncon,
                                npwgts+to*ncon, nvwgt, ubvec);

        if (zero_gain &&
            !going_home &&
            !better_balance_ft &&
            fit_in_from)
          continue;

        /*=====================================================================
        * If we got here, we can now move the vertex from 'from' to 'to' 
        *======================================================================*/
        graph->mincut -= mynbrs[k].ed-myrinfo->id;

        /* Update where, weight, and ID/ED information of the vertex you moved */
        raxpy(ncon,  1.0, nvwgt, 1, npwgts+to*ncon,   1);
        raxpy(ncon, -1.0, nvwgt, 1, npwgts+from*ncon, 1);
        where[i] = to;
        myrinfo->ed += myrinfo->id-mynbrs[k].ed;
        gk_SWAP(myrinfo->id, mynbrs[k].ed, tmp);

        if (mynbrs[k].ed == 0) 
          mynbrs[k] = mynbrs[--myrinfo->nnbrs];
        else
          mynbrs[k].pid = from;

        /* Update the degrees of adjacent vertices */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          ii = adjncy[j];
          me = where[ii];

          myrinfo = graph->ckrinfo+ii;
          if (myrinfo->inbr == -1) {
            myrinfo->inbr  = cnbrpoolGetNext(ctrl, xadj[ii+1]-xadj[ii]+1);
            myrinfo->nnbrs = 0;
          }
          mynbrs = ctrl->cnbrpool + myrinfo->inbr;

          if (me == from) {
            INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);
          }
          else {
            if (me == to) {
              INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);
            }
          }

          /* Remove contribution of the ed from 'from' */
          if (me != from) {
            for (k=0; k<myrinfo->nnbrs; k++) {
              if (mynbrs[k].pid == from) {
                if (mynbrs[k].ed == adjwgt[j]) 
                  mynbrs[k] = mynbrs[--myrinfo->nnbrs];
                else
                  mynbrs[k].ed -= adjwgt[j];
                break;
              }
            }
          }

          /* Add contribution of the ed to 'to' */
          if (me != to) {
            for (k=0; k<myrinfo->nnbrs; k++) {
              if (mynbrs[k].pid == to) {
                mynbrs[k].ed += adjwgt[j];
                break;
              }
            }
            if (k == myrinfo->nnbrs) {
              mynbrs[k].pid = to;
              mynbrs[k].ed  = adjwgt[j];
              myrinfo->nnbrs++;
            }
          }
        }
        nmoves++;
      }
    }

    if (graph->mincut == oldcut)
      break;
  }

  WCOREPOP;

  return;
}



/*************************************************************************
* This function checks if the vertex weights of two vertices are below 
* a given set of values
**************************************************************************/
idx_t AreAllHVwgtsBelow(idx_t ncon, real_t alpha, real_t *vwgt1, real_t beta, 
          real_t *vwgt2, real_t *limit)
{
  idx_t i;

  for (i=0; i<ncon; i++)
    if (alpha*vwgt1[i] + beta*vwgt2[i] > limit[i])
      return 0;

  return 1;
}


/*************************************************************************
* This function computes the load imbalance over all the constrains
* For now assume that we just want balanced partitionings
**************************************************************************/ 
void ComputeHKWayLoadImbalance(idx_t ncon, idx_t nparts, real_t *npwgts, real_t *lbvec)
{
  idx_t i, j;
  real_t max;

  for (i=0; i<ncon; i++) {
    max = 0.0;
    for (j=0; j<nparts; j++) {
      if (npwgts[j*ncon+i] > max)
        max = npwgts[j*ncon+i];
    }

    lbvec[i] = max*nparts;
  }
}


/**************************************************************
*  This subroutine remaps a partitioning on a single processor
**************************************************************/
void SerialRemap(ctrl_t *ctrl, graph_t *graph, idx_t nparts, 
         idx_t *base, idx_t *scratch, idx_t *remap, real_t *tpwgts)
{
  idx_t i, ii, j, k;
  idx_t nvtxs, nmapped, max_mult;
  idx_t from, to, current_from, smallcount, bigcount;
  ikv_t *flowto, *bestflow;
  i2kv_t *sortvtx;
  idx_t *vsize, *htable, *map, *rowmap;

  WCOREPUSH;

  nvtxs = graph->nvtxs;
  vsize = graph->vsize;
  max_mult = gk_min(MAX_NPARTS_MULTIPLIER, nparts);

  sortvtx      = (i2kv_t *)wspacemalloc(ctrl, nvtxs*sizeof(i2kv_t));
  flowto       = ikvwspacemalloc(ctrl, nparts);
  bestflow     = ikvwspacemalloc(ctrl, nparts*max_mult);
  map = htable = iset(2*nparts, -1, iwspacemalloc(ctrl, 2*nparts));
  rowmap = map+nparts;

  for (i=0; i<nvtxs; i++) {
    sortvtx[i].key1 = base[i];
    sortvtx[i].key2 = vsize[i];
    sortvtx[i].val = i;
  }

  qsort((void *)sortvtx, (size_t)nvtxs, (size_t)sizeof(i2kv_t), SSMIncKeyCmp);

  for (j=0; j<nparts; j++) {
    flowto[j].key = 0;
    flowto[j].val = j;
  }

  /* this step has nparts*nparts*log(nparts) computational complexity */
  bigcount = smallcount = current_from = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = sortvtx[ii].val;
    from = base[i];
    to = scratch[i];

    if (from > current_from) {
      /* reset the hash table */
      for (j=0; j<smallcount; j++)
        htable[flowto[j].val] = -1;
      ASSERT(isum(nparts, htable, 1) == -nparts);

      ikvsorti(smallcount, flowto);

      for (j=0; j<gk_min(smallcount, max_mult); j++, bigcount++) {
        bestflow[bigcount].key = flowto[j].key;
        bestflow[bigcount].val = current_from*nparts+flowto[j].val;
      }

      smallcount = 0;
      current_from = from;
    }

    if (htable[to] == -1) {
      htable[to] = smallcount;
      flowto[smallcount].key = -vsize[i];
      flowto[smallcount].val = to;
      smallcount++;
    }
    else {
      flowto[htable[to]].key += -vsize[i];
    }
  }

  /* reset the hash table */
  for (j=0; j<smallcount; j++)
    htable[flowto[j].val] = -1;
  ASSERT(isum(nparts, htable, 1) == -nparts);

  ikvsorti(smallcount, flowto);

  for (j=0; j<gk_min(smallcount, max_mult); j++, bigcount++) {
    bestflow[bigcount].key = flowto[j].key;
    bestflow[bigcount].val = current_from*nparts+flowto[j].val;
  }
  ikvsorti(bigcount, bestflow);

  ASSERT(isum(nparts, map, 1) == -nparts);
  ASSERT(isum(nparts, rowmap, 1) == -nparts);
  nmapped = 0;

  /* now make as many assignments as possible */
  for (ii=0; ii<bigcount; ii++) {
    i = bestflow[ii].val;
    j = i % nparts;  /* to */
    k = i / nparts;  /* from */

    if (map[j] == -1 && rowmap[k] == -1 && SimilarTpwgts(tpwgts, graph->ncon, j, k)) {
      map[j] = k;
      rowmap[k] = j;
      nmapped++;
    }

    if (nmapped == nparts)
      break;
  }


  /* remap the rest */
  /* it may help try remapping to the same label first */
  if (nmapped < nparts) {
    for (j=0; j<nparts && nmapped<nparts; j++) {
      if (map[j] == -1) {
        for (ii=0; ii<nparts; ii++) {
          i = (j+ii) % nparts;
          if (rowmap[i] == -1 && SimilarTpwgts(tpwgts, graph->ncon, i, j)) {
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
  if (nmapped < nparts)
    for (i=0; i<nparts; i++)
      map[i] = i;

  for (i=0; i<nvtxs; i++)
    remap[i] = map[remap[i]];

  WCOREPOP;
}


/*************************************************************************
*  This is a comparison function for Serial Remap
**************************************************************************/
int SSMIncKeyCmp(const void *fptr, const void *sptr)
{
  i2kv_t *first, *second;

  first = (i2kv_t *)(fptr);
  second = (i2kv_t *)(sptr);

  if (first->key1 > second->key1)
    return 1;

  if (first->key1 < second->key1)
     return -1;

  if (first->key2 < second->key2)
    return 1;

  if (first->key2 > second->key2)
     return -1;

   return 0;
}


/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
void Mc_Serial_FM_2WayRefine(ctrl_t *ctrl, graph_t *graph, real_t *tpwgts, idx_t npasses)
{
  idx_t i, ii, j, k;
  idx_t kwgt, nvtxs, ncon, nbnd, nswaps, from, to, pass, limit, tmp, cnum;
  idx_t *xadj, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind;
  idx_t *moved, *swaps, *qnum;
  real_t *nvwgt, *npwgts, *mindiff, *tmpdiff, origbal, minbal, newbal;
  rpq_t **parts[2];
  idx_t higain, mincut, initcut, newcut, mincutorder;
  real_t *rtpwgts;
  idx_t mype;

  WCOREPUSH;

  gkMPI_Comm_rank(MPI_COMM_WORLD, &mype);

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  nvwgt  = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;
  id     = graph->sendind;
  ed     = graph->recvind;
  npwgts = graph->gnpwgts;
  bndptr = graph->sendptr;
  bndind = graph->recvptr;

  mindiff  = rwspacemalloc(ctrl, ncon);
  tmpdiff  = rwspacemalloc(ctrl, ncon);
  rtpwgts  = rwspacemalloc(ctrl, 2*ncon);
  parts[0] = (rpq_t **)wspacemalloc(ctrl, sizeof(rpq_t *)*ncon);
  parts[1] = (rpq_t **)wspacemalloc(ctrl, sizeof(rpq_t *)*ncon);

  moved   = iwspacemalloc(ctrl, nvtxs);
  swaps   = iwspacemalloc(ctrl, nvtxs);
  qnum    = iwspacemalloc(ctrl, nvtxs);

  limit = gk_min(gk_max(0.01*nvtxs, 25), 150);

  /* Initialize the queues */
  for (i=0; i<ncon; i++) {
    parts[0][i] = rpqCreate(nvtxs);
    parts[1][i] = rpqCreate(nvtxs);
  }
  for (i=0; i<nvtxs; i++)
    qnum[i] = rargmax(ncon, nvwgt+i*ncon);

  origbal = Serial_Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);

  for (i=0; i<ncon; i++) {
    rtpwgts[i] = origbal*tpwgts[i];
    rtpwgts[ncon+i] = origbal*tpwgts[ncon+i];
  }

  iset(nvtxs, -1, moved);
  for (pass=0; pass<npasses; pass++) { /* Do a number of passes */
    for (i=0; i<ncon; i++) {
      rpqReset(parts[0][i]);
      rpqReset(parts[1][i]);
    }

    mincutorder = -1;
    newcut = mincut = initcut = graph->mincut;
    for (i=0; i<ncon; i++)
      mindiff[i] = fabs(tpwgts[i]-npwgts[i]);
    minbal = Serial_Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);

    /* Insert boundary nodes in the priority queues */
    nbnd = graph->gnvtxs;
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[ii];
      rpqInsert(parts[where[i]][qnum[i]], i, (real_t)(ed[i]-id[i]));
    }

    for (nswaps=0; nswaps<nvtxs; nswaps++) {
      Serial_SelectQueue(ncon, npwgts, rtpwgts, &from, &cnum, parts);
      to = (from+1)%2;

      if (from == -1 || (higain = rpqGetTop(parts[from][cnum])) == -1)
        break;

      raxpy(ncon,  1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon,   1);
      raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);

      newcut -= (ed[higain]-id[higain]);
      newbal = Serial_Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);

      if ((newcut < mincut && newbal-origbal <= .00001) ||
          (newcut == mincut && (newbal < minbal ||
                                (newbal == minbal && 
                                 Serial_BetterBalance(ncon, npwgts, tpwgts, mindiff, tmpdiff))))) {
        mincut = newcut;
        minbal = newbal;
        mincutorder = nswaps;
        for (i=0; i<ncon; i++)
          mindiff[i] = fabs(tpwgts[i]-npwgts[i]);
      }
      else if (nswaps-mincutorder > limit) { /* We hit the limit, undo last move */
        newcut += (ed[higain]-id[higain]);
        raxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
        raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
        break;
      }

      where[higain] = to;
      moved[higain] = nswaps;
      swaps[nswaps] = higain;

      /**************************************************************
      * Update the id[i]/ed[i] values of the affected nodes
      ***************************************************************/
      gk_SWAP(id[higain], ed[higain], tmp);
      if (ed[higain] == 0 && xadj[higain] < xadj[higain+1])
        BNDDelete(nbnd, bndind,  bndptr, higain);

      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];

        kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
        INC_DEC(id[k], ed[k], kwgt);

        /* Update its boundary information and queue position */
        if (bndptr[k] != -1) { /* If k was a boundary vertex */
          if (ed[k] == 0) { /* Not a boundary vertex any more */
            BNDDelete(nbnd, bndind, bndptr, k);
            if (moved[k] == -1)  /* Remove it if in the queues */
              rpqDelete(parts[where[k]][qnum[k]], k);
          }
          else { /* If it has not been moved, update its position in the queue */
            if (moved[k] == -1)
              rpqUpdate(parts[where[k]][qnum[k]], k, (real_t)(ed[k]-id[k]));
          }
        }
        else {
          if (ed[k] > 0) {  /* It will now become a boundary vertex */
            BNDInsert(nbnd, bndind, bndptr, k);
            if (moved[k] == -1)
              rpqInsert(parts[where[k]][qnum[k]], k, (real_t)(ed[k]-id[k]));
          }
        }
      }
    }

    /****************************************************************
    * Roll back computations
    *****************************************************************/
    for (i=0; i<nswaps; i++)
      moved[swaps[i]] = -1;  /* reset moved array */
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];

      to = where[higain] = (where[higain]+1)%2;
     gk_SWAP(id[higain], ed[higain], tmp);
      if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
        BNDDelete(nbnd, bndind,  bndptr, higain);
      else if (ed[higain] > 0 && bndptr[higain] == -1)
        BNDInsert(nbnd, bndind,  bndptr, higain);

      raxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
      raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+((to+1)%2)*ncon, 1);
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];

        kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
        INC_DEC(id[k], ed[k], kwgt);

        if (bndptr[k] != -1 && ed[k] == 0)
          BNDDelete(nbnd, bndind, bndptr, k);
        if (bndptr[k] == -1 && ed[k] > 0)
          BNDInsert(nbnd, bndind, bndptr, k);
      }
    }

    graph->mincut = mincut;
    graph->gnvtxs = nbnd;

    if (mincutorder == -1 || mincut == initcut)
      break;
  }

  for (i=0; i<ncon; i++) {
    rpqDestroy(parts[0][i]);
    rpqDestroy(parts[1][i]);
  }

  WCOREPOP;

  return;
}


/*************************************************************************
* This function selects the partition number and the queue from which
* we will move vertices out
**************************************************************************/
void Serial_SelectQueue(idx_t ncon, real_t *npwgts, real_t *tpwgts, idx_t *from, 
         idx_t *cnum, rpq_t **queues[2])
{
  idx_t i, part;
  real_t maxgain=0.0;
  real_t max = -1.0, maxdiff=0.0;
  idx_t mype;
  gkMPI_Comm_rank(MPI_COMM_WORLD, &mype);

  *from = -1;
  *cnum = -1;

  /* First determine the side and the queue, irrespective of the presence of nodes */
  for (part=0; part<2; part++) {
    for (i=0; i<ncon; i++) {
      if (npwgts[part*ncon+i]-tpwgts[part*ncon+i] >= maxdiff) {
        maxdiff = npwgts[part*ncon+i]-tpwgts[part*ncon+i];
        *from = part;
        *cnum = i;
      }
    }
  }

  if (*from != -1 && rpqLength(queues[*from][*cnum]) == 0) {
    /* The desired queue is empty, select a node from that side anyway */
    for (i=0; i<ncon; i++) {
      if (rpqLength(queues[*from][i]) > 0) {
        max = npwgts[(*from)*ncon + i];
        *cnum = i;
        break;
      }
    }

    for (i++; i<ncon; i++) {
      if (npwgts[(*from)*ncon + i] > max && rpqLength(queues[*from][i]) > 0) {
        max = npwgts[(*from)*ncon + i];
        *cnum = i;
      }
    }
  }


  /* Check to see if you can focus on the cut */
  if (maxdiff <= 0.0 || *from == -1) {
    maxgain = -100000.0;

    for (part=0; part<2; part++) {
      for (i=0; i<ncon; i++) {
        if (rpqLength(queues[part][i]) > 0 &&
            rpqSeeTopKey(queues[part][i]) > maxgain) {
          maxgain = rpqSeeTopKey(queues[part][i]);
          *from = part;
          *cnum = i;
        }
      }
    }
  }

  return;
}


/*************************************************************************
* This function checks if the balance achieved is better than the diff
* For now, it uses a 2-norm measure
**************************************************************************/
idx_t Serial_BetterBalance(idx_t ncon, real_t *npwgts, real_t *tpwgts, 
          real_t *diff, real_t *tmpdiff)
{
  idx_t i;

  for (i=0; i<ncon; i++)
    tmpdiff[i] = fabs(tpwgts[i]-npwgts[i]);

  return rnorm2(ncon, tmpdiff, 1) < rnorm2(ncon, diff, 1);
}


/*************************************************************************
* This function computes the load imbalance over all the constrains
**************************************************************************/
real_t Serial_Compute2WayHLoadImbalance(idx_t ncon, real_t *npwgts, real_t *tpwgts)
{
  idx_t i;
  real_t max=0.0, temp;

  for (i=0; i<ncon; i++) {
    if (tpwgts[i] == 0.0)
      temp = 0.0;
    else
      temp = fabs(tpwgts[i]-npwgts[i])/tpwgts[i];
    max = (max < temp ? temp : max);
  }
  return 1.0+max;
}



/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
void Mc_Serial_Balance2Way(ctrl_t *ctrl, graph_t *graph, real_t *tpwgts, real_t lbfactor)
{
  idx_t i, ii, j, k, kwgt, nvtxs, ncon, nbnd, nswaps, from, to, limit, tmp, cnum;
  idx_t *xadj, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind;
  idx_t *moved, *swaps, *qnum;
  real_t *nvwgt, *npwgts, *mindiff, *tmpdiff, origbal, minbal, newbal;
  rpq_t **parts[2];
  idx_t higain, mincut, newcut, mincutorder;
  idx_t *qsizes[2];

  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  nvwgt  = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;
  id     = graph->sendind;
  ed     = graph->recvind;
  npwgts = graph->gnpwgts;
  bndptr = graph->sendptr;
  bndind = graph->recvptr;

  mindiff   = rwspacemalloc(ctrl, ncon);
  tmpdiff   = rwspacemalloc(ctrl, ncon);
  parts[0]  = (rpq_t **)wspacemalloc(ctrl, sizeof(rpq_t *)*ncon);
  parts[1]  = (rpq_t **)wspacemalloc(ctrl, sizeof(rpq_t *)*ncon);
  qsizes[0] = iset(ncon, 0, iwspacemalloc(ctrl, ncon));
  qsizes[1] = iset(ncon, 0, iwspacemalloc(ctrl, ncon));

  moved = iwspacemalloc(ctrl, nvtxs);
  swaps = iwspacemalloc(ctrl, nvtxs);
  qnum  = iwspacemalloc(ctrl, nvtxs);

  limit = gk_min(gk_max(0.01*nvtxs, 15), 100);

  /* Initialize the queues */
  for (i=0; i<ncon; i++) {
    parts[0][i] = rpqCreate(nvtxs);
    parts[1][i] = rpqCreate(nvtxs);
  }

  for (i=0; i<nvtxs; i++) {
    qnum[i] = rargmax(ncon, nvwgt+i*ncon);
    qsizes[qnum[i]][where[i]]++;
  }

  for (from=0; from<2; from++) {
    for (j=0; j<ncon; j++) {
      if (qsizes[j][from] == 0) {
        for (i=0; i<nvtxs; i++) {
          if (where[i] != from)
            continue;

          k = rargmax2(ncon, nvwgt+i*ncon);
          if (k == j &&
               qsizes[qnum[i]][from] > qsizes[j][from] &&
               nvwgt[i*ncon+qnum[i]] < 1.3*nvwgt[i*ncon+j]) {
            qsizes[qnum[i]][from]--;
            qsizes[j][from]++;
            qnum[i] = j;
          }
        }
      }
    }
  }


  for (i=0; i<ncon; i++)
    mindiff[i] = fabs(tpwgts[i]-npwgts[i]);
  minbal = origbal = Serial_Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);
  newcut = mincut = graph->mincut;
  mincutorder = -1;

  iset(nvtxs, -1, moved);

  /* Insert all nodes in the priority queues */
  nbnd = graph->gnvtxs;
  for (i=0; i<nvtxs; i++) 
    rpqInsert(parts[where[i]][qnum[i]], i, (real_t)(ed[i]-id[i]));

  for (nswaps=0; nswaps<nvtxs; nswaps++) {
    if (minbal < lbfactor)
      break;

    Serial_SelectQueue(ncon, npwgts, tpwgts, &from, &cnum, parts);
    to = (from+1)%2;

    if (from == -1 || (higain = rpqGetTop(parts[from][cnum])) == -1)
      break;

    raxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
    newcut -= (ed[higain]-id[higain]);
    newbal = Serial_Compute2WayHLoadImbalance(ncon, npwgts, tpwgts);

    if (newbal < minbal || (newbal == minbal &&
        (newcut < mincut || (newcut == mincut &&
          Serial_BetterBalance(ncon, npwgts, tpwgts, mindiff, tmpdiff))))) {
      mincut = newcut;
      minbal = newbal;
      mincutorder = nswaps;
      for (i=0; i<ncon; i++)
        mindiff[i] = fabs(tpwgts[i]-npwgts[i]);
    }
    else if (nswaps-mincutorder > limit) { /* We hit the limit, undo last move */
      newcut += (ed[higain]-id[higain]);
      raxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);
      raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
      break;
    }

    where[higain] = to;
    moved[higain] = nswaps;
    swaps[nswaps] = higain;

    /**************************************************************
    * Update the id[i]/ed[i] values of the affected nodes
    ***************************************************************/
    gk_SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      /* Update the queue position */
      if (moved[k] == -1)
        rpqUpdate(parts[where[k]][qnum[k]], k, (real_t)(ed[k]-id[k]));

      /* Update its boundary information */
      if (ed[k] == 0 && bndptr[k] != -1)
        BNDDelete(nbnd, bndind, bndptr, k);
      else if (ed[k] > 0 && bndptr[k] == -1)
        BNDInsert(nbnd, bndind, bndptr, k);
    }
  }


  /****************************************************************
  * Roll back computations
  *****************************************************************/
  for (nswaps--; nswaps>mincutorder; nswaps--) {
    higain = swaps[nswaps];

    to = where[higain] = (where[higain]+1)%2;
    gk_SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    else if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    raxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+((to+1)%2)*ncon, 1);
    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      if (bndptr[k] != -1 && ed[k] == 0)
        BNDDelete(nbnd, bndind, bndptr, k);
      if (bndptr[k] == -1 && ed[k] > 0)
        BNDInsert(nbnd, bndind, bndptr, k);
    }
  }

  graph->mincut = mincut;
  graph->gnvtxs = nbnd;


  for (i=0; i<ncon; i++) {
    rpqDestroy(parts[0][i]);
    rpqDestroy(parts[1][i]);
  }

  WCOREPOP;

  return;
}



/*************************************************************************
* This function balances two partitions by moving the highest gain
* (including negative gain) vertices to the other domain.
* It is used only when tha unbalance is due to non contigous
* subdomains. That is, the are no boundary vertices.
* It moves vertices from the domain that is overweight to the one that
* is underweight.
**************************************************************************/
void Mc_Serial_Init2WayBalance(ctrl_t *ctrl, graph_t *graph, real_t *tpwgts)
{
  idx_t i, ii, j, k;
  idx_t kwgt, nvtxs, nbnd, ncon, nswaps, from, to, cnum, tmp;
  idx_t *xadj, *adjncy, *adjwgt, *where, *id, *ed, *bndptr, *bndind;
  idx_t *qnum;
  real_t *nvwgt, *npwgts;
  rpq_t **parts[2];
  idx_t higain, mincut;

  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  nvwgt  = graph->nvwgt;
  adjwgt = graph->adjwgt;
  where  = graph->where;
  id     = graph->sendind;
  ed     = graph->recvind;
  npwgts = graph->gnpwgts;
  bndptr = graph->sendptr;
  bndind = graph->recvptr;

  parts[0] = (rpq_t **)wspacemalloc(ctrl, sizeof(rpq_t *)*ncon);
  parts[1] = (rpq_t **)wspacemalloc(ctrl, sizeof(rpq_t *)*ncon);

  qnum = iwspacemalloc(ctrl, nvtxs);

  /* This is called for initial partitioning so we know from where to pick nodes */
  from = 1;
  to = (from+1)%2;

  for (i=0; i<ncon; i++) {
    parts[0][i] = rpqCreate(nvtxs);
    parts[1][i] = rpqCreate(nvtxs);
  }

  /* Compute the queues in which each vertex will be assigned to */
  for (i=0; i<nvtxs; i++)
    qnum[i] = rargmax(ncon, nvwgt+i*ncon);

  /* Insert the nodes of the proper partition in the appropriate priority queue */
  for (i=0; i<nvtxs; i++) {
    if (where[i] == from) {
      if (ed[i] > 0)
        rpqInsert(parts[0][qnum[i]], i, (real_t)(ed[i]-id[i]));
      else
        rpqInsert(parts[1][qnum[i]], i, (real_t)(ed[i]-id[i]));
    }
  }

  mincut = graph->mincut;
  nbnd   = graph->gnvtxs;
  for (nswaps=0; nswaps<nvtxs; nswaps++) {
    if (Serial_AreAnyVwgtsBelow(ncon, 1.0, npwgts+from*ncon, 0.0, nvwgt, tpwgts+from*ncon))
      break;

    if ((cnum = Serial_SelectQueueOneWay(ncon, npwgts, tpwgts, from, parts)) == -1)
      break;


    if ((higain = rpqGetTop(parts[0][cnum])) == -1)
      higain = rpqGetTop(parts[1][cnum]);

    mincut -= (ed[higain]-id[higain]);
    raxpy(ncon, 1.0, nvwgt+higain*ncon, 1, npwgts+to*ncon, 1);
    raxpy(ncon, -1.0, nvwgt+higain*ncon, 1, npwgts+from*ncon, 1);

    where[higain] = to;

    /**************************************************************
    * Update the id[i]/ed[i] values of the affected nodes
    ***************************************************************/
    gk_SWAP(id[higain], ed[higain], tmp);
    if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
      BNDDelete(nbnd, bndind,  bndptr, higain);
    if (ed[higain] > 0 && bndptr[higain] == -1)
      BNDInsert(nbnd, bndind,  bndptr, higain);

    for (j=xadj[higain]; j<xadj[higain+1]; j++) {
      k = adjncy[j];

      kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]);
      INC_DEC(id[k], ed[k], kwgt);

      /* Update the queue position */
      if (where[k] == from) {
        if (ed[k] > 0 && bndptr[k] == -1) {  /* It moves in boundary */
          rpqDelete(parts[1][qnum[k]], k);
          rpqInsert(parts[0][qnum[k]], k, (real_t)(ed[k]-id[k]));
        }
        else { /* It must be in the boundary already */
          rpqUpdate(parts[0][qnum[k]], k, (real_t)(ed[k]-id[k]));
        }
      }

      /* Update its boundary information */
      if (ed[k] == 0 && bndptr[k] != -1)
        BNDDelete(nbnd, bndind, bndptr, k);
      else if (ed[k] > 0 && bndptr[k] == -1)
        BNDInsert(nbnd, bndind, bndptr, k);
    }
  }

  graph->mincut = mincut;
  graph->gnvtxs = nbnd;

  for (i=0; i<ncon; i++) {
    rpqDestroy(parts[0][i]);
    rpqDestroy(parts[1][i]);
  }

  WCOREPOP;
}


/*************************************************************************
* This function selects the partition number and the queue from which
* we will move vertices out
**************************************************************************/
idx_t Serial_SelectQueueOneWay(idx_t ncon, real_t *npwgts, real_t *tpwgts, 
          idx_t from, rpq_t **queues[2])
{
  idx_t i, cnum=-1;
  real_t max=0.0;

  for (i=0; i<ncon; i++) {
    if (npwgts[from*ncon+i]-tpwgts[from*ncon+i] >= max &&
        rpqLength(queues[0][i]) + rpqLength(queues[1][i]) > 0) {
      max = npwgts[from*ncon+i]-tpwgts[i];
      cnum = i;
    }
  }

  return cnum;
}


/*************************************************************************
* This function computes the initial id/ed
**************************************************************************/
void Mc_Serial_Compute2WayPartitionParams(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, me, nvtxs, ncon, nbnd, mincut;
  idx_t *xadj, *adjncy, *adjwgt;
  real_t *nvwgt, *npwgts;
  idx_t *id, *ed, *where;
  idx_t *bndptr, *bndind;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  nvwgt  = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  npwgts = rset(2*ncon, 0.0, graph->gnpwgts);
  id     = iset(nvtxs, 0, graph->sendind);
  ed     = iset(nvtxs, 0, graph->recvind);
  bndptr = iset(nvtxs, -1, graph->sendptr);
  bndind = graph->recvptr;

  /*------------------------------------------------------------
  / Compute now the id/ed degrees
  /------------------------------------------------------------*/
  nbnd = mincut = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    raxpy(ncon, 1.0, nvwgt+i*ncon, 1, npwgts+me*ncon, 1);

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[adjncy[j]])
        id[i] += adjwgt[j];
      else
        ed[i] += adjwgt[j];
    }

    if (ed[i] > 0 || xadj[i] == xadj[i+1]) {
      mincut += ed[i];
      BNDInsert(nbnd, bndind, bndptr, i);
    }
  }

  graph->mincut = mincut/2;
  graph->gnvtxs = nbnd;

}


/*************************************************************************
* This function checks if the vertex weights of two vertices are below
* a given set of values
**************************************************************************/
idx_t Serial_AreAnyVwgtsBelow(idx_t ncon, real_t alpha, real_t *vwgt1, real_t beta, real_t *vwgt2, real_t *limit)
{
  idx_t i;

  for (i=0; i<ncon; i++)
    if (alpha*vwgt1[i] + beta*vwgt2[i] < limit[i])
      return 1;

  return 0;
}


/*************************************************************************
*  This function computes the edge-cut of a serial graph.
**************************************************************************/
idx_t ComputeSerialEdgeCut(graph_t *graph)
{
  idx_t i, j;
  idx_t cut = 0;

  for (i=0; i<graph->nvtxs; i++) {
    for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
      if (graph->where[i] != graph->where[graph->adjncy[j]])
        cut += graph->adjwgt[j];
  }
  graph->mincut = cut/2;

  return graph->mincut;
}


/*************************************************************************
*  This function computes the TotalV of a serial graph.
**************************************************************************/
idx_t ComputeSerialTotalV(graph_t *graph, idx_t *home)
{
  idx_t i;
  idx_t totalv = 0;

  for (i=0; i<graph->nvtxs; i++)
    if (graph->where[i] != home[i])
      totalv += (graph->vsize == NULL) ? graph->vwgt[i] : graph->vsize[i];

  return totalv;
}


