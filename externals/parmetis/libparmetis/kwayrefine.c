/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * medge_refine.c
 *
 * This file contains code that performs the k-way refinement
 *
 * Started 3/1/96
 * George
 *
 * $Id: kwayrefine.c 10987 2011-10-31 14:42:33Z karypis $
 */

#include <parmetislib.h>

#define ProperSide(c, from, other) \
              (((c) == 0 && (from)-(other) < 0) || ((c) == 1 && (from)-(other) > 0))


/*************************************************************************
* This function projects a partition.
**************************************************************************/
void ProjectPartition(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, nvtxs, nnbrs = -1, firstvtx, cfirstvtx;
  idx_t *match, *cmap, *where, *cwhere;
  idx_t *peind, *slens = NULL, *rlens = NULL;
  ikv_t *rcand, *scand = NULL;
  graph_t *cgraph;

  WCOREPUSH;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ProjectTmr));

  cgraph    = graph->coarser;
  cwhere    = cgraph->where;
  cfirstvtx = cgraph->vtxdist[ctrl->mype];

  nvtxs    = graph->nvtxs;
  match    = graph->match;
  cmap     = graph->cmap;
  where    = graph->where = imalloc(nvtxs+graph->nrecv, "ProjectPartition: graph->where");
  firstvtx = graph->vtxdist[ctrl->mype];


  if (graph->match_type == PARMETIS_MTYPE_GLOBAL) {  /* Only if global matching is on */
    /*------------------------------------------------------------
    / Start the transmission of the remote where information 
    /------------------------------------------------------------*/
    nnbrs = graph->nnbrs;
    peind = graph->peind;
    slens = graph->slens;
    rlens = graph->rlens;
    rcand = graph->rcand;

    scand = ikvwspacemalloc(ctrl, slens[nnbrs]);

    /* Issue the receives first */
    for (i=0; i<nnbrs; i++) {
      if (slens[i+1]-slens[i] > 0) /* Issue a receive only if you are getting something */
        gkMPI_Irecv((void *)(scand+slens[i]), 2*(slens[i+1]-slens[i]), IDX_T, 
            peind[i], 1, ctrl->comm, ctrl->rreq+i);
    }

#ifdef DEBUG_PROJECT
    PrintPairs(ctrl, rlens[nnbrs], rcand, "rcand"); 
#endif

    /* Put the where[rcand[].key] into the val field */
    for (i=0; i<rlens[nnbrs]; i++) {
      PASSERT(ctrl, rcand[i].val >= 0 && rcand[i].val < cgraph->nvtxs);
      rcand[i].val = cwhere[rcand[i].val];
    }

#ifdef DEBUG_PROJECT
    PrintPairs(ctrl, rlens[nnbrs], rcand, "rcand");
    PrintVector(ctrl, nvtxs, firstvtx, cmap, "cmap");
#endif

    /* Issue the sends next */
    for (i=0; i<nnbrs; i++) {
      if (rlens[i+1]-rlens[i] > 0) /* Issue a send only if you are sending something */
        gkMPI_Isend((void *)(rcand+rlens[i]), 2*(rlens[i+1]-rlens[i]), IDX_T, 
            peind[i], 1, ctrl->comm, ctrl->sreq+i);
    }
  }

  /*------------------------------------------------------------
  / Project local vertices first
  /------------------------------------------------------------*/
  for (i=0; i<nvtxs; i++) {
    if (match[i] >= KEEP_BIT) {
      PASSERT(ctrl, cmap[i]-cfirstvtx>=0 && cmap[i]-cfirstvtx<cgraph->nvtxs);
      where[i] = cwhere[cmap[i]-cfirstvtx];
    }
  }

  if (graph->match_type == PARMETIS_MTYPE_GLOBAL) {  /* Only if global matching is on */
    /*------------------------------------------------------------
    / Wait for the nonblocking operations to finish
    /------------------------------------------------------------*/
    for (i=0; i<nnbrs; i++) {
      if (rlens[i+1]-rlens[i] > 0)  
        gkMPI_Wait(ctrl->sreq+i, &ctrl->status);
    }
    for (i=0; i<nnbrs; i++) {
      if (slens[i+1]-slens[i] > 0)  
        gkMPI_Wait(ctrl->rreq+i, &ctrl->status);
    }

#ifdef DEBUG_PROJECT
    PrintPairs(ctrl, slens[nnbrs], scand, "scand"); 
#endif

    /*------------------------------------------------------------
    / Project received vertices now
    /------------------------------------------------------------*/
    for (i=0; i<slens[nnbrs]; i++) {
      PASSERTP(ctrl, scand[i].key-firstvtx>=0 && scand[i].key-firstvtx<graph->nvtxs, 
          (ctrl, "%"PRIDX" %"PRIDX" %"PRIDX"\n", scand[i].key, firstvtx, graph->nvtxs));
      where[scand[i].key-firstvtx] = scand[i].val;
    }
  }

  FreeGraph(graph->coarser);
  graph->coarser = NULL;

  WCOREPOP;

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ProjectTmr));
}


/*************************************************************************
* This function computes the initial id/ed 
**************************************************************************/
void ComputePartitionParams(ctrl_t *ctrl, graph_t *graph)
{
  idx_t h, i, j, k, nvtxs, ncon, firstvtx, lastvtx;
  idx_t *xadj, *adjncy, *adjwgt, *vtxdist, *where;
  real_t *lnpwgts, *gnpwgts;
  ckrinfo_t *myrinfo;
  cnbr_t *mynbrs;
  idx_t me, other;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayInitTmr));

  nvtxs   = graph->nvtxs;
  ncon    = graph->ncon;
  vtxdist = graph->vtxdist;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  where   = graph->where;

  graph->ckrinfo = (ckrinfo_t *)gk_malloc(sizeof(ckrinfo_t)*nvtxs, "CPP: ckrinfo");
  memset(graph->ckrinfo, 0, sizeof(ckrinfo_t)*nvtxs);

  lnpwgts = graph->lnpwgts = rsmalloc(ctrl->nparts*ncon, 0.0, "CPP: lnpwgts");
  gnpwgts = graph->gnpwgts = rmalloc(ctrl->nparts*ncon, "CPP: gnpwgts");

  firstvtx = vtxdist[ctrl->mype];
  lastvtx  = vtxdist[ctrl->mype+1];

  /* Send/Receive the where information of interface vertices */
  CommInterfaceData(ctrl, graph, where, where+nvtxs); 

#ifdef DEBUG_COMPUTEPPARAM
  rprintf(ctrl, "where sum: %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"\n", 
      GlobalSESum(ctrl, nvtxs+graph->nrecv), 
      GlobalSESum(ctrl, isum(nvtxs+graph->nrecv, where, 1)),
      GlobalSESum(ctrl, isum(nvtxs+1, xadj, 1)),
      GlobalSESum(ctrl, isum(graph->nedges, adjncy, 1)),
      GlobalSESum(ctrl, isum(graph->nedges, adjwgt, 1))
      );

  PrintVector(ctrl, nvtxs+graph->nrecv, 0, where, "where");
#endif


  /* Compute now the id/ed degrees */
  graph->lmincut = 0;
  for (i=0; i<nvtxs; i++) {
    me      = where[i];
    myrinfo = graph->ckrinfo+i;

    for (h=0; h<ncon; h++)
      lnpwgts[me*ncon+h] += graph->nvwgt[i*ncon+h];

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[adjncy[j]])
        myrinfo->id += adjwgt[j];
      else
        myrinfo->ed += adjwgt[j];
    }


    if (myrinfo->ed > 0) {  /* Time to do some serious work */
      graph->lmincut += myrinfo->ed;

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
          PASSERT(ctrl, myrinfo->nnbrs <= xadj[i+1]-xadj[i]);
        }
      }
    }
    else {
      myrinfo->inbr = -1;
    }
  }

#ifdef DEBUG_COMPUTEPPARAM
  PrintVector(ctrl, ctrl->nparts*ncon, 0, lnpwgts, "lnpwgts");
#endif

  /* Finally, sum-up the partition weights */
  gkMPI_Allreduce((void *)lnpwgts, (void *)gnpwgts, ctrl->nparts*ncon, 
      REAL_T, MPI_SUM, ctrl->comm);

  graph->mincut = GlobalSESum(ctrl, graph->lmincut)/2;

#ifdef DEBUG_COMPUTEPPARAM
  PrintVector(ctrl, ctrl->nparts*ncon, 0, gnpwgts, "gnpwgts");
#endif

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayInitTmr));
}


/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void KWayFM(ctrl_t *ctrl, graph_t *graph, idx_t npasses)
{
  idx_t from, to = -1, oldto, oldcut, mydomain, yourdomain, imbalanced, overweight;
  idx_t h, i, ii, iii, j, k, c;
  idx_t pass, nvtxs, nedges, ncon;
  idx_t nmoves, nmoved, nswaps;
  idx_t me, firstvtx, lastvtx, yourlastvtx;
  idx_t npes = ctrl->npes, mype = ctrl->mype, nparts = ctrl->nparts;
  idx_t nlupd, nsupd, nnbrs, nchanged;
  idx_t *xadj, *adjncy, *adjwgt, *vtxdist;
  idx_t *where, *tmp_where, *moved, *oldEDs;
  real_t *lnpwgts, *gnpwgts, *ognpwgts, *pgnpwgts, *movewgts, *overfill;
  idx_t *update, *supdate, *rupdate, *pe_updates;
  idx_t *changed, *perm, *pperm, *htable;
  idx_t *peind, *recvptr, *sendptr;
  ikv_t *swchanges, *rwchanges;
  ckrinfo_t *myrinfo;
  cnbr_t *mynbrs;
  real_t *lbvec, *nvwgt, *badmaxpwgt, *ubvec, *tpwgts, lbavg, ubavg;
  idx_t *nupds_pe;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayTmr));
  WCOREPUSH;

  /*************************/
  /* set up common aliases */
  /*************************/
  nvtxs  = graph->nvtxs;
  nedges = graph->nedges;
  ncon   = graph->ncon;

  vtxdist = graph->vtxdist;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;

  firstvtx = vtxdist[mype];
  lastvtx  = vtxdist[mype+1];

  where   = graph->where;
  lnpwgts = graph->lnpwgts;
  gnpwgts = graph->gnpwgts;

  ubvec   = ctrl->ubvec;
  tpwgts  = ctrl->tpwgts;

  nnbrs   = graph->nnbrs;
  peind   = graph->peind;
  recvptr = graph->recvptr;
  sendptr = graph->sendptr;

  /************************************/
  /* set up important data structures */
  /************************************/
  lbvec        = rwspacemalloc(ctrl, ncon);
  badmaxpwgt   = rwspacemalloc(ctrl, nparts*ncon);
  movewgts     = rwspacemalloc(ctrl, nparts*ncon);
  ognpwgts     = rwspacemalloc(ctrl, nparts*ncon);
  pgnpwgts     = rwspacemalloc(ctrl, nparts*ncon);
  overfill     = rwspacemalloc(ctrl, nparts*ncon);

  pperm        = iwspacemalloc(ctrl, nparts);
  nupds_pe     = iwspacemalloc(ctrl, npes);

  oldEDs       = iwspacemalloc(ctrl, nvtxs);
  changed      = iwspacemalloc(ctrl, nvtxs);
  perm         = iwspacemalloc(ctrl, nvtxs);
  update       = iwspacemalloc(ctrl, nvtxs);
  moved        = iwspacemalloc(ctrl, nvtxs);
  htable       = iset(nvtxs+graph->nrecv, 0, iwspacemalloc(ctrl, nvtxs+graph->nrecv));

  rwchanges    = ikvwspacemalloc(ctrl, graph->nrecv);
  swchanges    = ikvwspacemalloc(ctrl, graph->nsend);
  supdate      = iwspacemalloc(ctrl, graph->nrecv);
  rupdate      = iwspacemalloc(ctrl, graph->nsend);

  tmp_where    = iwspacemalloc(ctrl, nvtxs+graph->nrecv);

  for (i=0; i<nparts; i++) {
    for (h=0; h<ncon; h++) 
      badmaxpwgt[i*ncon+h] = ubvec[h]*tpwgts[i*ncon+h];
  }

  icopy(nvtxs+graph->nrecv, where, tmp_where);

  /* this will record the overall external degrees of the vertices
     prior to a inner refinement iteration in order to allow for
     the proper updating of the lmincut */
  for (i=0; i<nvtxs; i++)
    oldEDs[i] = graph->ckrinfo[i].ed;


  /*********************************************************/
  /* perform a small number of passes through the vertices */
  /*********************************************************/
  for (nswaps=0, pass=0; pass<npasses; pass++) {
    if (mype == 0)
      RandomPermute(nparts, pperm, 1);

    gkMPI_Bcast((void *)pperm, nparts, IDX_T, 0, ctrl->comm);
    FastRandomPermute(nvtxs, perm, 1);
    oldcut = graph->mincut;

    /* check to see if the partitioning is imbalanced */
    ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
    ubavg = ravg(ncon, ubvec);
    lbavg = ravg(ncon, lbvec);
    imbalanced = (lbavg > ubavg) ? 1 : 0;

    for (c=0; c<2; c++) {
      rcopy(ncon*nparts, gnpwgts, ognpwgts);
      rset(ncon*nparts, 0.0, movewgts);
      nmoved = 0;

      /**********************************************/
      /* PASS ONE -- record stats for desired moves */
      /**********************************************/
      for (iii=0; iii<nvtxs; iii++) {
        i     = perm[iii];
        from  = tmp_where[i];
        nvwgt = graph->nvwgt+i*ncon;

        for (h=0; h<ncon; h++) {
          if (rabs(nvwgt[h]-gnpwgts[from*ncon+h]) < SMALLFLOAT)
            break;
        }

        if (h < ncon) 
          continue;

        /* check for a potential improvement */
        myrinfo = graph->ckrinfo + i;
        if (myrinfo->ed == 0 || myrinfo->ed < myrinfo->id) 
          continue;

        PASSERT(ctrl, myrinfo->inbr != -1);
        mynbrs = ctrl->cnbrpool + myrinfo->inbr;

        for (k=myrinfo->nnbrs-1; k>=0; k--) {
          to = mynbrs[k].pid;
          if (ProperSide(c, pperm[from], pperm[to])) {
            for (h=0; h<ncon; h++) {
              if (gnpwgts[to*ncon+h]+nvwgt[h] > badmaxpwgt[to*ncon+h] && nvwgt[h] > 0.0)
                break;
            }

            if (h == ncon)
              break;
          }
        }

        /* break out if you did not find a candidate */
        if (k < 0)
          continue;

        oldto = to;
        for (j=k-1; j>=0; j--) {
          to = mynbrs[j].pid;
          if (ProperSide(c, pperm[from], pperm[to])) {
            for (h=0; h<ncon; h++) {
              if (gnpwgts[to*ncon+h]+nvwgt[h] > badmaxpwgt[to*ncon+h] && nvwgt[h] > 0.0)
                break;
            }

            if (h == ncon) {
              if (mynbrs[j].ed > mynbrs[k].ed ||
                  (mynbrs[j].ed == mynbrs[k].ed &&
                   IsHBalanceBetterTT(ncon,gnpwgts+oldto*ncon,gnpwgts+to*ncon,nvwgt,ubvec))){
                oldto = to;
                k = j;
              }
            }
          }
        }
        to = oldto;

        if (mynbrs[k].ed > myrinfo->id ||
            (mynbrs[k].ed == myrinfo->id &&
            (imbalanced ||  graph->level > 3  || iii % 8 == 0) &&
            IsHBalanceBetterFT(ncon,gnpwgts+from*ncon,gnpwgts+to*ncon,nvwgt,ubvec))) {

          /****************************************/
          /* Update tmp arrays of the moved vertex */
          /****************************************/
          tmp_where[i] = to;
          moved[nmoved++] = i;
          for (h=0; h<ncon; h++) {
            INC_DEC(lnpwgts[to*ncon+h],  lnpwgts[from*ncon+h],  nvwgt[h]);
            INC_DEC(gnpwgts[to*ncon+h],  gnpwgts[from*ncon+h],  nvwgt[h]);
            INC_DEC(movewgts[to*ncon+h], movewgts[from*ncon+h], nvwgt[h]);
          }

          myrinfo->ed += myrinfo->id-mynbrs[k].ed;
          gk_SWAP(myrinfo->id, mynbrs[k].ed, j);
          if (mynbrs[k].ed == 0) 
            mynbrs[k] = mynbrs[--myrinfo->nnbrs];
          else 
            mynbrs[k].pid = from;

          /* Update the degrees of adjacent vertices */
          for (j=xadj[i]; j<xadj[i+1]; j++) {
            /* no need to bother about vertices on different pe's */
            if (adjncy[j] >= nvtxs)
              continue;

            me = adjncy[j];
            mydomain = tmp_where[me];

            myrinfo = graph->ckrinfo+me;
            if (myrinfo->inbr == -1) {
              myrinfo->inbr  = cnbrpoolGetNext(ctrl, xadj[me+1]-xadj[me]+1);
              myrinfo->nnbrs = 0;
            }
            mynbrs = ctrl->cnbrpool + myrinfo->inbr;

            if (mydomain == from) {
              INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);
            }
            else {
              if (mydomain == to) {
                INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);
              }
            }

            /* Remove contribution from the .ed of 'from' */
            if (mydomain != from) {
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

            /* Add contribution to the .ed of 'to' */
            if (mydomain != to) {
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
        }
      }

      /******************************************/
      /* Let processors know the subdomain wgts */
      /* if all proposed moves commit.          */
      /******************************************/
      gkMPI_Allreduce((void *)lnpwgts, (void *)pgnpwgts, nparts*ncon, REAL_T, 
          MPI_SUM, ctrl->comm);

      /**************************/
      /* compute overfill array */
      /**************************/
      overweight = 0;
      for (j=0; j<nparts; j++) {
        for (h=0; h<ncon; h++) {
          if (pgnpwgts[j*ncon+h] > ognpwgts[j*ncon+h]) {
            overfill[j*ncon+h] = (pgnpwgts[j*ncon+h]-badmaxpwgt[j*ncon+h]) /
                                 (pgnpwgts[j*ncon+h]-ognpwgts[j*ncon+h]);
          }
          else {
            overfill[j*ncon+h] = 0.0;
          }

          overfill[j*ncon+h]  = gk_max(overfill[j*ncon+h], 0.0);
          overfill[j*ncon+h] *= movewgts[j*ncon+h];

          if (overfill[j*ncon+h] > 0.0)
            overweight = 1;

          PASSERTP(ctrl, ognpwgts[j*ncon+h] <= badmaxpwgt[j*ncon+h] ||
                  pgnpwgts[j*ncon+h] <= ognpwgts[j*ncon+h],
                  (ctrl, "%.4"PRREAL" %.4"PRREAL" %.4"PRREAL"\n", ognpwgts[j*ncon+h],
                   badmaxpwgt[j*ncon+h], pgnpwgts[j*ncon+h]));
        }
      }

      /****************************************************/
      /* select moves to undo according to overfill array */
      /****************************************************/
      if (overweight == 1) {
        for (iii=0; iii<nmoved; iii++) {
          i     = moved[iii];
          oldto = tmp_where[i];
          nvwgt = graph->nvwgt+i*ncon;

          myrinfo = graph->ckrinfo + i;
          mynbrs  = ctrl->cnbrpool + myrinfo->inbr;
          PASSERT(ctrl, myrinfo->nnbrs == 0 || myrinfo->inbr != -1);

          for (k=0; k<myrinfo->nnbrs; k++) {
            if (mynbrs[k].pid == where[i])
              break;
          }

          for (h=0; h<ncon; h++) {
            if (nvwgt[h] > 0.0 && overfill[oldto*ncon+h] > nvwgt[h]/4.0)
              break;
          }

          /**********************************/
          /* nullify this move if necessary */
          /**********************************/
          if (k != myrinfo->nnbrs && h != ncon) {
            moved[iii] = -1;
            from       = oldto;
            to         = where[i];

            for (h=0; h<ncon; h++) 
              overfill[oldto*ncon+h] =gk_max(overfill[oldto*ncon+h]-nvwgt[h], 0.0);

            tmp_where[i] = to;
            myrinfo->ed += myrinfo->id-mynbrs[k].ed;
            gk_SWAP(myrinfo->id, mynbrs[k].ed, j);
            if (mynbrs[k].ed == 0) 
              mynbrs[k] = mynbrs[--myrinfo->nnbrs];
            else 
              mynbrs[k].pid = from;

            for (h=0; h<ncon; h++) {
              INC_DEC(lnpwgts[to*ncon+h],  lnpwgts[from*ncon+h],  nvwgt[h]);
            }

            /* Update the degrees of adjacent vertices */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              /* no need to bother about vertices on different pe's */
              if (adjncy[j] >= nvtxs)
                continue;

              me = adjncy[j];
              mydomain = tmp_where[me];

              myrinfo = graph->ckrinfo+me;
              if (myrinfo->inbr == -1) {
                myrinfo->inbr  = cnbrpoolGetNext(ctrl, xadj[me+1]-xadj[me]+1);
                myrinfo->nnbrs = 0;
              }
              mynbrs = ctrl->cnbrpool + myrinfo->inbr;

              if (mydomain == from) {
                INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);
              }
              else {
                if (mydomain == to) {
                  INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);
                }
              }

              /* Remove contribution from the .ed of 'from' */
              if (mydomain != from) {
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

              /* Add contribution to the .ed of 'to' */
              if (mydomain != to) {
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
          }
        }
      }


      /*************************************************/
      /* PASS TWO -- commit the remainder of the moves */
      /*************************************************/
      nlupd = nsupd = nmoves = nchanged = 0;
      for (iii=0; iii<nmoved; iii++) {
        if ((i = moved[iii]) == -1)
          continue;

        where[i] = tmp_where[i];

        /* Make sure to update the vertex information */
        if (htable[i] == 0) {
          /* make sure you do the update */
          htable[i] = 1;
          update[nlupd++] = i;
        }

        /* Put the vertices adjacent to i into the update array */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          if (htable[k] == 0) {
            htable[k] = 1;
            if (k<nvtxs)
              update[nlupd++] = k;
            else
              supdate[nsupd++] = k;
          }
        }
        nmoves++;
        nswaps++;

        if (graph->pexadj[i+1]-graph->pexadj[i] > 0)
          changed[nchanged++] = i;
      }

      /* Tell interested pe's the new where[] info for the interface vertices */
      CommChangedInterfaceData(ctrl, graph, nchanged, changed, where, 
          swchanges, rwchanges); 


      IFSET(ctrl->dbglvl, DBG_RMOVEINFO, rprintf(ctrl, 
            "\t[%"PRIDX" %"PRIDX"], [%.4"PRREAL"], [%"PRIDX" %"PRIDX" %"PRIDX"]\n", 
            pass, c, badmaxpwgt[0], GlobalSESum(ctrl, nmoves),
            GlobalSESum(ctrl, nsupd), GlobalSESum(ctrl, nlupd)));


      /*-------------------------------------------------------------
      / Time to communicate with processors to send the vertices
      / whose degrees need to be update.
      /-------------------------------------------------------------*/
      /* Issue the receives first */
      for (i=0; i<nnbrs; i++) {
        gkMPI_Irecv((void *)(rupdate+sendptr[i]), sendptr[i+1]-sendptr[i], IDX_T,
            peind[i], 1, ctrl->comm, ctrl->rreq+i);
      }

      /* Issue the sends next. This needs some preporcessing */
      for (i=0; i<nsupd; i++) {
        htable[supdate[i]] = 0;
        supdate[i] = graph->imap[supdate[i]];
      }
      isorti(nsupd, supdate);

      for (j=i=0; i<nnbrs; i++) {
        yourlastvtx = vtxdist[peind[i]+1];
        for (k=j; k<nsupd && supdate[k] < yourlastvtx; k++); 
        gkMPI_Isend((void *)(supdate+j), k-j, IDX_T, peind[i], 1, ctrl->comm, ctrl->sreq+i);
        j = k;
      }

      /* OK, now get into the loop waiting for the send/recv operations to finish */
      gkMPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
      for (i=0; i<nnbrs; i++) 
        gkMPI_Get_count(ctrl->statuses+i, IDX_T, nupds_pe+i);
      gkMPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


      /*-------------------------------------------------------------
      / Place the recieved to-be updated vertices into update[] 
      /-------------------------------------------------------------*/
      for (i=0; i<nnbrs; i++) {
        pe_updates = rupdate+sendptr[i];
        for (j=0; j<nupds_pe[i]; j++) {
          k = pe_updates[j];
          if (htable[k-firstvtx] == 0) {
            htable[k-firstvtx] = 1;
            update[nlupd++] = k-firstvtx;
          }
        }
      }


      /*-------------------------------------------------------------
      / Update the rinfo of the vertices in the update[] array
      /-------------------------------------------------------------*/
      for (ii=0; ii<nlupd; ii++) {
        i = update[ii];
        PASSERT(ctrl, htable[i] == 1);

        htable[i] = 0;

        mydomain = where[i];
        myrinfo  = graph->ckrinfo+i;

        if (myrinfo->inbr == -1) 
          myrinfo->inbr  = cnbrpoolGetNext(ctrl, xadj[i+1]-xadj[i]+1);
        mynbrs = ctrl->cnbrpool + myrinfo->inbr;

        graph->lmincut -= oldEDs[i];
        myrinfo->nnbrs  = 0;
        myrinfo->id     = 0;
        myrinfo->ed     = 0;

        for (j=xadj[i]; j<xadj[i+1]; j++) {
          yourdomain = where[adjncy[j]];
          if (mydomain != yourdomain) {
            myrinfo->ed += adjwgt[j];

            for (k=0; k<myrinfo->nnbrs; k++) {
              if (mynbrs[k].pid == yourdomain) {
                mynbrs[k].ed += adjwgt[j];
                break;
              }
            }
            if (k == myrinfo->nnbrs) {
              mynbrs[k].pid = yourdomain;
              mynbrs[k].ed  = adjwgt[j];
              myrinfo->nnbrs++;
            }
            PASSERT(ctrl, myrinfo->nnbrs <= xadj[i+1]-xadj[i]);
          }
          else {
            myrinfo->id += adjwgt[j];
          }
        }
        graph->lmincut += myrinfo->ed;
        oldEDs[i]       = myrinfo->ed;  /* for the next iteration */
      }

      /* finally, sum-up the partition weights */
      gkMPI_Allreduce((void *)lnpwgts, (void *)gnpwgts, nparts*ncon,
          REAL_T, MPI_SUM, ctrl->comm);
    }
    graph->mincut = GlobalSESum(ctrl, graph->lmincut)/2;

    IFSET(ctrl->dbglvl, DBG_RMOVEINFO, 
        rprintf(ctrl, "\t\tcut: %"PRIDX"\n", graph->mincut));

    if (graph->mincut == oldcut)
      break;
  }

  WCOREPOP;
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayTmr));
}


/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void KWayBalance(ctrl_t *ctrl, graph_t *graph, idx_t npasses)
{
  idx_t h, i, ii, iii, j, k, c;
  idx_t pass, nvtxs, nedges, ncon;
  idx_t nmoves, nmoved, nswaps;
  idx_t me, firstvtx, lastvtx, yourlastvtx;
  idx_t from, to = -1, oldto, oldcut, mydomain, yourdomain, imbalanced;
  idx_t npes = ctrl->npes, mype = ctrl->mype, nparts = ctrl->nparts;
  idx_t nlupd, nsupd, nnbrs, nchanged;
  idx_t *xadj, *adjncy, *adjwgt, *vtxdist;
  idx_t *where, *tmp_where, *moved, *oldEDs;
  real_t *lnpwgts, *gnpwgts;
  idx_t *update, *supdate, *rupdate, *pe_updates;
  idx_t *changed, *perm, *pperm, *htable;
  idx_t *peind, *recvptr, *sendptr;
  ikv_t *swchanges, *rwchanges;
  ckrinfo_t *myrinfo;
  cnbr_t *mynbrs;
  real_t *lbvec, *nvwgt, *badmaxpwgt, *ubvec, *tpwgts, lbavg, ubavg;
  idx_t *nupds_pe;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayTmr));
  WCOREPUSH;

  /*************************/
  /* set up common aliases */
  /*************************/
  nvtxs   = graph->nvtxs;
  nedges  = graph->nedges;
  ncon    = graph->ncon;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  where   = graph->where;
  lnpwgts = graph->lnpwgts;
  gnpwgts = graph->gnpwgts;

  nnbrs   = graph->nnbrs;
  peind   = graph->peind;
  recvptr = graph->recvptr;
  sendptr = graph->sendptr;

  vtxdist  = graph->vtxdist;
  firstvtx = vtxdist[mype];
  lastvtx  = vtxdist[mype+1];

  ubvec  = ctrl->ubvec;
  tpwgts = ctrl->tpwgts;


  /************************************/
  /* set up important data structures */
  /************************************/
  lbvec        = rwspacemalloc(ctrl, ncon);
  badmaxpwgt   = rwspacemalloc(ctrl, nparts*ncon);

  pperm        = iwspacemalloc(ctrl, nparts);
  nupds_pe     = iwspacemalloc(ctrl, npes);

  oldEDs       = iwspacemalloc(ctrl, nvtxs);
  changed      = iwspacemalloc(ctrl, nvtxs);
  perm         = iwspacemalloc(ctrl, nvtxs);
  update       = iwspacemalloc(ctrl, nvtxs);
  moved        = iwspacemalloc(ctrl, nvtxs);
  htable       = iset(nvtxs+graph->nrecv, 0, iwspacemalloc(ctrl, nvtxs+graph->nrecv));

  rwchanges    = ikvwspacemalloc(ctrl, graph->nrecv);
  swchanges    = ikvwspacemalloc(ctrl, graph->nsend);
  supdate      = iwspacemalloc(ctrl, graph->nrecv);
  rupdate      = iwspacemalloc(ctrl, graph->nsend);

  tmp_where    = iwspacemalloc(ctrl, nvtxs+graph->nrecv);

  for (i=0; i<nparts; i++) {
    for (h=0; h<ncon; h++) 
      badmaxpwgt[i*ncon+h] = ubvec[h]*tpwgts[i*ncon+h];
  }

  icopy(nvtxs+graph->nrecv, where, tmp_where);

  /* this will record the overall external degrees of the vertices
     prior to a inner refinement iteration in order to allow for
     the proper updating of the lmincut */
  for (i=0; i<nvtxs; i++)
    oldEDs[i] = graph->ckrinfo[i].ed;


  /*********************************************************/
  /* perform a small number of passes through the vertices */
  /*********************************************************/
  for (nswaps=0, pass=0; pass<npasses; pass++) {
    oldcut = graph->mincut;
    if (mype == 0)
      RandomPermute(nparts, pperm, 1);
    gkMPI_Bcast((void *)pperm, nparts, IDX_T, 0, ctrl->comm);
    FastRandomPermute(nvtxs, perm, 1);

    /* check to see if the partitioning is imbalanced */
    ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
    ubavg = ravg(ncon, ubvec);
    lbavg = ravg(ncon, lbvec);
    imbalanced = (lbavg > ubavg) ? 1 : 0;

    for (c=0; c<2; c++) {
      nmoved = 0;

      /**********************************************/
      /* PASS ONE -- record stats for desired moves */
      /**********************************************/
      for (iii=0; iii<nvtxs; iii++) {
        i = perm[iii];
        from = tmp_where[i];
        nvwgt = graph->nvwgt+i*ncon;

        for (h=0; h<ncon; h++) {
          if (fabs(nvwgt[h]-gnpwgts[from*ncon+h]) < SMALLFLOAT)
            break;
        }

        if (h < ncon) 
          continue;

        /* check for a potential improvement */
        myrinfo = graph->ckrinfo+i;
        if (myrinfo->ed == 0 || myrinfo->ed < myrinfo->id) 
          continue;

        PASSERT(ctrl, myrinfo->inbr != -1);
        mynbrs = ctrl->cnbrpool + myrinfo->inbr;

        for (k=myrinfo->nnbrs-1; k>=0; k--) {
          to = mynbrs[k].pid;
          if (ProperSide(c, pperm[from], pperm[to]) &&
              IsHBalanceBetterFT(ncon, gnpwgts+from*ncon, gnpwgts+to*ncon, nvwgt, ubvec)) {
            break;
          }
        }

        /* break out if you did not find a candidate */
        if (k < 0)
          continue;

        oldto = to;
        for (j=k-1; j>=0; j--) {
          to = mynbrs[j].pid;
          if (ProperSide(c, pperm[from], pperm[to]) &&
              IsHBalanceBetterTT(ncon, gnpwgts+oldto*ncon, gnpwgts+to*ncon, nvwgt, ubvec)){
            oldto = to;
            k = j;
          }
        }
        to = oldto;

        if (iii % npes == 0) {
          /****************************************/
          /* Update tmp arrays of the moved vertex */
          /****************************************/
          tmp_where[i] = to;
          moved[nmoved++] = i;
          for (h=0; h<ncon; h++) {
            INC_DEC(lnpwgts[to*ncon+h], lnpwgts[from*ncon+h], nvwgt[h]);
            INC_DEC(gnpwgts[to*ncon+h], gnpwgts[from*ncon+h], nvwgt[h]);
          }

          myrinfo->ed += myrinfo->id-mynbrs[k].ed;
          gk_SWAP(myrinfo->id, mynbrs[k].ed, j);
          if (mynbrs[k].ed == 0) 
            mynbrs[k] = mynbrs[--myrinfo->nnbrs];
          else
            mynbrs[k].pid = from;

          /* Update the degrees of adjacent vertices */
          for (j=xadj[i]; j<xadj[i+1]; j++) {
            /* no need to bother about vertices on different pe's */
            if (adjncy[j] >= nvtxs)
              continue;

            me = adjncy[j];
            mydomain = tmp_where[me];

            myrinfo = graph->ckrinfo+me;
            if (myrinfo->inbr == -1) {
              myrinfo->inbr  = cnbrpoolGetNext(ctrl, xadj[me+1]-xadj[me]+1);
              myrinfo->nnbrs = 0;
            }
            mynbrs = ctrl->cnbrpool + myrinfo->inbr;

            if (mydomain == from) {
              INC_DEC(myrinfo->ed, myrinfo->id, adjwgt[j]);
            }
            else {
              if (mydomain == to) {
                INC_DEC(myrinfo->id, myrinfo->ed, adjwgt[j]);
              }
            }

            /* Remove contribution from the .ed of 'from' */
            if (mydomain != from) {
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

            /* Add contribution to the .ed of 'to' */
            if (mydomain != to) {
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
        }
      }

      /*************************************************/
      /* PASS TWO -- commit the remainder of the moves */
      /*************************************************/
      nlupd = nsupd = nmoves = nchanged = 0;
      for (iii=0; iii<nmoved; iii++) {
        i = moved[iii];
        if (i == -1)
          continue;

        where[i] = tmp_where[i];

        /* Make sure to update the vertex information */
        if (htable[i] == 0) {
          /* make sure you do the update */
          htable[i] = 1;
          update[nlupd++] = i;
        }

        /* Put the vertices adjacent to i into the update array */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          if (htable[k] == 0) {
            htable[k] = 1;
            if (k<nvtxs)
              update[nlupd++] = k;
            else
              supdate[nsupd++] = k;
          }
        }
        nmoves++;
        nswaps++;

        if (graph->pexadj[i+1]-graph->pexadj[i] > 0)
          changed[nchanged++] = i;
      }

      /* Tell interested pe's the new where[] info for the interface vertices */
      CommChangedInterfaceData(ctrl, graph, nchanged, changed, where, swchanges, 
          rwchanges); 

      IFSET(ctrl->dbglvl, DBG_RMOVEINFO, rprintf(ctrl, 
            "\t[%"PRIDX" %"PRIDX"], [%.4"PRREAL"],  [%"PRIDX" %"PRIDX" %"PRIDX"]\n",
            pass, c, badmaxpwgt[0], GlobalSESum(ctrl, nmoves),
            GlobalSESum(ctrl, nsupd), GlobalSESum(ctrl, nlupd)));

      /*-------------------------------------------------------------
      / Time to communicate with processors to send the vertices
      / whose degrees need to be update.
      /-------------------------------------------------------------*/
      /* Issue the receives first */
      for (i=0; i<nnbrs; i++) {
        gkMPI_Irecv((void *)(rupdate+sendptr[i]), sendptr[i+1]-sendptr[i], IDX_T,
            peind[i], 1, ctrl->comm, ctrl->rreq+i);
      }

      /* Issue the sends next. This needs some preporcessing */
      for (i=0; i<nsupd; i++) {
        htable[supdate[i]] = 0;
        supdate[i] = graph->imap[supdate[i]];
      }
      isorti(nsupd, supdate);

      for (j=i=0; i<nnbrs; i++) {
        yourlastvtx = vtxdist[peind[i]+1];
        for (k=j; k<nsupd && supdate[k] < yourlastvtx; k++); 
        gkMPI_Isend((void *)(supdate+j), k-j, IDX_T, peind[i], 1, ctrl->comm, ctrl->sreq+i);
        j = k;
      }

      /* OK, now get into the loop waiting for the send/recv operations to finish */
      gkMPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
      for (i=0; i<nnbrs; i++) 
        gkMPI_Get_count(ctrl->statuses+i, IDX_T, nupds_pe+i);
      gkMPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


      /*-------------------------------------------------------------
      / Place the recieved to-be updated vertices into update[] 
      /-------------------------------------------------------------*/
      for (i=0; i<nnbrs; i++) {
        pe_updates = rupdate+sendptr[i];
        for (j=0; j<nupds_pe[i]; j++) {
          k = pe_updates[j];
          if (htable[k-firstvtx] == 0) {
            htable[k-firstvtx] = 1;
            update[nlupd++] = k-firstvtx;
          }
        }
      }


      /*-------------------------------------------------------------
      / Update the rinfo of the vertices in the update[] array
      /-------------------------------------------------------------*/
      for (ii=0; ii<nlupd; ii++) {
        i = update[ii];
        PASSERT(ctrl, htable[i] == 1);

        htable[i] = 0;

        mydomain = where[i];
        myrinfo  = graph->ckrinfo+i;
        if (myrinfo->inbr == -1)
          myrinfo->inbr  = cnbrpoolGetNext(ctrl, xadj[i+1]-xadj[i]+1);
        mynbrs = ctrl->cnbrpool + myrinfo->inbr;

        graph->lmincut -= oldEDs[i];
        myrinfo->nnbrs  = 0;
        myrinfo->id     = 0;
        myrinfo->ed     = 0;

        for (j=xadj[i]; j<xadj[i+1]; j++) {
          yourdomain = where[adjncy[j]];
          if (mydomain != yourdomain) {
            myrinfo->ed += adjwgt[j];

            for (k=0; k<myrinfo->nnbrs; k++) {
              if (mynbrs[k].pid == yourdomain) {
                mynbrs[k].ed += adjwgt[j];
                break;
              }
            }
            if (k == myrinfo->nnbrs) {
              mynbrs[k].pid = yourdomain;
              mynbrs[k].ed  = adjwgt[j];
              myrinfo->nnbrs++;
            }
            PASSERT(ctrl, myrinfo->nnbrs <= xadj[i+1]-xadj[i]);
          }
          else {
            myrinfo->id += adjwgt[j];
          }
        }
        graph->lmincut += myrinfo->ed;
        oldEDs[i]       = myrinfo->ed;  /* for the next iteration */
      }

      /* finally, sum-up the partition weights */
      gkMPI_Allreduce((void *)lnpwgts, (void *)gnpwgts, nparts*ncon,
          REAL_T, MPI_SUM, ctrl->comm);
    }
    graph->mincut = GlobalSESum(ctrl, graph->lmincut)/2;

    if (graph->mincut == oldcut)
      break;
  }

/*
  gnswaps = GlobalSESum(ctrl, nswaps);
  if (mype == 0)
    printf("niters: %"PRIDX", nswaps: %"PRIDX"\n", pass+1, gnswaps);
*/

  WCOREPOP;
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayTmr));
}


