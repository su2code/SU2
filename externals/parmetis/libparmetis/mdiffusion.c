/* * Copyright 1997, Regents of the University of Minnesota
 *
 * mdiffusion.c
 *
 * This file contains code that performs mc-diffusion
 *
 * Started 9/16/99
 * George
 *
 * $Id: mdiffusion.c 10542 2011-07-11 16:56:22Z karypis $
 */

#include <parmetislib.h>

#define PE	-1

/*************************************************************************
* This function is the entry point of the initial partitioning algorithm.
* This algorithm assembles the graph to all the processors and preceed
* serially.
**************************************************************************/
idx_t Mc_Diffusion(ctrl_t *ctrl, graph_t *graph, idx_t *vtxdist, idx_t *where, 
          idx_t *home, idx_t npasses)
{
  idx_t h, i, j;
  idx_t nvtxs, nedges, ncon, pass, iter, domain, processor;
  idx_t nparts, mype, npes, nlinks, me, you, wsize;
  idx_t nvisited, nswaps = -1, tnswaps, done, alldone = -1;
  idx_t *rowptr, *colind, *diff_where, *sr_where, *ehome, *map, *rmap;
  idx_t *pack, *unpack, *match, *proc2sub, *sub2proc;
  idx_t *visited, *gvisited;
  real_t *transfer, *npwgts, maxdiff, minflow, maxflow;
  real_t lbavg, oldlbavg, ubavg, *lbvec;
  real_t *diff_flows, *sr_flows;
  real_t diff_lbavg, sr_lbavg, diff_cost, sr_cost;
  idx_t *rbuffer, *sbuffer; 
  idx_t *rcount, *rdispl;
  real_t *solution, *load, *workspace;
  matrix_t matrix;
  graph_t *egraph;

  if (graph->ncon > 3)
    return 0;

  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  nedges = graph->nedges;
  ncon   = graph->ncon;

  nparts = ctrl->nparts;
  mype   = ctrl->mype;
  npes   = ctrl->npes;
  ubavg  = ravg(ncon, ctrl->ubvec);

  /* initialize variables and allocate memory */
  lbvec      = rwspacemalloc(ctrl, ncon);
  diff_flows = rwspacemalloc(ctrl, ncon);
  sr_flows   = rwspacemalloc(ctrl, ncon);

  load                       = rwspacemalloc(ctrl, nparts);
  solution                   = rwspacemalloc(ctrl, nparts);
  npwgts = graph->gnpwgts    = rwspacemalloc(ctrl, ncon*nparts);
  matrix.values              = rwspacemalloc(ctrl, nedges);
  transfer = matrix.transfer = rwspacemalloc(ctrl, ncon*nedges);

  proc2sub               = iwspacemalloc(ctrl, gk_max(nparts, npes*2));
  sub2proc               = iwspacemalloc(ctrl, nparts);
  match                  = iwspacemalloc(ctrl, nparts);
  rowptr = matrix.rowptr = iwspacemalloc(ctrl, nparts+1);
  colind = matrix.colind = iwspacemalloc(ctrl, nedges);

  rcount = iwspacemalloc(ctrl, npes);
  rdispl = iwspacemalloc(ctrl, npes+1);

  pack       = iwspacemalloc(ctrl, nvtxs);
  unpack     = iwspacemalloc(ctrl, nvtxs);
  rbuffer    = iwspacemalloc(ctrl, nvtxs);
  sbuffer    = iwspacemalloc(ctrl, nvtxs);
  map        = iwspacemalloc(ctrl, nvtxs);
  rmap       = iwspacemalloc(ctrl, nvtxs);
  diff_where = iwspacemalloc(ctrl, nvtxs);
  ehome      = iwspacemalloc(ctrl, nvtxs);


  wsize = gk_max(sizeof(real_t)*nparts*6, sizeof(idx_t)*(nvtxs+nparts*2+1));
  workspace = (real_t *)gk_malloc(wsize, "Mc_Diffusion: workspace");

  graph->ckrinfo = (ckrinfo_t *)gk_malloc(nvtxs*sizeof(ckrinfo_t), "Mc_Diffusion: rinfo");


  /* construct subdomain connectivity graph */
  matrix.nrows = nparts;
  SetUpConnectGraph(graph, &matrix, (idx_t *)workspace);
  nlinks = (matrix.nnzs-nparts) / 2;

  visited  = iwspacemalloc(ctrl, matrix.nnzs);
  gvisited = iwspacemalloc(ctrl, matrix.nnzs);

  for (pass=0; pass<npasses; pass++) {
    rset(matrix.nnzs*ncon, 0.0, transfer);
    iset(matrix.nnzs, 0, gvisited);
    iset(matrix.nnzs, 0, visited);
    iter = nvisited = 0;

    /* compute ncon flow solutions */
    for (h=0; h<ncon; h++) {
      rset(nparts, 0.0, solution);
      ComputeLoad(graph, nparts, load, ctrl->tpwgts, h);

      lbvec[h] = (rmax(nparts, load)+1.0/nparts) * (real_t)nparts;

      ConjGrad2(&matrix, load, solution, 0.001, workspace);
      ComputeTransferVector(ncon, &matrix, solution, transfer, h);
    }

    oldlbavg = ravg(ncon, lbvec);
    tnswaps = 0;
    maxdiff = 0.0;
    for (i=0; i<nparts; i++) {
      for (j=rowptr[i]; j<rowptr[i+1]; j++) {
        maxflow = rmax(ncon, transfer+j*ncon);
        minflow = rmin(ncon, transfer+j*ncon);
        maxdiff = (maxflow - minflow > maxdiff) ? maxflow - minflow : maxdiff;
      }
    }

    while (nvisited < nlinks) {
      /* compute independent sets of subdomains */
      iset(gk_max(nparts, npes*2), UNMATCHED, proc2sub);
      CSR_Match_SHEM(&matrix, match, proc2sub, gvisited, ncon);

      /* set up the packing arrays */
      iset(nparts, UNMATCHED, sub2proc);
      for (i=0; i<npes*2; i++) {
        if (proc2sub[i] == UNMATCHED)
          break;

        sub2proc[proc2sub[i]] = i/2;
      }

      iset(npes, 0, rcount);
      for (i=0; i<nvtxs; i++) {
        domain = where[i];
        processor = sub2proc[domain];
        if (processor != UNMATCHED) 
          rcount[processor]++;
      }

      rdispl[0] = 0;
      for (i=1; i<npes+1; i++)
        rdispl[i] = rdispl[i-1] + rcount[i-1];

      iset(nvtxs, UNMATCHED, unpack);
      for (i=0; i<nvtxs; i++) {
        domain = where[i];
        processor = sub2proc[domain];
        if (processor != UNMATCHED) 
          unpack[rdispl[processor]++] = i;
      }

      SHIFTCSR(i, npes, rdispl);

      iset(nvtxs, UNMATCHED, pack);
      for (i=0; i<rdispl[npes]; i++) {
        ASSERT(unpack[i] != UNMATCHED);
        domain = where[unpack[i]];
        processor = sub2proc[domain];
        if (processor != UNMATCHED) 
          pack[unpack[i]] = i;
      }

      /* Compute the flows */
      if (proc2sub[mype*2] != UNMATCHED) {
        me  = proc2sub[2*mype];
        you = proc2sub[2*mype+1];
        ASSERT(me != you);

        for (j=rowptr[me]; j<rowptr[me+1]; j++) {
          if (colind[j] == you) {
            visited[j] = 1;
            rcopy(ncon, transfer+j*ncon, diff_flows);
            break;
          }
        }

        for (j=rowptr[you]; j<rowptr[you+1]; j++) {
          if (colind[j] == me) {
            visited[j] = 1;
            for (h=0; h<ncon; h++) {
              if (transfer[j*ncon+h] > 0.0)
                diff_flows[h] = -1.0 * transfer[j*ncon+h];
            }
            break;
          }
        } 

        nswaps = 1;
        rcopy(ncon, diff_flows, sr_flows);

        iset(nvtxs, 0, sbuffer);
        for (i=0; i<nvtxs; i++) {
          if (where[i] == me || where[i] == you)
            sbuffer[i] = 1;
        }

        egraph = ExtractGraph(ctrl, graph, sbuffer, map, rmap);

        if (egraph != NULL) {
          icopy(egraph->nvtxs, egraph->where, diff_where);
          for (j=0; j<egraph->nvtxs; j++)
            ehome[j] = home[map[j]];
 
          RedoMyLink(ctrl, egraph, ehome, me, you, sr_flows, &sr_cost, &sr_lbavg);

          if (ncon <= 4) {
            sr_where      = egraph->where;
            egraph->where = diff_where;

            nswaps = BalanceMyLink(ctrl, egraph, ehome, me, you, diff_flows, maxdiff, 
                         &diff_cost, &diff_lbavg, 1.0/(real_t)nvtxs);

            if ((sr_lbavg < diff_lbavg &&
                (diff_lbavg >= ubavg-1.0 || sr_cost == diff_cost)) ||
                (sr_lbavg < ubavg-1.0 && sr_cost < diff_cost)) {
              for (i=0; i<egraph->nvtxs; i++)
                where[map[i]] = sr_where[i];
            }
            else {
              for (i=0; i<egraph->nvtxs; i++)
                where[map[i]] = diff_where[i];
            }
          }
          else {
            for (i=0; i<egraph->nvtxs; i++)
              where[map[i]] = egraph->where[i];
          }

          gk_free((void **)&egraph->xadj, &egraph->nvwgt, &egraph->adjncy, &egraph, LTERM);
        }

        /* Pack the flow data */
        iset(nvtxs, UNMATCHED, sbuffer);
        for (i=0; i<nvtxs; i++) {
          domain = where[i];
          if (domain == you || domain == me) 
            sbuffer[pack[i]] = where[i];
        }
      }

      /* Broadcast the flow data */
      gkMPI_Allgatherv((void *)&sbuffer[rdispl[mype]], rcount[mype], IDX_T, 
          (void *)rbuffer, rcount, rdispl, IDX_T, ctrl->comm);

      /* Unpack the flow data */
      for (i=0; i<rdispl[npes]; i++) {
        if (rbuffer[i] != UNMATCHED) 
          where[unpack[i]] = rbuffer[i];
      }


      /* Do other stuff */
      gkMPI_Allreduce((void *)visited, (void *)gvisited, matrix.nnzs,
          IDX_T, MPI_MAX, ctrl->comm);
      nvisited = isum(matrix.nnzs, gvisited, 1)/2;
      tnswaps += GlobalSESum(ctrl, nswaps);

      if (iter++ == NGD_PASSES)
        break;
    }

    /* perform serial refinement */
    Mc_ComputeSerialPartitionParams(ctrl, graph, nparts);
    Mc_SerialKWayAdaptRefine(ctrl, graph, nparts, home, ctrl->ubvec, 10);

    /* check for early breakout */
    for (h=0; h<ncon; h++) {
      lbvec[h] = (real_t)(nparts) *
        npwgts[rargmax_strd(nparts,npwgts+h,ncon)*ncon+h];
    }
    lbavg = ravg(ncon, lbvec);

    done = 0;
    if (tnswaps == 0 || lbavg >= oldlbavg || lbavg <= ubavg + 0.035)
      done = 1;

    alldone = GlobalSEMax(ctrl, done);
    if (alldone == 1)
      break;
  }

  /* ensure that all subdomains have at least one vertex */
/*
  iset(nparts, 0, match);
  for (i=0; i<nvtxs; i++)
    match[where[i]]++;

  done = 0;
  while (done == 0) {
    done = 1;

    me = iargmin(nparts, match);  
    if (match[me] == 0) {
      if (ctrl->mype == PE) printf("WARNING: empty subdomain %"PRIDX" in Mc_Diffusion\n", me);
      you = iargmax(nparts, match);  
      for (i=0; i<nvtxs; i++) {
        if (where[i] == you) {
          where[i] = me;
          match[you]--;
          match[me]++;
          done = 0;
          break;
        }
      }
    }
  }
*/
 
  /* now free memory and return */
  gk_free((void **)&workspace, (void **)&graph->ckrinfo, LTERM);
  graph->gnpwgts = NULL;
  graph->ckrinfo = NULL;

  WCOREPOP;

  return 0;
}


/*************************************************************************
* This function extracts a subgraph from a graph given an indicator array.
**************************************************************************/
graph_t *ExtractGraph(ctrl_t *ctrl, graph_t *graph, idx_t *indicator,
            idx_t *map, idx_t *rmap)
{
  idx_t h, i, j;
  idx_t nvtxs, envtxs, enedges, ncon;
  idx_t vtx, count;
  idx_t *xadj, *vsize, *adjncy, *adjwgt, *where;
  idx_t *exadj, *evsize, *eadjncy, *eadjwgt, *ewhere;
  real_t *nvwgt, *envwgt;
  graph_t *egraph;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  nvwgt  = graph->nvwgt;
  vsize  = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  count = 0;
  for (i=0; i<nvtxs; i++) {
    if (indicator[i] == 1) {
      map[count] = i;
      rmap[i] = count;
      count++;
    }
  }

  if (count == 0) 
    return NULL;

  /*******************/
  /* allocate memory */
  /*******************/
  egraph = CreateGraph();
  envtxs = egraph->nvtxs = count;
  egraph->ncon = graph->ncon;

  exadj  = egraph->xadj  = imalloc(envtxs*3+1, "exadj");
  ewhere = egraph->where = exadj + envtxs + 1;
  evsize = egraph->vsize = exadj + 2*envtxs + 1;

  envwgt = egraph->nvwgt = rmalloc(envtxs*ncon, "envwgt");

  /************************************************/
  /* compute xadj, where, nvwgt, and vsize arrays */
  /************************************************/
  iset(envtxs+1, 0, exadj);
  for (i=0; i<envtxs; i++) {
    vtx = map[i];

    ewhere[i] = where[vtx];
    for (h=0; h<ncon; h++)
      envwgt[i*ncon+h] = nvwgt[vtx*ncon+h];

    if (ctrl->partType == ADAPTIVE_PARTITION || ctrl->partType == REFINE_PARTITION) 
      evsize[i] = vsize[vtx];

    for (j=xadj[vtx]; j<xadj[vtx+1]; j++)
      if (indicator[adjncy[j]] == 1)
        exadj[i]++;

  }
  MAKECSR(i, envtxs, exadj);

  /************************************/
  /* compute adjncy and adjwgt arrays */
  /************************************/
  enedges = egraph->nedges = exadj[envtxs];
  eadjncy = egraph->adjncy = imalloc(enedges*2, "eadjncy");
  eadjwgt = egraph->adjwgt = eadjncy + enedges;

  for (i=0; i<envtxs; i++) {
    vtx = map[i];
    for (j=xadj[vtx]; j<xadj[vtx+1]; j++) {
      if (indicator[adjncy[j]] == 1) {
        eadjncy[exadj[i]] = rmap[adjncy[j]];
        eadjwgt[exadj[i]++] = adjwgt[j];
      }
    }
  }

  for (i=envtxs; i>0; i--)
    exadj[i] = exadj[i-1];
  exadj[0] = 0;

  return egraph;
}
