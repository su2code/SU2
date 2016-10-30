/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * wavefrontK.c 
 *
 * This file contains code for the initial directed diffusion at the coarsest
 * graph
 *
 * Started 5/19/97, Kirk, George
 *
 * $Id: diffutil.c 13946 2013-03-30 15:51:45Z karypis $
 *
 */

#include <parmetislib.h>


/*************************************************************************
*  This function computes the load for each subdomain
**************************************************************************/
void SetUpConnectGraph(graph_t *graph, matrix_t *matrix, idx_t *workspace)
{
  idx_t i, ii, j, jj, k, l;
  idx_t nvtxs, nrows;
  idx_t *xadj, *adjncy, *where;
  idx_t *rowptr, *colind;
  idx_t *pcounts, *perm, *marker;
  real_t *values;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  where  = graph->where;

  nrows  = matrix->nrows;
  rowptr = matrix->rowptr;
  colind = matrix->colind;
  values = matrix->values;

  perm    = workspace;
  marker  = iset(nrows, -1, workspace+nvtxs);
  pcounts = iset(nrows+1, 0, workspace+nvtxs+nrows);

  for (i=0; i<nvtxs; i++)
    pcounts[where[i]]++;
  MAKECSR(i, nrows, pcounts);

  for (i=0; i<nvtxs; i++)
    perm[pcounts[where[i]]++] = i;

  for (i=nrows; i>0; i--)
    pcounts[i] = pcounts[i-1];
  pcounts[0] = 0;

  /************************/
  /* Construct the matrix */
  /************************/
  rowptr[0] = k = 0;
  for (ii=0; ii<nrows; ii++) {
    colind[k++] = ii;
    marker[ii] = ii;

    for (jj=pcounts[ii]; jj<pcounts[ii+1]; jj++) {
      i = perm[jj];
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        l = where[adjncy[j]];
        if (marker[l] != ii) {
          colind[k] = l;
          values[k++] = -1.0;
          marker[l] = ii;
        }
      }
    }
    values[rowptr[ii]] = (real_t)(k-rowptr[ii]-1);
    rowptr[ii+1] = k;
  }
  matrix->nnzs = rowptr[nrows];

  return;
}


/*************************************************************************
* This function computes movement statistics for adaptive refinement
* schemes
**************************************************************************/
void Mc_ComputeMoveStatistics(ctrl_t *ctrl, graph_t *graph, idx_t *nmoved, idx_t *maxin, idx_t *maxout)
{
  idx_t i, nvtxs, nparts, myhome;
  idx_t *vwgt, *where;
  idx_t *lend, *gend, *lleft, *gleft, *lstart, *gstart;

  nvtxs = graph->nvtxs;
  vwgt = graph->vwgt;
  where = graph->where;
  nparts = ctrl->nparts;

  lstart = ismalloc(nparts, 0, "ComputeMoveStatistics: lstart");
  gstart = ismalloc(nparts, 0, "ComputeMoveStatistics: gstart");
  lleft = ismalloc(nparts, 0, "ComputeMoveStatistics: lleft");
  gleft = ismalloc(nparts, 0, "ComputeMoveStatistics: gleft");
  lend = ismalloc(nparts, 0, "ComputeMoveStatistics: lend");
  gend = ismalloc(nparts, 0, "ComputeMoveStatistics: gend");

  for (i=0; i<nvtxs; i++) {
    myhome = (ctrl->ps_relation == PARMETIS_PSR_COUPLED) ? ctrl->mype : graph->home[i];
    lstart[myhome] += (graph->vsize == NULL) ? 1 : graph->vsize[i];
    lend[where[i]] += (graph->vsize == NULL) ? 1 : graph->vsize[i];
    if (where[i] != myhome)
      lleft[myhome] += (graph->vsize == NULL) ? 1 : graph->vsize[i];
  }

  /* PrintVector(ctrl, ctrl->npes, 0, lend, "Lend: "); */

  gkMPI_Allreduce((void *)lstart, (void *)gstart, nparts, IDX_T, MPI_SUM, ctrl->comm);
  gkMPI_Allreduce((void *)lleft, (void *)gleft, nparts, IDX_T, MPI_SUM, ctrl->comm);
  gkMPI_Allreduce((void *)lend, (void *)gend, nparts, IDX_T, MPI_SUM, ctrl->comm);

  *nmoved = isum(nparts, gleft, 1);
  *maxout = imax(nparts, gleft);
  for (i=0; i<nparts; i++)
    lstart[i] = gend[i]+gleft[i]-gstart[i];
  *maxin = imax(nparts, lstart);

  gk_free((void **)&lstart, (void **)&gstart, (void **)&lleft, (void **)&gleft, (void **)&lend, (void **)&gend, LTERM);
}


/*************************************************************************
*  This function computes the TotalV of a serial graph.
**************************************************************************/
idx_t Mc_ComputeSerialTotalV(graph_t *graph, idx_t *home)
{
  idx_t i;
  idx_t totalv = 0;

  for (i=0; i<graph->nvtxs; i++) {
    if (graph->where[i] != home[i])
      totalv += (graph->vsize == NULL) ? graph->vwgt[i*graph->ncon] : graph->vsize[i];
  }

  return totalv;
}


/*************************************************************************
*  This function computes the load for each subdomain
**************************************************************************/
void ComputeLoad(graph_t *graph, idx_t nparts, real_t *load, real_t *tpwgts, idx_t index)
{
  idx_t i;
  idx_t nvtxs, ncon;
  idx_t *where;
  real_t *nvwgt;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  where = graph->where;
  nvwgt = graph->nvwgt;

  rset(nparts, 0.0, load);

  for (i=0; i<nvtxs; i++)
    load[where[i]] += nvwgt[i*ncon+index];

  ASSERT(fabs(rsum(nparts, load, 1)-1.0) < 0.001);

  for (i=0; i<nparts; i++) {
    load[i] -= tpwgts[i*ncon+index];
  }

  return;
}


/*************************************************************************
* This function implements the CG solver used during the directed diffusion
**************************************************************************/
void ConjGrad2(matrix_t *A, real_t *b, real_t *x, real_t tol, real_t *workspace)
{
  idx_t i, k, n;
  real_t *p, *r, *q, *z, *M;
  real_t alpha, beta, rho, rho_1 = -1.0, error, bnrm2, tmp;
  idx_t *rowptr, *colind;
  real_t *values;

  n      = A->nrows;
  rowptr = A->rowptr;
  colind = A->colind;
  values = A->values;

  /* Initial Setup */
  p = workspace;
  r = workspace + n;
  q = workspace + 2*n;
  z = workspace + 3*n;
  M = workspace + 4*n;

  for (i=0; i<n; i++) {
    x[i] = 0.0;
    if (values[rowptr[i]] != 0.0)
      M[i] = 1.0/values[rowptr[i]];
    else
      M[i] = 0.0;
  }

  /* r = b - Ax */
  mvMult2(A, x, r);
  for (i=0; i<n; i++)
    r[i] = b[i]-r[i];

  bnrm2 = rnorm2(n, b, 1);
  if (bnrm2 > 0.0) {
    error = rnorm2(n, r, 1) / bnrm2;

    if (error > tol) {
      /* Begin Iterations */
      for (k=0; k<n; k++) {
        for (i=0; i<n; i++)
          z[i] = r[i]*M[i];

        rho = rdot(n, r, 1, z, 1);

        if (k == 0)
          rcopy(n, z, p);
        else {
          if (rho_1 != 0.0)
            beta = rho/rho_1;
          else
            beta = 0.0;
          for (i=0; i<n; i++)
            p[i] = z[i] + beta*p[i];
        }

        mvMult2(A, p, q); /* q = A*p */

        tmp = rdot(n, p, 1, q, 1);
        if (tmp != 0.0)
          alpha = rho/tmp;
        else
          alpha = 0.0;
        raxpy(n,  alpha, p, 1, x, 1);   /* x = x + alpha*p */
        raxpy(n, -alpha, q, 1, r, 1);   /* r = r - alpha*q */
        error = rnorm2(n, r, 1) / bnrm2;
        if (error < tol)
          break;

        rho_1 = rho;
      }
    }
  }
}


/*************************************************************************
* This function performs Matrix-Vector multiplication
**************************************************************************/
void mvMult2(matrix_t *A, real_t *v, real_t *w)
{
  idx_t i, j;

  for (i = 0; i < A->nrows; i++)
    w[i] = 0.0;

  for (i = 0; i < A->nrows; i++)
    for (j = A->rowptr[i]; j < A->rowptr[i+1]; j++)
      w[i] += A->values[j] * v[A->colind[j]];

  return;
  }


/*************************************************************************
* This function sets up the transfer vectors
**************************************************************************/
void ComputeTransferVector(idx_t ncon, matrix_t *matrix, real_t *solution,
  real_t *transfer, idx_t index)
{
  idx_t j, k;
  idx_t nrows;
  idx_t *rowptr, *colind;

  nrows = matrix->nrows;
  rowptr = matrix->rowptr;
  colind = matrix->colind;

  for (j=0; j<nrows; j++) {
    for (k=rowptr[j]+1; k<rowptr[j+1]; k++) {
      if (solution[j] > solution[colind[k]]) {
        transfer[k*ncon+index] = solution[j] - solution[colind[k]];
      }
      else {
        transfer[k*ncon+index] = 0.0;
      }
    }
  }
}

