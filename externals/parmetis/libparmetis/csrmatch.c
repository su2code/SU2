/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * csrmatch.c
 *
 * This file contains the code that computes matchings
 *
 * Started 7/23/97
 * George
 *
 * $Id: csrmatch.c 10057 2011-06-02 13:44:44Z karypis $
 *
 */

#include <parmetislib.h>


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void CSR_Match_SHEM(matrix_t *matrix, idx_t *match, idx_t *mlist,
          idx_t *skip, idx_t ncon)
{
  idx_t h, i, ii, j;
  idx_t nrows, edge, maxidx, count;
  real_t maxwgt;
  idx_t *rowptr, *colind;
  real_t *transfer;
  rkv_t *links;

  nrows    = matrix->nrows;
  rowptr   = matrix->rowptr;
  colind   = matrix->colind;
  transfer = matrix->transfer;

  iset(nrows, UNMATCHED, match);

  links = rkvmalloc(nrows, "links");
  for (i=0; i<nrows; i++) {
    links[i].key = 0.0;
    links[i].val = i;
    for (j=rowptr[i]; j<rowptr[i+1]; j++) {
      for (h=0; h<ncon; h++) {
        if (links[i].key < fabs(transfer[j*ncon+h]))
          links[i].key = fabs(transfer[j*ncon+h]);
      }
    }
  }

  rkvsortd(nrows, links);

  for (count=0, ii=0; ii<nrows; ii++) {
    i = links[ii].val;

    if (match[i] == UNMATCHED) {
      maxidx = i;
      maxwgt = 0.0;

      /* Find a heavy-edge matching */
      for (j=rowptr[i]; j<rowptr[i+1]; j++) {
        edge = colind[j];
        if (match[edge] == UNMATCHED && edge != i && skip[j] == 0) {
          for (h=0; h<ncon; h++)
            if (maxwgt < fabs(transfer[j*ncon+h]))
              break;

          if (h != ncon) {
            maxwgt = fabs(transfer[j*ncon+h]);
            maxidx = edge;
          }
        }
      }

      if (maxidx != i) {
        match[i] = maxidx;
        match[maxidx] = i;
        mlist[count++] = gk_max(i, maxidx);
        mlist[count++] = gk_min(i, maxidx);
      }
    }
  }

  gk_free((void **)&links, LTERM);
}

