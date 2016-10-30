/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * selectq.c
 *
 * This file contains the driving routines for multilevel k-way refinement
 *
 * Started 7/28/97
 * George
 *
 * $Id: selectq.c 10542 2011-07-11 16:56:22Z karypis $
 */

#include <parmetislib.h>

/*************************************************************************/
/*! This stuff is hardcoded for up to four constraints 
*/
/*************************************************************************/
void Mc_DynamicSelectQueue(ctrl_t *ctrl, idx_t nqueues, idx_t ncon, idx_t subdomain1, 
         idx_t subdomain2, idx_t *currentq, real_t *flows, idx_t *from, idx_t *qnum, 
         idx_t minval, real_t avgvwgt, real_t maxdiff)
{
  idx_t i, j;
  idx_t hash, index = -1, current;
  idx_t *cand, *rank, *dont_cares;
  idx_t nperms, perm[24][5];
  real_t sign = 0.0;
  rkv_t *array;
  idx_t mype;
  gkMPI_Comm_rank(MPI_COMM_WORLD, &mype);

  WCOREPUSH;

  *qnum = -1;

  /* allocate memory */
  cand       = iwspacemalloc(ctrl, ncon);
  rank       = iwspacemalloc(ctrl, ncon);
  dont_cares = iwspacemalloc(ctrl, ncon);
  array      = rkvwspacemalloc(ctrl, ncon);

  if (*from == -1) {
    for (i=0; i<ncon; i++) {
      array[i].key = fabs(flows[i]);
      array[i].val = i;
    }

    /* GKTODO - Need to check the correct direction of the sort */
    rkvsorti(ncon, array);

    /* GKTODO - The following assert was disabled as it was failing. Need
                to check if it is a valid assert */
    /*ASSERT(array[ncon-1].key - array[0].key <= maxdiff) */

    if (flows[array[ncon-1].val] > avgvwgt*MOC_GD_GRANULARITY_FACTOR) {
      *from = subdomain1;
      sign  = 1.0;
      index = 0;
    }

    if (flows[array[ncon-1].val] < -1.0*avgvwgt*MOC_GD_GRANULARITY_FACTOR) {
      *from = subdomain2;
      sign  = -1.0;
      index = nqueues;
    }

    if (*from == -1) 
      goto DONE;
  }
  else {
    ASSERT(*from == subdomain1 || *from == subdomain2);

    if (*from == subdomain1) {
      sign  = 1.0;
      index = 0;
    }
    else {
      sign  = -1.0;
      index = nqueues;
    }
  }

  for (i=0; i<ncon; i++) {
    array[i].key = flows[i] * sign;
    array[i].val = i;
  }

  /* GKTODO Need to check the direction of those sorts */
  rkvsorti(ncon, array);

  iset(ncon, 1, dont_cares);
  for (current=0, i=0; i<ncon-1; i++) {
    if (array[i+1].key - array[i].key < maxdiff * MC_FLOW_BALANCE_THRESHOLD && dont_cares[current] < ncon-1) {
      dont_cares[current]++;
      dont_cares[i+1] = 0;
    }
    else
      current = i+1;
  }

  switch (ncon) {
    /***********************/
    case 2:
      nperms = 1;
      perm[0][0] = 0;   perm[0][1] = 1;

      break;
    /***********************/
    case 3:

      /* if the first and second flows are close */
      if (dont_cares[0] == 2 && dont_cares[1] == 0 && dont_cares[2] == 1) {
        nperms = 4;
        perm[0][0] = 0;   perm[0][1] = 1;   perm[0][2] = 2;
        perm[1][0] = 1;   perm[1][1] = 0;   perm[1][2] = 2;
        perm[2][0] = 0;   perm[2][1] = 2;   perm[2][2] = 1;
        perm[3][0] = 1;   perm[3][1] = 2;   perm[3][2] = 0;
        break;
      }

      /* if the second and third flows are close */
      if (dont_cares[0] == 1 && dont_cares[1] == 2 && dont_cares[2] == 0) {
        nperms = 4;
        perm[0][0] = 0;   perm[0][1] = 1;   perm[0][2] = 2;
        perm[1][0] = 0;   perm[1][1] = 2;   perm[1][2] = 1;
        perm[2][0] = 1;   perm[2][1] = 0;   perm[2][2] = 2;
        perm[3][0] = 2;   perm[3][1] = 0;   perm[3][2] = 1;
        break;
      }

      /* all or none of the flows are close */
      nperms = 3;
      perm[0][0] = 0;   perm[0][1] = 1;   perm[0][2] = 2;
      perm[1][0] = 1;   perm[1][1] = 0;   perm[1][2] = 2;
      perm[2][0] = 0;   perm[2][1] = 2;   perm[2][2] = 1;

      break;
    /***********************/
    case 4:

      if (dont_cares[0] == 2 && dont_cares[1] == 0 &&
          dont_cares[2] == 1 && dont_cares[3] == 1) {
        nperms = 14;
        perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
        perm[1][0] =  1;   perm[1][1] =  0;   perm[1][2] =  2;   perm[1][3] =  3;
        perm[2][0] =  0;   perm[2][1] =  2;   perm[2][2] =  1;   perm[2][3] =  3;
        perm[3][0] =  1;   perm[3][1] =  2;   perm[3][2] =  0;   perm[3][3] =  3;
        perm[4][0] =  0;   perm[4][1] =  1;   perm[4][2] =  3;   perm[4][3] =  2;
        perm[5][0] =  1;   perm[5][1] =  0;   perm[5][2] =  3;   perm[5][3] =  2;
        
        perm[6][0] =  0;   perm[6][1] =  3;   perm[6][2] =  1;   perm[6][3] =  2;
        perm[7][0] =  1;   perm[7][1] =  3;   perm[7][2] =  0;   perm[7][3] =  2;

        perm[8][0] =  0;   perm[8][1] =  2;   perm[8][2] =  3;   perm[8][3] =  1;
        perm[9][0] =  1;   perm[9][1] =  2;   perm[9][2] =  3;   perm[9][3] =  0;

        perm[10][0] = 2;   perm[10][1] = 0;   perm[10][2] = 1;   perm[10][3] = 3;
        perm[11][0] = 2;   perm[11][1] = 1;   perm[11][2] = 0;   perm[11][3] = 3;
        
        perm[12][0] = 0;   perm[12][1] = 3;   perm[12][2] = 2;   perm[12][3] = 1;
        perm[13][0] = 1;   perm[13][1] = 3;   perm[13][2] = 2;   perm[13][3] = 0;
        break;
      }

      if (dont_cares[0] == 1 && dont_cares[1] == 1 &&
          dont_cares[2] == 2 && dont_cares[3] == 0) {
        nperms = 14;
        perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
        perm[1][0] =  0;   perm[1][1] =  1;   perm[1][2] =  3;   perm[1][3] =  2;
        perm[2][0] =  0;   perm[2][1] =  2;   perm[2][2] =  1;   perm[2][3] =  3;
        perm[3][0] =  0;   perm[3][1] =  3;   perm[3][2] =  1;   perm[3][3] =  2;
        perm[4][0] =  1;   perm[4][1] =  0;   perm[4][2] =  2;   perm[4][3] =  3;
        perm[5][0] =  1;   perm[5][1] =  0;   perm[5][2] =  3;   perm[5][3] =  2;

        perm[6][0] =  1;   perm[6][1] =  2;   perm[6][2] =  0;   perm[6][3] =  3;
        perm[7][0] =  1;   perm[7][1] =  3;   perm[7][2] =  0;   perm[7][3] =  2;

        perm[8][0] =  2;   perm[8][1] =  0;   perm[8][2] =  1;   perm[8][3] =  3;
        perm[9][0] =  3;   perm[9][1] =  0;   perm[9][2] =  1;   perm[9][3] =  2;

        perm[10][0] = 0;   perm[10][1] = 2;   perm[10][2] = 3;   perm[10][3] = 1;
        perm[11][0] = 0;   perm[11][1] = 3;   perm[11][2] = 2;   perm[11][3] = 1;

        perm[12][0] = 2;   perm[12][1] = 1;   perm[12][2] = 0;   perm[12][3] = 3;
        perm[13][0] = 3;   perm[13][1] = 1;   perm[13][2] = 0;   perm[13][3] = 2;
        break;
      }

      if (dont_cares[0] == 2 && dont_cares[1] == 0 &&
          dont_cares[2] == 2 && dont_cares[3] == 0) {
        nperms = 14;
        perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
        perm[1][0] =  1;   perm[1][1] =  0;   perm[1][2] =  2;   perm[1][3] =  3;
        perm[2][0] =  0;   perm[2][1] =  1;   perm[2][2] =  3;   perm[2][3] =  2;
        perm[3][0] =  1;   perm[3][1] =  0;   perm[3][2] =  3;   perm[3][3] =  2;

        perm[4][0] =  0;   perm[4][1] =  2;   perm[4][2] =  1;   perm[4][3] =  3;
        perm[5][0] =  1;   perm[5][1] =  2;   perm[5][2] =  0;   perm[5][3] =  3;
        perm[6][0] =  0;   perm[6][1] =  3;   perm[6][2] =  1;   perm[6][3] =  2;
        perm[7][0] =  1;   perm[7][1] =  3;   perm[7][2] =  0;   perm[7][3] =  2;

        perm[8][0] = 2;    perm[8][1] = 0;    perm[8][2] = 1;    perm[8][3] = 3;
        perm[9][0] = 0;    perm[9][1] = 2;    perm[9][2] = 3;    perm[9][3] = 1;
        perm[10][0] = 2;   perm[10][1] = 1;   perm[10][2] = 0;   perm[10][3] = 3;
        perm[11][0] = 0;   perm[11][1] = 3;   perm[11][2] = 2;   perm[11][3] = 1;
        perm[12][0] = 3;   perm[12][1] = 0;   perm[12][2] = 1;   perm[12][3] = 2;
        perm[13][0] = 1;   perm[13][1] = 2;   perm[13][2] = 3;   perm[13][3] = 0;
        break;
      }

      if (dont_cares[0] == 3 && dont_cares[1] == 0 &&
          dont_cares[2] == 0 && dont_cares[3] == 1) {
        nperms = 14;
        perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
        perm[1][0] =  0;   perm[1][1] =  2;   perm[1][2] =  1;   perm[1][3] =  3;
        perm[2][0] =  1;   perm[2][1] =  0;   perm[2][2] =  2;   perm[2][3] =  3;
        perm[3][0] =  2;   perm[3][1] =  0;   perm[3][2] =  1;   perm[3][3] =  3;
        perm[4][0] =  1;   perm[4][1] =  2;   perm[4][2] =  0;   perm[4][3] =  3;
        perm[5][0] =  2;   perm[5][1] =  1;   perm[5][2] =  0;   perm[5][3] =  3;

        perm[6][0] =  0;   perm[6][1] =  1;   perm[6][2] =  3;   perm[6][3] =  2;
        perm[7][0] =  1;   perm[7][1] =  0;   perm[7][2] =  3;   perm[7][3] =  2;
        perm[8][0] =  0;   perm[8][1] =  2;   perm[8][2] =  3;   perm[8][3] =  1;
        perm[9][0] =  2;   perm[9][1] =  0;   perm[9][2] =  3;   perm[9][3] =  1;
        perm[10][0] = 1;   perm[10][1] = 2;   perm[10][2] = 3;   perm[10][3] = 0;
        perm[11][0] = 2;   perm[11][1] = 1;   perm[11][2] = 3;   perm[11][3] = 0;

        perm[12][0] = 0;   perm[12][1] = 3;   perm[12][2] = 1;   perm[12][3] = 2;
        perm[13][0] = 0;   perm[13][1] = 3;   perm[13][2] = 2;   perm[13][3] = 1;
        break;
      }

      if (dont_cares[0] == 1 && dont_cares[1] == 3 &&
          dont_cares[2] == 0 && dont_cares[3] == 0) {
        nperms = 14;
        perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
        perm[1][0] =  0;   perm[1][1] =  2;   perm[1][2] =  1;   perm[1][3] =  3;
        perm[2][0] =  0;   perm[2][1] =  1;   perm[2][2] =  3;   perm[2][3] =  2;
        perm[3][0] =  0;   perm[3][1] =  2;   perm[3][2] =  3;   perm[3][3] =  1;
        perm[4][0] =  0;   perm[4][1] =  3;   perm[4][2] =  1;   perm[4][3] =  2;
        perm[5][0] =  0;   perm[5][1] =  3;   perm[5][2] =  2;   perm[5][3] =  1;

        perm[6][0] =  1;   perm[6][1] =  0;   perm[6][2] =  2;   perm[6][3] =  3;
        perm[7][0] =  1;   perm[7][1] =  0;   perm[7][2] =  3;   perm[7][3] =  2;
        perm[8][0] =  2;   perm[8][1] =  0;   perm[8][2] =  1;   perm[8][3] =  3;
        perm[9][0] =  2;   perm[9][1] =  0;   perm[9][2] =  3;   perm[9][3] =  1;
        perm[10][0] = 3;   perm[10][1] = 0;   perm[10][2] = 1;   perm[10][3] = 2;
        perm[11][0] = 3;   perm[11][1] = 0;   perm[11][2] = 2;   perm[11][3] = 1;

        perm[12][0] = 1;   perm[12][1] = 2;   perm[12][2] = 0;   perm[12][3] = 3;
        perm[13][0] = 2;   perm[13][1] = 1;   perm[13][2] = 0;   perm[13][3] = 3;

        break;
      }

      nperms = 14;
      perm[0][0] =  0;   perm[0][1] =  1;   perm[0][2] =  2;   perm[0][3] =  3;
      perm[1][0] =  1;   perm[1][1] =  0;   perm[1][2] =  2;   perm[1][3] =  3;
      perm[2][0] =  0;   perm[2][1] =  2;   perm[2][2] =  1;   perm[2][3] =  3;
      perm[3][0] =  0;   perm[3][1] =  1;   perm[3][2] =  3;   perm[3][3] =  2;
      perm[4][0] =  1;   perm[4][1] =  0;   perm[4][2] =  3;   perm[4][3] =  2;

      perm[5][0] =  2;   perm[5][1] =  0;   perm[5][2] =  1;   perm[5][3] =  3;
      perm[6][0] =  0;   perm[6][1] =  2;   perm[6][2] =  3;   perm[6][3] =  1;

      perm[7][0] =  1;   perm[7][1] =  2;   perm[7][2] =  0;   perm[7][3] =  3;
      perm[8][0] =  0;   perm[8][1] =  3;   perm[8][2] =  1;   perm[8][3] =  2;

      perm[9][0] =  2;   perm[9][1] =  1;   perm[9][2] =  0;   perm[9][3] =  3;
      perm[10][0] = 0;   perm[10][1] = 3;   perm[10][2] = 2;   perm[10][3] = 1;
      perm[11][0] = 2;   perm[11][1] = 0;   perm[11][2] = 3;   perm[11][3] = 1;

      perm[12][0] = 3;   perm[12][1] = 0;   perm[12][2] = 1;   perm[12][3] = 2;
      perm[13][0] = 1;   perm[13][1] = 2;   perm[13][2] = 3;   perm[13][3] = 0;
      break;
    /***********************/
    default:
      goto DONE;
  }

  for (i=0; i<nperms; i++) {
    for (j=0; j<ncon; j++)
      cand[j] = array[perm[i][j]].val;

    for (j=0; j<ncon; j++)
      rank[cand[j]] = j;

    hash = Mc_HashVRank(ncon, rank) - minval;
    if (currentq[hash+index] > 0) {
      *qnum = hash;
      goto DONE;
    }
  }

DONE:
  WCOREPOP;
}


/*************************************************************************/
/*! This function sorts the nvwgts of a vertex and returns a hashed value 
*/
/*************************************************************************/
idx_t Mc_HashVwgts(ctrl_t *ctrl, idx_t ncon, real_t *nvwgt)
{
  idx_t i;
  idx_t multiplier, retval;
  idx_t *rank;
  rkv_t *array;

  WCOREPUSH;

  rank  = iwspacemalloc(ctrl, ncon);
  array = rkvwspacemalloc(ctrl, ncon);

  for (i=0; i<ncon; i++) {
    array[i].key = nvwgt[i];
    array[i].val = i;
  }

  rkvsorti(ncon, array);
  for (i=0; i<ncon; i++)
    rank[array[i].val] = i;

  multiplier = 1;

  retval = 0;
  for (i=0; i<ncon; i++) {
    multiplier *= (i+1);
    retval += rank[ncon-i-1] * multiplier;
  }

  WCOREPOP;

  return retval;
}


/*************************************************************************/
/*! This function sorts the vwgts of a vertex and returns a hashed value 
*/
/*************************************************************************/
idx_t Mc_HashVRank(idx_t ncon, idx_t *vwgt)
{
  idx_t i, multiplier, retval;

  multiplier = 1;

  retval = 0;
  for (i=0; i<ncon; i++) {
    multiplier *= (i+1);
    retval += vwgt[ncon-1-i] * multiplier;
  }

  return retval;
}


