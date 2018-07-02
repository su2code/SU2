/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * Started 9/28/95
 * George
 *
 * $Id: util.c 10057 2011-06-02 13:44:44Z karypis $
 */

#include <parmetislib.h>


/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void myprintf(ctrl_t *ctrl, char *f_str,...)
{
  va_list argp;

  fprintf(stdout, "[%2"PRIDX"] ", ctrl->mype);

  va_start(argp, f_str);
  vfprintf(stdout, f_str, argp);
  va_end(argp);

  if (strlen(f_str) == 0 || f_str[strlen(f_str)-1] != '\n')
    fprintf(stdout,"\n");
  fflush(stdout);

}


/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void rprintf(ctrl_t *ctrl, char *f_str,...)
{
  va_list argp;

  if (ctrl->mype == 0) {
    va_start(argp, f_str);
    vfprintf(stdout, f_str, argp);
    va_end(argp);
  }

  fflush(stdout);

  gkMPI_Barrier(ctrl->comm);

}


/*************************************************************************
* This function does a binary search on an array for a key and returns
* the index
**************************************************************************/
idx_t BSearch(idx_t n, idx_t *array, idx_t key)
{
  idx_t a=0, b=n, c;

  while (b-a > 8) {
    c = (a+b)>>1;
    if (array[c] > key)
      b = c;
    else
      a = c;
  }

  for (c=a; c<b; c++) {
    if (array[c] == key)
      return c;
  }

  errexit("Key %"PRIDX" not found!\n", key);

  return 0;
}


/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i 
**************************************************************************/
void RandomPermute(idx_t n, idx_t *p, idx_t flag)
{
  idx_t i, u, v;
  idx_t tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  for (i=0; i<n; i++) {
    v = RandomInRange(n);
    u = RandomInRange(n);
   gk_SWAP(p[v], p[u], tmp);
  }
}


/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i 
**************************************************************************/
void FastRandomPermute(idx_t n, idx_t *p, idx_t flag)
{
  idx_t i, u, v;
  idx_t tmp;

  /* this is for very small arrays */
  if (n < 25) {
    RandomPermute(n, p, flag);
    return;
  }

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  for (i=0; i<n; i+=8) {
    v = RandomInRange(n-4);
    u = RandomInRange(n-4);
   gk_SWAP(p[v], p[u], tmp);
   gk_SWAP(p[v+1], p[u+1], tmp);
   gk_SWAP(p[v+2], p[u+2], tmp);
   gk_SWAP(p[v+3], p[u+3], tmp);
  }
}

/*************************************************************************
* This function returns true if the a is a power of 2
**************************************************************************/
idx_t ispow2(idx_t a)
{
  for (; a%2 != 1; a = a>>1);
  return (a > 1 ? 0 : 1);
}

/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
idx_t log2Int(idx_t a)
{
  idx_t i;

  for (i=1; a > 1; i++, a = a>>1);
  return i-1;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
size_t rargmax_strd(size_t n, real_t *x, size_t incx)
{
  size_t i, max=0;

  n *= incx;
  for (i=incx; i<n; i+=incx)
    max = (x[i] > x[max] ? i : max);

  return max/incx;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
size_t rargmin_strd(size_t n, real_t *x, size_t incx)
{
  size_t i, min=0;

  n *= incx;
  for (i=incx; i<n; i+=incx)
    min = (x[i] < x[min] ? i : min);

  return min/incx;
}


/*************************************************************************
* These functions return the index of the almost maximum element in a vector
**************************************************************************/
size_t rargmax2(size_t n, real_t *x)
{
  size_t i, max1, max2;

  if (x[0] > x[1]) {
    max1 = 0;
    max2 = 1;
  }
  else {
    max1 = 1;
    max2 = 0;
  }

  for (i=2; i<n; i++) {
    if (x[i] > x[max1]) {
      max2 = max1;
      max1 = i;
    }
    else if (x[i] > x[max2])
      max2 = i;
  }

  return max2;
}


/*************************************************************************
* This function returns the average value of an array
**************************************************************************/
real_t ravg(size_t n, real_t *x)
{
  size_t i;
  real_t retval = 0.0;

  for (i=0; i<n; i++)
    retval += x[i];

  return retval/n;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
real_t rfavg(size_t n, real_t *x)
{
  size_t i;
  real_t total = 0.0;

  if (n == 0)
    return 0.0;

  for (i=0; i<n; i++)
    total += fabs(x[i]);

  return total/n;
}


/*************************************************************************
* This function checks if v+u2 provides a better balance in the weight
* vector that v+u1
**************************************************************************/
real_t BetterVBalance(idx_t ncon, real_t *vwgt, real_t *u1wgt, real_t *u2wgt)
{
  idx_t i;
  real_t sum1, sum2, diff1, diff2;

  if (ncon == 1)
    return u1wgt[0] - u1wgt[0];

  sum1 = sum2 = 0.0;
  for (i=0; i<ncon; i++) {
    sum1 += vwgt[i]+u1wgt[i];
    sum2 += vwgt[i]+u2wgt[i];
  }
  sum1 = sum1/(1.0*ncon);
  sum2 = sum2/(1.0*ncon);

  diff1 = diff2 = 0.0;
  for (i=0; i<ncon; i++) {
    diff1 += fabs(sum1 - (vwgt[i]+u1wgt[i]));
    diff2 += fabs(sum2 - (vwgt[i]+u2wgt[i]));
  }

  return diff1 - diff2;

}


/*************************************************************************
* This function checks if the pairwise balance of the between the two
* partitions will improve by moving the vertex v from pfrom to pto,
* subject to the target partition weights of tfrom, and tto respectively
**************************************************************************/
idx_t IsHBalanceBetterFT(idx_t ncon, real_t *pfrom, real_t *pto, real_t *nvwgt, real_t *ubvec)
{
  idx_t i;
  real_t blb1=0.0, alb1=0.0, sblb=0.0, salb=0.0;
  real_t blb2=0.0, alb2=0.0;
  real_t temp;

  for (i=0; i<ncon; i++) {
    temp =gk_max(pfrom[i], pto[i])/ubvec[i];
    if (blb1 < temp) {
      blb2 = blb1;
      blb1 = temp;
    }
    else if (blb2 < temp)
      blb2 = temp;
    sblb += temp;

    temp =gk_max(pfrom[i]-nvwgt[i], pto[i]+nvwgt[i])/ubvec[i];
    if (alb1 < temp) {
      alb2 = alb1;
      alb1 = temp;
    }
    else if (alb2 < temp)
      alb2 = temp;
    salb += temp;
  }

  if (alb1 < blb1)
    return 1;
  if (blb1 < alb1)
    return 0;
  if (alb2 < blb2)
    return 1;
  if (blb2 < alb2)
    return 0;

  return salb < sblb;

}


/*************************************************************************
* This function checks if it will be better to move a vertex to pt2 than
* to pt1 subject to their target weights of tt1 and tt2, respectively
* This routine takes into account the weight of the vertex in question
**************************************************************************/
idx_t IsHBalanceBetterTT(idx_t ncon, real_t *pt1, real_t *pt2, real_t *nvwgt, real_t *ubvec)
{
  idx_t i;
  real_t m11=0.0, m12=0.0, m21=0.0, m22=0.0, sm1=0.0, sm2=0.0, temp;

  for (i=0; i<ncon; i++) {
    temp = (pt1[i]+nvwgt[i])/ubvec[i];
    if (m11 < temp) {
      m12 = m11;
      m11 = temp;
    }
    else if (m12 < temp)
      m12 = temp;
    sm1 += temp;
    temp = (pt2[i]+nvwgt[i])/ubvec[i];
    if (m21 < temp) {
      m22 = m21;
      m21 = temp;
    }
    else if (m22 < temp)
      m22 = temp;
    sm2 += temp;
  }
  if (m21 < m11)
    return 1;
  if (m21 > m11)
    return 0;
  if (m22 < m12)
    return 1;
  if (m22 > m12)
    return 0;

  return sm2 < sm1;
}


/*************************************************************************
* This function computes the top three values of a real_t array
**************************************************************************/
void GetThreeMax(idx_t n, real_t *x, idx_t *first, idx_t *second, idx_t *third)
{
  idx_t i;

  if (n <= 0) {
    *first = *second = *third = -1;
    return;
  }

  *second = *third = -1;
  *first = 0;

  for (i=1; i<n; i++) {
    if (x[i] > x[*first]) {
      *third = *second;
      *second = *first;
      *first = i;
      continue;
    }

    if (*second == -1 || x[i] > x[*second]) {
      *third = *second;
      *second = i;
      continue;
    }

    if (*third == -1 || x[i] > x[*third])
      *third = i;
  }

  return;
}
