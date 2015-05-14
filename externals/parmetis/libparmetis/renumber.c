/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mgrsetup.c
 *
 * This file contain various graph setting up routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: renumber.c 10531 2011-07-09 21:58:13Z karypis $
 *
 */

#include <parmetislib.h>




/*************************************************************************
* This function changes the numbering from 1 to 0 or 0 to 1
**************************************************************************/
void ChangeNumbering(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *part, idx_t npes, idx_t mype, idx_t from)
{
  idx_t i, nvtxs;

  nvtxs = vtxdist[mype+1]-vtxdist[mype];

  if (from == 1) {  /* Change it from 1 to 0 */
    for (i=0; i<npes+1; i++)
      vtxdist[i]--;

    for (i=0; i<nvtxs+1; i++) 
      xadj[i]--;
    for (i=0; i<xadj[nvtxs]; i++) 
      adjncy[i]--;
  }
  else {  /* Change it from 0 to 1 */
    for (i=0; i<npes+1; i++) 
      vtxdist[i]++;

    for (i=0; i<xadj[nvtxs]; i++) 
      adjncy[i]++; 
    for (i=0; i<nvtxs+1; i++) 
      xadj[i]++; 

    for (i=0; i<nvtxs; i++)
      part[i]++;

  }
}


/*************************************************************************
* This function changes the numbering from 1 to 0 or 0 to 1
**************************************************************************/
void ChangeNumberingMesh(idx_t *elmdist, idx_t *eptr, idx_t *eind, 
                         idx_t *xadj, idx_t *adjncy, idx_t *part, 
			 idx_t npes, idx_t mype, idx_t from)
{
  idx_t i, nelms;

  nelms = elmdist[mype+1]-elmdist[mype];

  if (from == 1) {  /* Change it from 1 to 0 */
    for (i=0; i<npes+1; i++)
      elmdist[i]--;

    for (i=0; i<nelms+1; i++) 
      eptr[i]--;
    for (i=0; i<eptr[nelms]; i++) 
      eind[i]--;
  }
  else {  /* Change it from 0 to 1 */
    for (i=0; i<npes+1; i++) 
      elmdist[i]++;

    for (i=0; i<eptr[nelms]; i++) 
      eind[i]++;
    for (i=0; i<nelms+1; i++) 
      eptr[i]++;

    for (i=0; i<xadj[nelms]; i++) 
      adjncy[i]++; 
    for (i=0; i<nelms+1; i++) 
      xadj[i]++; 

    if (part != NULL)
      for (i=0; i<nelms; i++)
        part[i]++;
  }
}


