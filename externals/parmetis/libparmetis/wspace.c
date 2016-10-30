/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * memory.c
 *
 * This file contains routines that deal with memory allocation
 *
 * Started 2/24/96
 * George
 *
 * $Id: wspace.c 10540 2011-07-11 15:42:13Z karypis $
 *
 */

#include <parmetislib.h>


/*************************************************************************/
/*! This function allocate various pools of memory */
/*************************************************************************/
void AllocateWSpace(ctrl_t *ctrl, size_t nwords)
{
  ctrl->mcore = gk_mcoreCreate(nwords*sizeof(idx_t));
}


/*************************************************************************/
/*! This function allocates refinement-specific memory for the workspace */
/*************************************************************************/
void AllocateRefinementWorkSpace(ctrl_t *ctrl, idx_t nbrpoolsize)
{
  ctrl->nbrpoolsize     = nbrpoolsize;
  ctrl->nbrpoolcpos     = 0;
  ctrl->nbrpoolreallocs = 0;

  ctrl->cnbrpool = (cnbr_t *)gk_malloc(ctrl->nbrpoolsize*sizeof(cnbr_t), 
                                  "AllocateRefinementWorkSpace: cnbrpool");

}


/*************************************************************************/
/*! This function de-allocate various pools of memory */
/**************************************************************************/
void FreeWSpace(ctrl_t *ctrl)
{
  ctrl->dbglvl = 0;
  gk_mcoreDestroy(&ctrl->mcore, (ctrl->dbglvl&DBG_INFO));

  if (ctrl->dbglvl&DBG_INFO) {
    printf(" nbrpool statistics [pe:%"PRIDX"]\n" 
           "        nbrpoolsize: %12zu   nbrpoolcpos: %12zu\n"
           "    nbrpoolreallocs: %12zu\n\n",
           ctrl->mype, ctrl->nbrpoolsize,  ctrl->nbrpoolcpos, 
           ctrl->nbrpoolreallocs);
  }

  gk_free((void **)&ctrl->cnbrpool, LTERM);
  ctrl->nbrpoolsize = 0;
  ctrl->nbrpoolcpos = 0;

}


/*************************************************************************/
/*! This function allocate space from the workspace/heap */
/*************************************************************************/
void *wspacemalloc(ctrl_t *ctrl, size_t nbytes)
{
  return gk_mcoreMalloc(ctrl->mcore, nbytes);
}


/*************************************************************************/
/*! This function allocate space from the core  */
/*************************************************************************/
idx_t *iwspacemalloc(ctrl_t *ctrl, size_t n)
{
  return (idx_t *)wspacemalloc(ctrl, n*sizeof(idx_t));
}

/*************************************************************************/
/*! This function resets the cnbrpool */
/*************************************************************************/
void cnbrpoolReset(ctrl_t *ctrl)
{
  ctrl->nbrpoolcpos = 0;
}


/*************************************************************************/
/*! This function gets the next free index from cnbrpool */
/*************************************************************************/
idx_t cnbrpoolGetNext(ctrl_t *ctrl, idx_t nnbrs)
{
  ctrl->nbrpoolcpos += nnbrs;

  if (ctrl->nbrpoolcpos > ctrl->nbrpoolsize) {
    ctrl->nbrpoolsize += gk_max(10*nnbrs, ctrl->nbrpoolsize/2);

    ctrl->cnbrpool = (cnbr_t *)gk_realloc(ctrl->cnbrpool,
                          ctrl->nbrpoolsize*sizeof(cnbr_t), "cnbrpoolGet: cnbrpool");
    ctrl->nbrpoolreallocs++;
  }

  return ctrl->nbrpoolcpos - nnbrs;
}

/*************************************************************************/
/*! This function allocate space from the core */
/*************************************************************************/
real_t *rwspacemalloc(ctrl_t *ctrl, size_t n)
{
  return (real_t *)wspacemalloc(ctrl, n*sizeof(real_t));
}


/*************************************************************************/
/*! This function allocate space from the core */
/*************************************************************************/
ikv_t *ikvwspacemalloc(ctrl_t *ctrl, size_t n)
{
  return (ikv_t *)wspacemalloc(ctrl, n*sizeof(ikv_t));
}


/*************************************************************************/
/*! This function allocate space from the core */
/*************************************************************************/
rkv_t *rkvwspacemalloc(ctrl_t *ctrl, size_t n)
{
  return (rkv_t *)wspacemalloc(ctrl, n*sizeof(rkv_t));
}


