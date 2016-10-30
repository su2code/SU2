/*!
 * Copyright 1997, Regents of the University of Minnesota
 *
 * \file
 * \brief Functions dealing with manipulating the ctrl_t structure
 *
 *
 * \date Started 10/19/1996
 * \author George Karypis
 * \version\verbatim $Id: ctrl.c 10592 2011-07-16 21:17:53Z karypis $ \endverbatime
 */

#include <parmetislib.h>



/*************************************************************************/
/*! This function sets the ctrl_t structure 
*/
/*************************************************************************/
ctrl_t *SetupCtrl(pmoptype_et optype, idx_t *options, idx_t ncon, idx_t nparts, 
            real_t *tpwgts, real_t *ubvec, MPI_Comm comm)
{
  idx_t i, j, defopts;
  ctrl_t *ctrl;

  ctrl = (ctrl_t *)gk_malloc(sizeof(ctrl_t), "SetupCtrl: ctrl");
  memset((void *)ctrl, 0, sizeof(ctrl_t));


  /* communicator-related info */
  MPI_Comm_dup(comm, &(ctrl->gcomm));
  ctrl->comm = ctrl->gcomm;
  ctrl->free_comm = 1;
  gkMPI_Comm_rank(ctrl->gcomm, &ctrl->mype);
  gkMPI_Comm_size(ctrl->gcomm, &ctrl->npes);


  /* options[]-related info */
  defopts = (options == NULL ? 1 : options[0] == 0);
  switch (optype) {
    case PARMETIS_OP_KMETIS:
    case PARMETIS_OP_GKMETIS:
      ctrl->partType    = STATIC_PARTITION;
      ctrl->ps_relation = -1;
      break;

    case PARMETIS_OP_GMETIS:
      break;

    case PARMETIS_OP_RMETIS:
      ctrl->partType    = REFINE_PARTITION;
      ctrl->ipc_factor  = 1000.0;
      ctrl->ps_relation = (defopts ? (ctrl->npes == nparts ? PARMETIS_PSR_COUPLED 
                                                           : PARMETIS_PSR_UNCOUPLED)
                                   : (ctrl->npes == nparts ? options[PMV3_OPTION_PSR]
                                                           : PARMETIS_PSR_UNCOUPLED));
      break;

    case PARMETIS_OP_AMETIS:
      ctrl->partType    = ADAPTIVE_PARTITION;
      ctrl->ps_relation = (defopts ? (ctrl->npes == nparts ? PARMETIS_PSR_COUPLED 
                                                           : PARMETIS_PSR_UNCOUPLED)
                                   : (ctrl->npes == nparts ? options[PMV3_OPTION_PSR]
                                                           : PARMETIS_PSR_UNCOUPLED));
      break;

    case PARMETIS_OP_OMETIS:
      /* This is handled directly by the code as its parameter passing does not
         conform to the options[] style. This will probably be changed once the
         changed have been debugged. */
      break;

    case PARMETIS_OP_M2DUAL:
      break;

    case PARMETIS_OP_MKMETIS:
      break;
  }
  ctrl->dbglvl = (defopts ? GLOBAL_DBGLVL : options[PMV3_OPTION_DBGLVL]);
  ctrl->seed   = (defopts ? GLOBAL_SEED   : options[PMV3_OPTION_SEED]);
  ctrl->sync   = GlobalSEMax(ctrl, ctrl->seed);
  ctrl->seed   = (ctrl->seed == 0 ? ctrl->mype : ctrl->seed*ctrl->mype);


  /* common info */
  ctrl->optype        = optype;
  ctrl->ncon          = ncon;    
  ctrl->nparts        = nparts;    
  ctrl->redist_factor = 1.0;
  ctrl->redist_base   = 1.0;

  /* setup tpwgts */
  ctrl->tpwgts = rmalloc(nparts*ncon, "SetupCtrl: tpwgts");
  if (tpwgts) {
    rcopy(nparts*ncon, tpwgts, ctrl->tpwgts);
  }
  else {
    for (i=0; i<nparts; i++) {
      for (j=0; j<ncon; j++)
        ctrl->tpwgts[i*ncon+j] = 1.0/nparts;
    }
  }

  /* setup ubvec */
  ctrl->ubvec = rsmalloc(ncon, UNBALANCE_FRACTION, "SetupCtrl: ubvec");
  if (ubvec)
    rcopy(ncon, ubvec, ctrl->ubvec);

  /* initialize the various timers */
  InitTimers(ctrl);

  /* initialize the random number generator */
  srand(ctrl->seed);

  return ctrl;
}


/*************************************************************************/
/*! This function computes the invtvwgts of a graph and stores them in ctrl
*/
/*************************************************************************/
void SetupCtrl_invtvwgts(ctrl_t *ctrl, graph_t *graph)
{
  idx_t j, ncon;

  ncon  = graph->ncon;

  ctrl->invtvwgts = rmalloc(ncon, "SetupCtrl_tvwgts: invtvwgts");

  for (j=0; j<ncon; j++) 
    ctrl->invtvwgts[j] = 1.0/GlobalSESum(ctrl, isum(graph->nvtxs, graph->vwgt+j, ncon));
    
}


/*************************************************************************/
/*! This function de-allocates memory allocated for the control structures 
*/
/*************************************************************************/
void FreeCtrl(ctrl_t **r_ctrl)
{
  ctrl_t *ctrl = *r_ctrl;

  FreeWSpace(ctrl);

  if (ctrl->free_comm)
    gkMPI_Comm_free(&(ctrl->gcomm));

  gk_free((void **)&ctrl->invtvwgts, 
      &ctrl->ubvec, &ctrl->tpwgts, 
      &ctrl->sreq, &ctrl->rreq, &ctrl->statuses,
      &ctrl, 
      LTERM);

  *r_ctrl = NULL;
}

