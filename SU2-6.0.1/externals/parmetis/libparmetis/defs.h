/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * defs.h
 *
 * This file contains constant definitions
 *
 * Started 8/27/94
 * George
 *
 * $Id: defs.h 10543 2011-07-11 19:32:24Z karypis $
 *
 */


#define GLOBAL_DBGLVL			0
#define GLOBAL_SEED			15

#define NUM_INIT_MSECTIONS              5

#define MC_FLOW_BALANCE_THRESHOLD       0.2
#define MOC_GD_GRANULARITY_FACTOR       1.0
#define RIP_SPLIT_FACTOR                8
#define MAX_NPARTS_MULTIPLIER		20

#define STATIC_PARTITION        1
#define ORDER_PARTITION		2
#define ADAPTIVE_PARTITION      3
#define REFINE_PARTITION        4

#define REDIST_WGT              2.0
#define MAXNVWGT_FACTOR         2.0

#define N_MOC_REDO_PASSES       10
#define N_MOC_GR_PASSES         8
#define NREMAP_PASSES           8
#define N_MOC_GD_PASSES         6
#define N_MOC_BAL_PASSES        4
#define NMATCH_PASSES           4

#define MAX_NCON_FOR_DIFFUSION  2
#define SMALLGRAPH              10000

#define LTERM                   (void **) 0     /* List terminator for GKfree() */

#define NGD_PASSES		20

#define PMV3_OPTION_DBGLVL	1
#define PMV3_OPTION_SEED	2
#define PMV3_OPTION_IPART	3
#define PMV3_OPTION_PSR		3

#define XYZ_XCOORD		1
#define XYZ_SPFILL		2

/* Type of initial vertex separator algorithms */
#define ISEP_EDGE		1
#define ISEP_NODE		2

#define UNMATCHED		-1
#define MAYBE_MATCHED		-2
#define TOO_HEAVY		-3


#define HTABLE_EMPTY    	-1

#define NGR_PASSES		4	/* Number of greedy refinement passes */
#define NIPARTS			8	/* Number of random initial partitions */
#define NLGR_PASSES		5	/* Number of GR refinement during IPartition */

#define SMALLFLOAT		0.000001


#define COARSEN_FRACTION	0.75	/* Node reduction between succesive coarsening levels */
#define COARSEN_FRACTION2	0.55	/* Node reduction between succesive coarsening levels */
#define UNBALANCE_FRACTION		1.05
#define ORDER_UNBALANCE_FRACTION	1.10

#define MAXVWGT_FACTOR		1.4


/* Debug Levels */
#define DBG_TIME	PARMETIS_DBGLVL_TIME 
#define DBG_INFO	PARMETIS_DBGLVL_INFO
#define DBG_PROGRESS   	PARMETIS_DBGLVL_PROGRESS
#define DBG_REFINEINFO	PARMETIS_DBGLVL_REFINEINFO
#define DBG_MATCHINFO	PARMETIS_DBGLVL_MATCHINFO
#define DBG_RMOVEINFO	PARMETIS_DBGLVL_RMOVEINFO
#define DBG_REMAP	PARMETIS_DBGLVL_REMAP
