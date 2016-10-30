/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * struct.h
 *
 * This file contains data structures for ILU routines.
 *
 * Started 9/26/95
 * George
 *
 * $Id: struct.h 10592 2011-07-16 21:17:53Z karypis $
 */


/*************************************************************************/
/*! This data structure stores cut-based k-way refinement info about an
 *     adjacent subdomain for a given vertex. */
/*************************************************************************/
typedef struct cnbr_t {
  idx_t pid;            /*!< The partition ID */
  idx_t ed;             /*!< The sum of the weights of the adjacent edges
                             that are incident on pid */
} cnbr_t;


/*************************************************************************
* The following data structure stores key-key-value triplets
**************************************************************************/
typedef struct i2kv_t {
  idx_t key1, key2;
  idx_t val;
} i2kv_t;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
typedef struct ckrinfo_t {
 idx_t id;              /*!< The internal degree of a vertex (sum of weights) */
 idx_t ed;              /*!< The total external degree of a vertex */
 idx_t nnbrs;           /*!< The number of neighboring subdomains */
 idx_t inbr;            /*!< The index in the cnbr_t array where the nnbrs list 
                             of neighbors is stored */
} ckrinfo_t;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct nrinfodef {
 idx_t edegrees[2];  
};

typedef struct nrinfodef NRInfoType;


/*************************************************************************
* The following data structure stores a sparse matrix in CSR format
* The diagonal entry is in the first position of each row.
**************************************************************************/
typedef struct matrix_t {
  idx_t nrows, nnzs;		/* Number of rows and nonzeros in the matrix */
  idx_t *rowptr;
  idx_t *colind;
  real_t *values;
  real_t *transfer;
} matrix_t;


/*************************************************************************
* This data structure holds the input graph
**************************************************************************/
typedef struct graph_t {
  idx_t gnvtxs, nvtxs, nedges, ncon, nobj;
  idx_t *xadj;		/* Pointers to the locally stored vertices */
  idx_t *vwgt;		/* Vertex weights */
  real_t *nvwgt;        /* Vertex weights */
  idx_t *vsize;		/* Vertex size */
  idx_t *adjncy;	/* Array that stores the adjacency lists of nvtxs */
  idx_t *adjwgt;	/* Array that stores the weights of the adjacency lists */
  idx_t *vtxdist;	/* Distribution of vertices */
  idx_t *home;		/* The initial partition of the vertex */

  /* used for not freeing application supplied arrays */
  idx_t free_vwgt;
  idx_t free_adjwgt;
  idx_t free_vsize;

  /* Coarsening structures */
  idx_t *match;
  idx_t *cmap;

  /* Used during initial partitioning */
  idx_t *label;

  /* Communication/Setup parameters */
  idx_t nnbrs;                  /*!< The number of neighboring processors */
  idx_t nrecv;                  /*!< The total number of remote vertices that need to 
                                     be received. nrecv == recvptr[nnbrs] */
  idx_t nsend;                  /*!< The total number of local vertices that need to 
                                     be sent. This corresponds to the communication 
                                     volume of each pe, in the sense that if a vertex 
                                     needs to be sent multiple times, it is accounted 
                                     in nsend. nsend == sendptr[nnbrs] */
  idx_t *peind;	                /*!< Array of size nnbrs storing the neighboring PEs */
  idx_t *sendptr, *sendind;     /*!< CSR format of the vertices that are sent to each
                                     of the neighboring processors */
  idx_t *recvptr, *recvind;     /*!< CSR format of the vertices that are received from
                                     each of the neighboring PEs. */
  idx_t *imap;			/*!< The inverse map of local to global indices */
  idx_t *pexadj, *peadjncy, 
        *peadjloc;	        /*!< CSR format of the PEs each vertex is adjancent to 
                                     along with the location in the sendind of the 
                                     non-local adjancent vertices */

  idx_t nlocal;			/*!< Number of interior vertices */
  idx_t *lperm;		        /*!< lperm[0:nlocal] points to interior vertices, 
                                     the rest are interface */

  /* Communication parameters for projecting the partition. 
   * These are computed during CreateCoarseGraph and used during projection 
   * Note that during projection, the meaning of received and sent is reversed! */
  idx_t *rlens, *slens;	/* Arrays of size nnbrs of how many vertices you are sending and receiving */
  ikv_t *rcand;


  /* Partition parameters */
  idx_t *where;
  idx_t *lpwgts, *gpwgts;
  real_t *lnpwgts, *gnpwgts;
  ckrinfo_t *ckrinfo;

  /* Node refinement information */
  idx_t nsep;  			/* The number of vertices in the separator */
  NRInfoType *nrinfo;
  idx_t *sepind;		/* The indices of the vertices in the separator */

  idx_t lmincut, mincut;

  idx_t level;
  idx_t match_type;
  idx_t edgewgt_type;

  struct graph_t *coarser, *finer;
} graph_t;



/*************************************************************************
* The following data type implements a timer
**************************************************************************/
typedef double timer;


/*************************************************************************
* The following structure stores information used by parallel kmetis
**************************************************************************/
typedef struct ctrl_t {
  pmoptype_et optype;           /*!< The operation being performed */
  idx_t mype, npes;		/* Info about the parallel system */
  idx_t ncon;                   /*!< The number of balancing constraints */ 
  idx_t CoarsenTo;		/* The # of vertices in the coarsest graph */
  idx_t dbglvl;			/* Controls the debuging output of the program */
  idx_t nparts;			/* The number of partitions */
  idx_t foldf;			/* What is the folding factor */
  idx_t mtype;                  /* The matching type */
  idx_t ipart;			/* The initial partitioning type */
  idx_t rtype;                  /* The refinement type */
  idx_t p_nseps;                /* The number of separators to compute at each 
                                   parallel bisection */
  idx_t s_nseps;                /* The number of separators to compute at each 
                                   serial bisection */
  real_t ubfrac;                /* The max/avg fraction for separator bisections */
  idx_t seed;			/* Random number seed */
  idx_t sync;			/* Random number seed */
  real_t *tpwgts;		/* Target subdomain weights */
  real_t *invtvwgts;            /* Per-constraint 1/total vertex weight */
  real_t *ubvec;                /* Per-constraint unbalance factor */

  idx_t partType;
  idx_t ps_relation;

  real_t redist_factor;
  real_t redist_base;
  real_t ipc_factor;
  real_t edge_size_ratio;
  matrix_t *matrix;

  idx_t free_comm;       /*!< Used to indicate if gcomm needs to be freed */
  MPI_Comm gcomm;        /*!< A copy of the application supplied communicator */
  MPI_Comm comm;	 /*!< The current communicator */
  idx_t ncommpes;        /*!< The maximum number of processors that a processor 
                              may need to communicate with. This determines the
                              size of the sreq/rreq/statuses arrays and is 
                              updated after every call to CommSetup() */
  MPI_Request *sreq;     /*!< MPI send requests */
  MPI_Request *rreq;     /*!< MPI receive requests */
  MPI_Status *statuses;  /*!< MPI status for p2p i-messages */
  MPI_Status status;

  /* workspace variables */
  gk_mcore_t *mcore;        /* GKlib's mcore */

  /* These are for use by the k-way refinement routines */
  size_t nbrpoolsize;      /*!< The number of cnbr_t entries that have been allocated */
  size_t nbrpoolcpos;      /*!< The position of the first free entry in the array */
  size_t nbrpoolreallocs;  /*!< The number of times the pool was resized */

  cnbr_t *cnbrpool;     /*!< The pool of cnbr_t entries to be used during refinement.
                             The size and current position of the pool is controlled
                             by nnbrs & cnbrs */


  /* Various Timers */
  timer TotalTmr, InitPartTmr, MatchTmr, ContractTmr, CoarsenTmr, RefTmr,
        SetupTmr, ProjectTmr, KWayInitTmr, KWayTmr, MoveTmr, RemapTmr, 
        SerialTmr, AuxTmr1, AuxTmr2, AuxTmr3, AuxTmr4, AuxTmr5, AuxTmr6;
} ctrl_t;



/*************************************************************************
* The following data structure stores a mesh.
**************************************************************************/
typedef struct mesh_t {
  idx_t etype;
  idx_t gnelms, gnns;
  idx_t nelms, nns;
  idx_t ncon;
  idx_t esize, gminnode;
  idx_t *elmdist;
  idx_t *elements;
  idx_t *elmwgt;
} mesh_t;

