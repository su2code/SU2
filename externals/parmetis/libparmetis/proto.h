/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h 10592 2011-07-16 21:17:53Z karypis $
 *
 */

/* ctrl.c */
ctrl_t *SetupCtrl(pmoptype_et optype, idx_t *options, idx_t ncon, idx_t nparts,
            real_t *tpwgts, real_t *ubvec, MPI_Comm comm);
void SetupCtrl_invtvwgts(ctrl_t *ctrl, graph_t *graph);
void FreeCtrl(ctrl_t **r_ctrl);



/* kmetis.c */
void Global_Partition(ctrl_t *, graph_t *);

/* mmetis.c */

/* gkmetis.c */

/* match.c */
void Match_Global(ctrl_t *, graph_t *);
void Match_Local(ctrl_t *, graph_t *);
void CreateCoarseGraph_Global(ctrl_t *, graph_t *, idx_t);
void CreateCoarseGraph_Local(ctrl_t *, graph_t *, idx_t);


/* initpart.c */
void InitPartition(ctrl_t *, graph_t *);
void KeepPart(ctrl_t *, graph_t *, idx_t *, idx_t);

/* kwayrefine.c */
void ProjectPartition(ctrl_t *, graph_t *);
void ComputePartitionParams(ctrl_t *, graph_t *);
void KWayFM(ctrl_t *, graph_t *, idx_t);
void KWayBalance(ctrl_t *, graph_t *, idx_t);


/* remap.c */
void ParallelReMapGraph(ctrl_t *, graph_t *);
void ParallelTotalVReMap(ctrl_t *, idx_t *, idx_t *, idx_t, idx_t);
idx_t SimilarTpwgts(real_t *, idx_t, idx_t, idx_t);

/* move.c */
graph_t *MoveGraph(ctrl_t *, graph_t *);
/* move.c */
void CheckMGraph(ctrl_t *, graph_t *); 
void ProjectInfoBack(ctrl_t *, graph_t *, idx_t *, idx_t *);
void FindVtxPerm(ctrl_t *, graph_t *, idx_t *);

/* wspace.c */
void AllocateWSpace(ctrl_t *ctrl, size_t nwords);
void AllocateRefinementWorkSpace(ctrl_t *ctrl, idx_t nbrpoolsize);
void FreeWSpace(ctrl_t *);
void *wspacemalloc(ctrl_t *ctrl, size_t nbytes);
idx_t *iwspacemalloc(ctrl_t *ctrl, size_t n);
real_t *rwspacemalloc(ctrl_t *ctrl, size_t n);
ikv_t *ikvwspacemalloc(ctrl_t *ctrl, size_t n);
rkv_t *rkvwspacemalloc(ctrl_t *ctrl, size_t n);
void cnbrpoolReset(ctrl_t *ctrl);
idx_t cnbrpoolGetNext(ctrl_t *ctrl, idx_t nnbrs);


/* ametis.c */
void Adaptive_Partition(ctrl_t *, graph_t *);

/* rmetis.c */


/* wave.c */
real_t WavefrontDiffusion(ctrl_t *, graph_t *, idx_t *);

/* balancemylink.c */
idx_t BalanceMyLink(ctrl_t *, graph_t *, idx_t *, idx_t, idx_t, real_t *, 
          real_t, real_t *, real_t *, real_t);

/* redomylink.c */
void RedoMyLink(ctrl_t *, graph_t *, idx_t *, idx_t, idx_t, real_t *, real_t *, real_t *);

/* initbalance.c */
void Balance_Partition(ctrl_t *, graph_t *);
graph_t *AssembleAdaptiveGraph(ctrl_t *, graph_t *);

/* mdiffusion.c */
idx_t Mc_Diffusion(ctrl_t *, graph_t *, idx_t *, idx_t *, idx_t *, idx_t);
graph_t *ExtractGraph(ctrl_t *, graph_t *, idx_t *, idx_t *, idx_t *);

/* diffutil.c */
void SetUpConnectGraph(graph_t *, matrix_t *, idx_t *);
void Mc_ComputeMoveStatistics(ctrl_t *, graph_t *, idx_t *, idx_t *, idx_t *);
 idx_t Mc_ComputeSerialTotalV(graph_t *, idx_t *);
void ComputeLoad(graph_t *, idx_t, real_t *, real_t *, idx_t);
void ConjGrad2(matrix_t *, real_t *, real_t *, real_t, real_t *);
void mvMult2(matrix_t *, real_t *, real_t *);
void ComputeTransferVector(idx_t, matrix_t *, real_t *, real_t *, idx_t);
idx_t ComputeSerialEdgeCut(graph_t *);
idx_t ComputeSerialTotalV(graph_t *, idx_t *);

/* akwayfm.c */
void KWayAdaptiveRefine(ctrl_t *, graph_t *, idx_t);

/* selectq.c */
void Mc_DynamicSelectQueue(ctrl_t *ctrl, idx_t nqueues, idx_t ncon, idx_t subdomain1,
         idx_t subdomain2, idx_t *currentq, real_t *flows, idx_t *from, idx_t *qnum,
         idx_t minval, real_t avgvwgt, real_t maxdiff);
idx_t Mc_HashVwgts(ctrl_t *ctrl, idx_t ncon, real_t *nvwgt);
idx_t Mc_HashVRank(idx_t ncon, idx_t *vwgt);


/* csrmatch.c */
void CSR_Match_SHEM(matrix_t *, idx_t *, idx_t *, idx_t *, idx_t);

/* serial.c */
void Mc_ComputeSerialPartitionParams(ctrl_t *ctrl, graph_t *, idx_t);
void Mc_SerialKWayAdaptRefine(ctrl_t *ctrl, graph_t *, idx_t, idx_t *, real_t *, idx_t);
idx_t AreAllHVwgtsBelow(idx_t, real_t, real_t *, real_t, real_t *, real_t *);
void ComputeHKWayLoadImbalance(idx_t, idx_t, real_t *, real_t *);
void SerialRemap(ctrl_t *ctrl, graph_t *, idx_t, idx_t *, idx_t *, idx_t *, real_t *);
int SSMIncKeyCmp(const void *, const void *);
void Mc_Serial_FM_2WayRefine(ctrl_t *ctrl, graph_t *, real_t *, idx_t);
void Serial_SelectQueue(idx_t, real_t *, real_t *, idx_t *, idx_t *, rpq_t **[2]);
idx_t Serial_BetterBalance(idx_t, real_t *, real_t *, real_t *, real_t *);
real_t Serial_Compute2WayHLoadImbalance(idx_t, real_t *, real_t *);
void Mc_Serial_Balance2Way(ctrl_t *ctrl, graph_t *, real_t *, real_t);
void Mc_Serial_Init2WayBalance(ctrl_t *ctrl, graph_t *, real_t *);
idx_t Serial_SelectQueueOneWay(idx_t, real_t *, real_t *, idx_t, rpq_t **[2]);
void Mc_Serial_Compute2WayPartitionParams(ctrl_t *ctrl, graph_t *);
idx_t Serial_AreAnyVwgtsBelow(idx_t, real_t, real_t *, real_t, real_t *, real_t *);

/* weird.c */
int CheckInputsPartKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
        idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts,
        real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part,
        MPI_Comm *comm);
int CheckInputsPartGeomKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
        idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ndims, real_t *xyz, 
        idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
        idx_t *edgecut, idx_t *part, MPI_Comm *comm);
int CheckInputsPartGeom(idx_t *vtxdist, idx_t *ndims, real_t *xyz, 
        idx_t *part, MPI_Comm *comm);
int CheckInputsAdaptiveRepart(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy,
        idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *wgtflag, 
        idx_t *numflag, idx_t *ncon, idx_t *nparts, real_t *tpwgts, 
        real_t *ubvec, real_t *ipc2redist, idx_t *options, idx_t *edgecut, 
        idx_t *part, MPI_Comm *comm);
int CheckInputsNodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, 
        idx_t *numflag, idx_t *options, idx_t *order, idx_t *sizes,
        MPI_Comm *comm);
int CheckInputsPartMeshKway(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt,
        idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *ncommon, idx_t *nparts,
        real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part,
        MPI_Comm *comm);
void PartitionSmallGraph(ctrl_t *, graph_t *);


/* mesh.c */

/* pspases.c */
graph_t *AssembleEntireGraph(ctrl_t *, idx_t *, idx_t *, idx_t *);

/* node_refine.c */
void AllocateNodePartitionParams(ctrl_t *, graph_t *);
void ComputeNodePartitionParams(ctrl_t *, graph_t *);
void UpdateNodePartitionParams(ctrl_t *, graph_t *);
void KWayNodeRefine_Greedy(ctrl_t *ctrl, graph_t *graph, idx_t npasses, real_t ubfrac);
void KWayNodeRefine2Phase(ctrl_t *ctrl, graph_t *graph, idx_t npasses, real_t ubfrac);
void KWayNodeRefineInterior(ctrl_t *ctrl, graph_t *graph, idx_t npasses, real_t ubfrac);
void PrintNodeBalanceInfo(ctrl_t *, idx_t, idx_t *, idx_t *, char *);


/* initmsection.c */
void InitMultisection(ctrl_t *, graph_t *);
graph_t *AssembleMultisectedGraph(ctrl_t *, graph_t *);


/* ometis.c */
void MultilevelOrder(ctrl_t *ctrl, graph_t *graph, idx_t *order, idx_t *sizes);
void Order_Partition_Multiple(ctrl_t *ctrl, graph_t *graph);
void Order_Partition(ctrl_t *ctrl, graph_t *graph, idx_t *nlevels, idx_t clevel);
void LabelSeparators(ctrl_t *, graph_t *, idx_t *, idx_t *, idx_t *, idx_t *);
void CompactGraph(ctrl_t *, graph_t *, idx_t *);
void LocalNDOrder(ctrl_t *, graph_t *, idx_t *, idx_t);


/* xyzpart.c */
void Coordinate_Partition(ctrl_t *, graph_t *, idx_t, real_t *, idx_t);
void IRBinCoordinates(ctrl_t *ctrl, graph_t *graph, idx_t ndims, real_t *xyz, 
         idx_t nbins, idx_t *bxyz);
void RBBinCoordinates(ctrl_t *ctrl, graph_t *graph, idx_t ndims, real_t *xyz, 
         idx_t nbins, idx_t *bxyz);
void SampleSort(ctrl_t *, graph_t *, ikv_t *);
void PseudoSampleSort(ctrl_t *, graph_t *, ikv_t *);


/* stat.c */
void ComputeSerialBalance(ctrl_t *, graph_t *, idx_t *, real_t *);
void ComputeParallelBalance(ctrl_t *, graph_t *, idx_t *, real_t *);
void Mc_PrintThrottleMatrix(ctrl_t *, graph_t *, real_t *);
void PrintPostPartInfo(ctrl_t *ctrl, graph_t *graph, idx_t movestats);
void ComputeMoveStatistics(ctrl_t *, graph_t *, idx_t *, idx_t *, idx_t *);

/* debug.c */
void PrintVector(ctrl_t *, idx_t, idx_t, idx_t *, char *);
void PrintVector2(ctrl_t *, idx_t, idx_t, idx_t *, char *);
void PrintPairs(ctrl_t *, idx_t, ikv_t *, char *);
void PrintGraph(ctrl_t *, graph_t *);
void PrintGraph2(ctrl_t *, graph_t *);
void PrintSetUpInfo(ctrl_t *ctrl, graph_t *graph);
void PrintTransferedGraphs(ctrl_t *, idx_t, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *);
void WriteMetisGraph(idx_t, idx_t *, idx_t *, idx_t *, idx_t *);


/* comm.c */
void CommSetup(ctrl_t *, graph_t *);
void CommUpdateNnbrs(ctrl_t *ctrl, idx_t nnbrs);
void CommInterfaceData(ctrl_t *ctrl, graph_t *graph, idx_t *data, idx_t *recvvector);
void CommChangedInterfaceData(ctrl_t *ctrl, graph_t *graph, idx_t nchanged,
         idx_t *changed, idx_t *data, ikv_t *sendpairs, ikv_t *recvpairs);
idx_t GlobalSEMax(ctrl_t *, idx_t);
idx_t GlobalSEMaxComm(MPI_Comm comm, idx_t value);
idx_t GlobalSEMin(ctrl_t *, idx_t);
idx_t GlobalSEMinComm(MPI_Comm comm, idx_t value);
idx_t GlobalSESum(ctrl_t *, idx_t);
idx_t GlobalSESumComm(MPI_Comm comm, idx_t value);
real_t GlobalSEMaxFloat(ctrl_t *, real_t);
real_t GlobalSEMinFloat(ctrl_t *, real_t);
real_t GlobalSESumFloat(ctrl_t *, real_t);

/* util.c */
void myprintf(ctrl_t *ctrl, char *f_str,...);
void rprintf(ctrl_t *ctrl, char *f_str,...);
void mypridx_tf(ctrl_t *, char *f_str,...);
void rpridx_tf(ctrl_t *, char *f_str,...);
idx_t BSearch(idx_t, idx_t *, idx_t);
void RandomPermute(idx_t, idx_t *, idx_t);
void FastRandomPermute(idx_t, idx_t *, idx_t);
idx_t ispow2(idx_t);
idx_t log2Int(idx_t);
void BucketSortKeysDec(idx_t, idx_t, idx_t *, idx_t *);
real_t BetterVBalance(idx_t, real_t *, real_t *, real_t *);
idx_t IsHBalanceBetterTT(idx_t, real_t *, real_t *, real_t *, real_t *);
idx_t IsHBalanceBetterFT(idx_t, real_t *, real_t *, real_t *, real_t *);
void GetThreeMax(idx_t, real_t *, idx_t *, idx_t *, idx_t *);
size_t rargmax_strd(size_t n, real_t *x, size_t incx);
size_t rargmin_strd(size_t n, real_t *x, size_t incx);
size_t rargmax2(size_t n, real_t *x);
real_t ravg(size_t n, real_t *x);
real_t rfavg(size_t n, real_t *x);

/* grsetup.c */
graph_t *SetupGraph(ctrl_t *ctrl, idx_t ncon, idx_t *vtxdist, idx_t *xadj,
             idx_t *vwgt, idx_t *vsize, idx_t *adjncy, idx_t *adjwgt, 
             idx_t wgtflag);
void SetupGraph_nvwgts(ctrl_t *ctrl, graph_t *graph);
graph_t *CreateGraph(void);
void InitGraph(graph_t *);
void FreeGraph(graph_t *graph);
void FreeNonGraphFields(graph_t *graph);
void FreeNonGraphNonSetupFields(graph_t *graph);
void FreeInitialGraphAndRemap(graph_t *graph);
void ChangeNumbering(idx_t *, idx_t *, idx_t *, idx_t *, idx_t, idx_t, idx_t);
void ChangeNumberingMesh(idx_t *elmdist, idx_t *eptr, idx_t *eind,
                         idx_t *xadj, idx_t *adjncy, idx_t *part,
			 idx_t npes, idx_t mype, idx_t from);

/* timer.c */
void InitTimers(ctrl_t *);
void PrintTimingInfo(ctrl_t *);
void PrintTimer(ctrl_t *, timer, char *);



/* parmetis.c */
void ChangeToFortranNumbering(idx_t *, idx_t *, idx_t *, idx_t, idx_t);


/* msetup.c */
mesh_t *SetUpMesh(idx_t *etype, idx_t *ncon, idx_t *elmdist, idx_t *elements,
      idx_t *elmwgt, idx_t *wgtflag, MPI_Comm *comm);
mesh_t *CreateMesh(void);
void InitMesh(mesh_t *mesh);


/* gkmpi.c */
int gkMPI_Comm_size(MPI_Comm comm, idx_t *size);
int gkMPI_Comm_rank(MPI_Comm comm, idx_t *rank);
int gkMPI_Get_count(MPI_Status *status, MPI_Datatype datatype,
        idx_t *count);
int gkMPI_Send(void *buf, idx_t count, MPI_Datatype datatype, idx_t dest,
        idx_t tag, MPI_Comm comm);
int gkMPI_Recv(void *buf, idx_t count, MPI_Datatype datatype,
        idx_t source, idx_t tag, MPI_Comm comm, MPI_Status *status);
int gkMPI_Isend(void *buf, idx_t count, MPI_Datatype datatype, idx_t dest,
        idx_t tag, MPI_Comm comm, MPI_Request *request);
int gkMPI_Irecv(void *buf, idx_t count, MPI_Datatype datatype,
        idx_t source, idx_t tag, MPI_Comm comm, MPI_Request *request);
int gkMPI_Wait(MPI_Request *request, MPI_Status *status);
int gkMPI_Waitall(idx_t count, MPI_Request *array_of_requests, 
        MPI_Status *array_of_statuses);
int gkMPI_Barrier(MPI_Comm comm);
int gkMPI_Bcast(void *buffer, idx_t count, MPI_Datatype datatype,
        idx_t root, MPI_Comm comm);
int gkMPI_Reduce(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, idx_t root, MPI_Comm comm);
int gkMPI_Allreduce(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int gkMPI_Scan(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int gkMPI_Allgather(void *sendbuf, idx_t sendcount,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, MPI_Comm comm);
int gkMPI_Alltoall(void *sendbuf, idx_t sendcount,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, MPI_Comm comm);
int gkMPI_Alltoallv(void *sendbuf, idx_t *sendcounts,
        idx_t *sdispls, MPI_Datatype sendtype, void *recvbuf, 
        idx_t *recvcounts, idx_t *rdispls, MPI_Datatype recvtype, 
        MPI_Comm comm);
int gkMPI_Allgatherv(void *sendbuf, idx_t sendcount, MPI_Datatype sendtype, 
        void *recvbuf, idx_t *recvcounts, idx_t *rdispls, 
        MPI_Datatype recvtype, MPI_Comm comm);
int gkMPI_Scatterv(void *sendbuf, idx_t *sendcounts, idx_t *sdispls,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, idx_t root, MPI_Comm comm);
int gkMPI_Gatherv(void *sendbuf, idx_t sendcount, MPI_Datatype sendtype,
        void *recvbuf, idx_t *recvcounts, idx_t *displs, MPI_Datatype recvtype,
        idx_t root, MPI_Comm comm);
int gkMPI_Comm_split(MPI_Comm comm, idx_t color, idx_t key,
        MPI_Comm *newcomm);
int gkMPI_Comm_free(MPI_Comm *comm);
int gkMPI_Finalize();



