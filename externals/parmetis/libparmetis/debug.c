/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * debug.c
 *
 * This file contains various functions that are used to display debuging 
 * information
 *
 * Started 10/20/96
 * George
 *
 * $Id: debug.c 10391 2011-06-23 19:00:08Z karypis $
 *
 */

#include <parmetislib.h>


/*************************************************************************
* This function prints a vector stored in each processor 
**************************************************************************/
void PrintVector(ctrl_t *ctrl, idx_t n, idx_t first, idx_t *vec, char *title)
{
  idx_t i, penum;

  for (penum=0; penum<ctrl->npes; penum++) {
    if (ctrl->mype == penum) {
      if (ctrl->mype == 0)
        fprintf(stdout, "%s\n", title);
      fprintf(stdout, "\t%3"PRIDX". ", ctrl->mype);
      for (i=0; i<n; i++)
        fprintf(stdout, "[%"PRIDX" %"PRIDX"] ", first+i, vec[i]);
      fprintf(stdout, "\n");
      fflush(stdout);
    }
    gkMPI_Barrier(ctrl->comm);
  }
}


/*************************************************************************
* This function prints a vector stored in each processor 
**************************************************************************/
void PrintVector2(ctrl_t *ctrl, idx_t n, idx_t first, idx_t *vec, char *title)
{
  idx_t i, penum;

  for (penum=0; penum<ctrl->npes; penum++) {
    if (ctrl->mype == penum) {
      if (ctrl->mype == 0)
        printf("%s\n", title);
      printf("\t%3"PRIDX". ", ctrl->mype);
      for (i=0; i<n; i++)
        printf("[%"PRIDX" %"PRIDX".%"PRIDX"] ", first+i, 
            (idx_t)(vec[i]>=KEEP_BIT ? 1 : 0), 
            (idx_t)(vec[i]>=KEEP_BIT ? vec[i]-KEEP_BIT : vec[i]));
      printf("\n");
      fflush(stdout);
    }
    gkMPI_Barrier(ctrl->comm);
  }
}


/*************************************************************************
* This function prints a vector stored in each processor 
**************************************************************************/
void PrintPairs(ctrl_t *ctrl, idx_t n, ikv_t *pairs, char *title)
{
  idx_t i, penum;

  for (penum=0; penum<ctrl->npes; penum++) {
    if (ctrl->mype == penum) {
      if (ctrl->mype == 0)
        printf("%s\n", title);
      printf("\t%3"PRIDX". ", ctrl->mype);
      for (i=0; i<n; i++)
        printf("[%"PRIDX" %"PRIDX", %"PRIDX"] ", i, pairs[i].key, pairs[i].val);
      printf("\n");
      fflush(stdout);
    }
    gkMPI_Barrier(ctrl->comm);
  }
}



/*************************************************************************
* This function prints the local portion of the graph stored at each 
* processor
**************************************************************************/
void PrintGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, penum;
  idx_t firstvtx;

  gkMPI_Barrier(ctrl->comm);

  firstvtx = graph->vtxdist[ctrl->mype];

  for (penum=0; penum<ctrl->npes; penum++) {
    if (ctrl->mype == penum) {
      printf("\t%"PRIDX"", penum);
      for (i=0; i<graph->nvtxs; i++) {
        if (i==0)
          printf("\t%2"PRIDX" %2"PRIDX"\t", firstvtx+i, graph->vwgt[i]);
        else
          printf("\t\t%2"PRIDX" %2"PRIDX"\t", firstvtx+i, graph->vwgt[i]);
        for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
          printf("[%"PRIDX" %"PRIDX"] ", graph->adjncy[j], graph->adjwgt[j]);
        printf("\n");
      }
      fflush(stdout);
    }
    gkMPI_Barrier(ctrl->comm);
  }
}


/*************************************************************************
* This function prints the local portion of the graph stored at each 
* processor along with degree information during refinement
**************************************************************************/
void PrintGraph2(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, penum;
  idx_t firstvtx;

  gkMPI_Barrier(ctrl->comm);

  firstvtx = graph->vtxdist[ctrl->mype];

  for (penum=0; penum<ctrl->npes; penum++) {
    if (ctrl->mype == penum) {
      printf("\t%"PRIDX"", penum);
      for (i=0; i<graph->nvtxs; i++) {
        if (i==0)
          printf("\t%2"PRIDX" %2"PRIDX" [%"PRIDX" %"PRIDX" %"PRIDX"]\t", firstvtx+i, graph->vwgt[i], graph->where[i], graph->ckrinfo[i].id, graph->ckrinfo[i].ed);
        else
          printf("\t\t%2"PRIDX" %2"PRIDX" [%"PRIDX" %"PRIDX" %"PRIDX"]\t", firstvtx+i, graph->vwgt[i], graph->where[i], graph->ckrinfo[i].id, graph->ckrinfo[i].ed);
        for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
          printf("[%"PRIDX" %"PRIDX"] ", graph->adjncy[j], graph->adjwgt[j]);
        printf("\n");
      }
      fflush(stdout);
    }
    gkMPI_Barrier(ctrl->comm);
  }
}


/*************************************************************************
* This function prints the information computed during setup
**************************************************************************/
void PrintSetUpInfo(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, penum;

  gkMPI_Barrier(ctrl->comm);

  for (penum=0; penum<ctrl->npes; penum++) {
    if (ctrl->mype == penum) {
      printf("PE: %"PRIDX", nnbrs: %"PRIDX"\n", ctrl->mype, graph->nnbrs);
      printf("\tSending...\n");
      for (i=0; i<graph->nnbrs; i++) {
        printf("\t\tTo: %"PRIDX": ", graph->peind[i]);
        for (j=graph->sendptr[i]; j<graph->sendptr[i+1]; j++)
          printf("%"PRIDX" ", graph->sendind[j]);
        printf("\n");
      }
      printf("\tReceiving...\n");
      for (i=0; i<graph->nnbrs; i++) {
        printf("\t\tFrom: %"PRIDX": ", graph->peind[i]);
        for (j=graph->recvptr[i]; j<graph->recvptr[i+1]; j++)
          printf("%"PRIDX" ", graph->recvind[j]);
        printf("\n");
      }
      printf("\n");
    }
    gkMPI_Barrier(ctrl->comm);
  }

}


/*************************************************************************
* This function prints information about the graphs that were sent/received
**************************************************************************/
void PrintTransferedGraphs(ctrl_t *ctrl, idx_t nnbrs, idx_t *peind, 
         idx_t *slens, idx_t *rlens, idx_t *sgraph, idx_t *rgraph)
{
  idx_t i, ii, jj, ll, penum;

  gkMPI_Barrier(ctrl->comm);
  for (penum=0; penum<ctrl->npes; penum++) {
    if (ctrl->mype == penum) {
      printf("PE: %"PRIDX", nnbrs: %"PRIDX"", ctrl->mype, nnbrs);
      for (ll=i=0; i<nnbrs; i++) {
        if (slens[i+1]-slens[i] > 0) {
          printf("\n\tTo %"PRIDX"\t", peind[i]);
          for (ii=slens[i]; ii<slens[i+1]; ii++) {
            printf("%"PRIDX" %"PRIDX" %"PRIDX", ", sgraph[ll], sgraph[ll+1], sgraph[ll+2]);
            for (jj=0; jj<sgraph[ll+1]; jj++)
              printf("[%"PRIDX" %"PRIDX"] ", sgraph[ll+3+2*jj], sgraph[ll+3+2*jj+1]);
            printf("\n\t\t");
            ll += 3+2*sgraph[ll+1];
          }
        }
      }

      for (ll=i=0; i<nnbrs; i++) {
        if (rlens[i+1]-rlens[i] > 0) {
          printf("\n\tFrom %"PRIDX"\t", peind[i]);
          for (ii=rlens[i]; ii<rlens[i+1]; ii++) {
            printf("%"PRIDX" %"PRIDX" %"PRIDX", ", rgraph[ll], rgraph[ll+1], rgraph[ll+2]);
            for (jj=0; jj<rgraph[ll+1]; jj++)
              printf("[%"PRIDX" %"PRIDX"] ", rgraph[ll+3+2*jj], rgraph[ll+3+2*jj+1]);
            printf("\n\t\t");
            ll += 3+2*rgraph[ll+1];
          }
        }
      }
      printf("\n");
    }
    gkMPI_Barrier(ctrl->comm);
  }

}


/*************************************************************************
* This function writes a graph in the format used by serial METIS
**************************************************************************/
void WriteMetisGraph(idx_t nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *adjwgt)
{
  idx_t i, j;
  FILE *fp;

  fp = fopen("test.graph", "w");

  fprintf(fp, "%"PRIDX" %"PRIDX" 11", nvtxs, xadj[nvtxs]/2);
  for (i=0; i<nvtxs; i++) {
    fprintf(fp, "\n%"PRIDX" ", vwgt[i]);
    for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fp, " %"PRIDX" %"PRIDX"", adjncy[j]+1, adjwgt[j]);
  }
  fclose(fp);
}

