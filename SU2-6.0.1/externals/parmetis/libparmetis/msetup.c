/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * msetup.c
 *
 * This file contain various routines for setting up a mesh
 *
 * Started 10/19/96
 * George
 *
 * $Id: msetup.c 10057 2011-06-02 13:44:44Z karypis $
 *
 */

#include <parmetislib.h>



/*************************************************************************
* This function setsup the ctrl_t structure
**************************************************************************/
mesh_t *SetUpMesh(idx_t *etype, idx_t *ncon, idx_t *elmdist, idx_t *elements,
  idx_t *elmwgt, idx_t *wgtflag, MPI_Comm *comm)
{
  mesh_t *mesh;
  idx_t i, npes, mype;
  idx_t esizes[5] = {-1, 3, 4, 8, 4};
  idx_t maxnode, gmaxnode, minnode, gminnode;

  gkMPI_Comm_size(*comm, &npes);
  gkMPI_Comm_rank(*comm, &mype);

  mesh = CreateMesh();
  mesh->elmdist = elmdist;
  mesh->gnelms = elmdist[npes];
  mesh->nelms = elmdist[mype+1]-elmdist[mype];
  mesh->elements = elements;
  mesh->elmwgt = elmwgt;
  mesh->etype = *etype;
  mesh->ncon = *ncon;
  mesh->esize = esizes[*etype];

  if (((*wgtflag)&1) == 0) {
    mesh->elmwgt = ismalloc(mesh->nelms*mesh->ncon, 1, "SetUpMesh: elmwgt");
  }

  minnode = imin(mesh->nelms*mesh->esize, elements);
  gkMPI_Allreduce((void *)&minnode, (void *)&gminnode, 1, IDX_T, MPI_MIN, *comm);
  for (i=0; i<mesh->nelms*mesh->esize; i++)
    elements[i] -= gminnode;
  mesh->gminnode = gminnode;

  maxnode = imax(mesh->nelms*mesh->esize, elements);
  gkMPI_Allreduce((void *)&maxnode, (void *)&gmaxnode, 1, IDX_T, MPI_MAX, *comm);
  mesh->gnns = gmaxnode+1;

  return mesh;
}

/*************************************************************************
* This function creates a mesh_t data structure and initializes
* the various fields
**************************************************************************/
mesh_t *CreateMesh(void)
{
  mesh_t *mesh;

  mesh = (mesh_t *)gk_malloc(sizeof(mesh_t), "CreateMesh: mesh");

  InitMesh(mesh);

  return mesh;
}

/*************************************************************************
* This function initializes the various fields of a mesh_t.
**************************************************************************/
void InitMesh(mesh_t *mesh)
{

  mesh->etype = -1;
  mesh->gnelms = -1;
  mesh->gnns = -1;
  mesh->nelms = -1;
  mesh->nns = -1;
  mesh->ncon = -1;
  mesh->esize = -1;
  mesh->gminnode = 0;
  mesh->elmdist = NULL;
  mesh->elements = NULL;
  mesh->elmwgt = NULL;

  return;
}

