void MVP1(ulong nPointDomain,ulong *row_ptr,ulong *col_ind, 
             CSysVector & prod,su2double *matrix,const CSysVector & vec) {
#define NVAR 1
  su2double tmp;
  ulong nnz,row_i,row_e,mat_begin,prod_begin,vec_begin,v_i,p_i;
  int i;
  row_i=0; 
  row_e=row_ptr[1];
  prod_begin=0;
  mat_begin=0;
  tmp=0.0;
  for (nnz = 0; nnz < row_ptr[nPointDomain]; nnz ++) {
   vec_begin=col_ind[nnz];
   if(nnz >= row_e) {
    row_i++;
    row_e=row_ptr[row_i+1];
    p_i=prod_begin;
    prod[p_i] = tmp;   // store calculated values 
    tmp=0.0;
    prod_begin++;
   }
    v_i=vec_begin;
    tmp += matrix[mat_begin]*vec[v_i];
    mat_begin++;
  }
  p_i=prod_begin;
  prod[p_i] = tmp;   // store last calculated values 
#undef NVAR
}

#pragma GCC optimize("4")
void MVP5(ulong nPointDomain,ulong *row_ptr,ulong *col_ind, 
             CSysVector & prod,su2double *matrix,const CSysVector & vec) {
#define NVAR 5
  su2double tmp[NVAR];
  ulong nnz,row_i,row_e,mat_begin,prod_begin,vec_begin,v_i,p_i;
  int i,iVar,jVar;
  row_i=0; 
  row_e=row_ptr[1];
  prod_begin=0;
  mat_begin=0;
  for(i=0;i<NVAR;i++) tmp[i]=0.;     // initialize
  for (nnz = 0; nnz < row_ptr[nPointDomain]; nnz ++) {
   vec_begin=col_ind[nnz]*NVAR;
   /* step ahead to new row in prod and matrix */
   if(nnz >= row_e) {
    row_i++;
    row_e=row_ptr[row_i+1];
    p_i=prod_begin;
#pragma unroll NVAR
    for(iVar=0;iVar<NVAR;iVar++) {// store calculated values and initialize tmp
        prod[p_i] = tmp[iVar];   
        p_i++;
        tmp[iVar]=0.;
    }
    prod_begin+=NVAR;
   }
    jVar=0;
    v_i=vec_begin;
#pragma unroll NVAR
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
    jVar=1;
    v_i=vec_begin;
#pragma unroll NVAR
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
    jVar=2;
    v_i=vec_begin;
#pragma unroll NVAR
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
    jVar=3;
    v_i=vec_begin;
#pragma unroll NVAR
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
    jVar=4;
    v_i=vec_begin;
#pragma unroll NVAR
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
  }
  p_i=prod_begin;
#pragma unroll NVAR
  for(iVar=0;iVar<NVAR;iVar++) {
      prod[p_i] = tmp[iVar];   // store last calculated values 
      p_i++;
  }
#undef NVAR
}

void MVP4(ulong nPointDomain,ulong *row_ptr,ulong *col_ind, 
             CSysVector & prod,su2double *matrix,const CSysVector & vec) {
#define NVAR 4
  su2double tmp[NVAR];
  ulong nnz,row_i,row_e,mat_begin,prod_begin,vec_begin,v_i,p_i;
  int i,iVar,jVar;
  row_i=0; 
  row_e=row_ptr[1];
  prod_begin=0;
  mat_begin=0;
  for(i=0;i<NVAR;i++) tmp[i]=0.;     // initialize
  for (nnz = 0; nnz < row_ptr[nPointDomain]; nnz ++) {
   vec_begin=col_ind[nnz]*NVAR;
   /* step ahead to new row in prod and matrix */
   if(nnz >= row_e) {
    row_i++;
    row_e=row_ptr[row_i+1];
    p_i=prod_begin;
    for(iVar=0;iVar<NVAR;iVar++) {// store calculated values and initialize tmp
        prod[p_i] = tmp[iVar];   
        p_i++;
        tmp[iVar]=0.;
    }
    prod_begin+=NVAR;
   }
    jVar=0;
    v_i=vec_begin;
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
    jVar=1;
    v_i=vec_begin;
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
    jVar=2;
    v_i=vec_begin;
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
    jVar=3;
    v_i=vec_begin;
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
  }
  p_i=prod_begin;
  for(iVar=0;iVar<NVAR;iVar++) {
      prod[p_i] = tmp[iVar];   // store last calculated values 
      p_i++;
  }
#undef NVAR
}

void MVP3(ulong nPointDomain,ulong *row_ptr,ulong *col_ind, 
             CSysVector & prod,su2double *matrix,const CSysVector & vec) {
#define NVAR 3
  su2double tmp[NVAR];
  ulong nnz,row_i,row_e,mat_begin,prod_begin,vec_begin,v_i,p_i;
  int i,iVar,jVar;
  row_i=0; 
  row_e=row_ptr[1];
  prod_begin=0;
  mat_begin=0;
  for(i=0;i<NVAR;i++) tmp[i]=0.;     // initialize
  for (nnz = 0; nnz < row_ptr[nPointDomain]; nnz ++) {
   vec_begin=col_ind[nnz]*NVAR;
   /* step ahead to new row in prod and matrix */
   if(nnz >= row_e) {
    row_i++;
    row_e=row_ptr[row_i+1];
    p_i=prod_begin;
    for(iVar=0;iVar<NVAR;iVar++) {// store calculated values and initialize tmp
        prod[p_i] = tmp[iVar];   
        p_i++;
        tmp[iVar]=0.;
    }
    prod_begin+=NVAR;
   }
    jVar=0;
    v_i=vec_begin;
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
    jVar=1;
    v_i=vec_begin;
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
    jVar=2;
    v_i=vec_begin;
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
  }
  p_i=prod_begin;
  for(iVar=0;iVar<NVAR;iVar++) {
      prod[p_i] = tmp[iVar];   // store last calculated values 
      p_i++;
  }
#undef NVAR
}

void MVP2(ulong nPointDomain,ulong *row_ptr,ulong *col_ind, 
             CSysVector & prod,su2double *matrix,const CSysVector & vec) {
#define NVAR 2
  su2double tmp[NVAR];
  ulong nnz,row_i,row_e,mat_begin,prod_begin,vec_begin,v_i,p_i;
  int i,iVar,jVar;
  row_i=0; 
  row_e=row_ptr[1];
  prod_begin=0;
  mat_begin=0;
  for(i=0;i<NVAR;i++) tmp[i]=0.;     // initialize
  for (nnz = 0; nnz < row_ptr[nPointDomain]; nnz ++) {
   vec_begin=col_ind[nnz]*NVAR;
   /* step ahead to new row in prod and matrix */
   if(nnz >= row_e) {
    row_i++;
    row_e=row_ptr[row_i+1];
    p_i=prod_begin;
    for(iVar=0;iVar<NVAR;iVar++) {// store calculated values and initialize tmp
        prod[p_i] = tmp[iVar];   
        p_i++;
        tmp[iVar]=0.;
    }
    prod_begin+=NVAR;
   }
    jVar=0;
    v_i=vec_begin;
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
    jVar=1;
    v_i=vec_begin;
    for(iVar=0;iVar<NVAR;iVar++) { 
     tmp[jVar] += matrix[mat_begin]*vec[v_i];
     v_i++; mat_begin++;
    }
  }
  p_i=prod_begin;
  for(iVar=0;iVar<NVAR;iVar++) {
      prod[p_i] = tmp[iVar];   // store last calculated values 
      p_i++;
  }
#undef NVAR
}
#pragma GCC reset_options
