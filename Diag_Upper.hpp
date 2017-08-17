#pragma GCC optimize("4")
void Diag_Upper5(int nVar, ulong row_i, su2double *vec, su2double *matrix,
                su2double *prod, ulong *row_ptr, ulong *col_ind) {
#define NVAR 5
 ulong mat_begin,index,j;
 int iVar;
 su2double *block;

 mat_begin = row_ptr[row_i] * NVAR * NVAR;
 for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
    if (col_ind[index] == row_i) {
      j=col_ind[index]*NVAR;
      block=&matrix[mat_begin];
      prod[0] = block[0*NVAR+0]*vec[j+0] + block[0*NVAR+1]*vec[j+1] +
         	block[0*NVAR+2]*vec[j+2] + block[0*NVAR+3]*vec[j+3] + 
		block[0*NVAR+4]*vec[j+4];
      prod[1] = block[1*NVAR+0]*vec[j+0] + block[1*NVAR+1]*vec[j+1] +
     		block[1*NVAR+2]*vec[j+2] + block[1*NVAR+3]*vec[j+3] + 
		block[1*NVAR+4]*vec[j+4];
      prod[2] = block[2*NVAR+0]*vec[j+0] + block[2*NVAR+1]*vec[j+1] +
		block[2*NVAR+2]*vec[j+2] + block[2*NVAR+3]*vec[j+3] + 
		block[2*NVAR+4]*vec[j+4];
      prod[3] = block[3*NVAR+0]*vec[j+0] + block[3*NVAR+1]*vec[j+1] +
		block[3*NVAR+2]*vec[j+2] + block[3*NVAR+3]*vec[j+3] + 
		block[3*NVAR+4]*vec[j+4];
      prod[4] = block[4*NVAR+0]*vec[j+0] + block[4*NVAR+1]*vec[j+1] +
		block[4*NVAR+2]*vec[j+2] + block[4*NVAR+3]*vec[j+3] + 
		block[4*NVAR+4]*vec[j+4];
    }
    if (col_ind[index] > row_i) {
      j=col_ind[index]*NVAR;
      block=&matrix[mat_begin];
      prod[0]-= block[0*NVAR+0]*vec[j+0] + block[0*NVAR+1]*vec[j+1] +
         	block[0*NVAR+2]*vec[j+2] + block[0*NVAR+3]*vec[j+3] + 
		block[0*NVAR+4]*vec[j+4];
      prod[1]-= block[1*NVAR+0]*vec[j+0] + block[1*NVAR+1]*vec[j+1] +
     		block[1*NVAR+2]*vec[j+2] + block[1*NVAR+3]*vec[j+3] + 
		block[1*NVAR+4]*vec[j+4];
      prod[2]-= block[2*NVAR+0]*vec[j+0] + block[2*NVAR+1]*vec[j+1] +
		block[2*NVAR+2]*vec[j+2] + block[2*NVAR+3]*vec[j+3] + 
		block[2*NVAR+4]*vec[j+4];
      prod[3]-= block[3*NVAR+0]*vec[j+0] + block[3*NVAR+1]*vec[j+1] +
		block[3*NVAR+2]*vec[j+2] + block[3*NVAR+3]*vec[j+3] + 
		block[3*NVAR+4]*vec[j+4];
      prod[4]-= block[4*NVAR+0]*vec[j+0] + block[4*NVAR+1]*vec[j+1] +
		block[4*NVAR+2]*vec[j+2] + block[4*NVAR+3]*vec[j+3] + 
		block[4*NVAR+4]*vec[j+4];
    }
    mat_begin += NVAR * NVAR;
 }
  
}
#undef NVAR

void Diag_Upper4(int nVar, ulong row_i, su2double *vec, su2double *matrix,
                su2double *prod, ulong *row_ptr, ulong *col_ind) {
#define NVAR 4
 ulong mat_begin,index,j;
 int iVar;
 su2double *block;

 mat_begin = row_ptr[row_i] * NVAR * NVAR;
 for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
    if (col_ind[index] == row_i) {
      j=col_ind[index]*NVAR;
      block=&matrix[mat_begin];
      prod[0] = block[0*NVAR+0]*vec[j+0] + block[0*NVAR+1]*vec[j+1] +
         	block[0*NVAR+2]*vec[j+2] + block[0*NVAR+3]*vec[j+3];
      prod[1] = block[1*NVAR+0]*vec[j+0] + block[1*NVAR+1]*vec[j+1] +
     		block[1*NVAR+2]*vec[j+2] + block[1*NVAR+3]*vec[j+3];
      prod[2] = block[2*NVAR+0]*vec[j+0] + block[2*NVAR+1]*vec[j+1] +
		block[2*NVAR+2]*vec[j+2] + block[2*NVAR+3]*vec[j+3];
      prod[3] = block[3*NVAR+0]*vec[j+0] + block[3*NVAR+1]*vec[j+1] +
		block[3*NVAR+2]*vec[j+2] + block[3*NVAR+3]*vec[j+3];
    }
    if (col_ind[index] > row_i) {
      j=col_ind[index]*NVAR;
      block=&matrix[mat_begin];
      prod[0]-= block[0*NVAR+0]*vec[j+0] + block[0*NVAR+1]*vec[j+1] +
         	block[0*NVAR+2]*vec[j+2] + block[0*NVAR+3]*vec[j+3];
      prod[1]-= block[1*NVAR+0]*vec[j+0] + block[1*NVAR+1]*vec[j+1] +
     		block[1*NVAR+2]*vec[j+2] + block[1*NVAR+3]*vec[j+3];
      prod[2]-= block[2*NVAR+0]*vec[j+0] + block[2*NVAR+1]*vec[j+1] +
		block[2*NVAR+2]*vec[j+2] + block[2*NVAR+3]*vec[j+3];
      prod[3]-= block[3*NVAR+0]*vec[j+0] + block[3*NVAR+1]*vec[j+1] +
		block[3*NVAR+2]*vec[j+2] + block[3*NVAR+3]*vec[j+3];
    }
    mat_begin += NVAR * NVAR;
 }
  
}
#undef NVAR

void Diag_Upper3(int nVar, ulong row_i, su2double *vec, su2double *matrix,
                su2double *prod, ulong *row_ptr, ulong *col_ind) {
#define NVAR 3
 ulong mat_begin,index,j;
 int iVar;
 su2double *block;

 mat_begin = row_ptr[row_i] * NVAR * NVAR;
 for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
    if (col_ind[index] == row_i) {
      j=col_ind[index]*NVAR;
      block=&matrix[mat_begin];
      prod[0] = block[0*NVAR+0]*vec[j+0] + block[0*NVAR+1]*vec[j+1] +
         	block[0*NVAR+2]*vec[j+2];
      prod[1] = block[1*NVAR+0]*vec[j+0] + block[1*NVAR+1]*vec[j+1] +
     		block[1*NVAR+2]*vec[j+2];
      prod[2] = block[2*NVAR+0]*vec[j+0] + block[2*NVAR+1]*vec[j+1] +
		block[2*NVAR+2]*vec[j+2];
    }
    if (col_ind[index] > row_i) {
      j=col_ind[index]*NVAR;
      block=&matrix[mat_begin];
      prod[0]-= block[0*NVAR+0]*vec[j+0] + block[0*NVAR+1]*vec[j+1] +
         	block[0*NVAR+2]*vec[j+2];
      prod[1]-= block[1*NVAR+0]*vec[j+0] + block[1*NVAR+1]*vec[j+1] +
     		block[1*NVAR+2]*vec[j+2];
      prod[2]-= block[2*NVAR+0]*vec[j+0] + block[2*NVAR+1]*vec[j+1] +
		block[2*NVAR+2]*vec[j+2];
    }
    mat_begin += NVAR * NVAR;
 }
  
}
#undef NVAR
#pragma GCC reset_options

void Diag_Upper2(int nVar, ulong row_i, su2double *vec, su2double *matrix,
                su2double *prod, ulong *row_ptr, ulong *col_ind) {
#define NVAR 2
 ulong mat_begin,index,j;
 int iVar;
 su2double *block;

 mat_begin = row_ptr[row_i] * NVAR * NVAR;
 for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
    if (col_ind[index] == row_i) {
      j=col_ind[index]*NVAR;
      block=&matrix[mat_begin];
      prod[0] = block[0*NVAR+0]*vec[j+0] + block[0*NVAR+1]*vec[j+1];
      prod[1] = block[1*NVAR+0]*vec[j+0] + block[1*NVAR+1]*vec[j+1];
    }
    if (col_ind[index] > row_i) {
      j=col_ind[index]*NVAR;
      block=&matrix[mat_begin];
      prod[0]-= block[0*NVAR+0]*vec[j+0] + block[0*NVAR+1]*vec[j+1];
      prod[1]-= block[1*NVAR+0]*vec[j+0] + block[1*NVAR+1]*vec[j+1];
    }
    mat_begin += NVAR * NVAR;
 }
  
}
#undef NVAR
void Diag_Upper1(int nVar, ulong row_i, su2double *vec, su2double *matrix,
                su2double *prod, ulong *row_ptr, ulong *col_ind) {
#define NVAR 1
 ulong mat_begin,index,j;
 int iVar;
 su2double *block;

 mat_begin = row_ptr[row_i] * NVAR * NVAR;
 for (index = row_ptr[row_i]; index < row_ptr[row_i+1]; index++) {
    if (col_ind[index] == row_i) {
      j=col_ind[index]*NVAR;
      block=&matrix[mat_begin];
      prod[0] = block[0*NVAR+0]*vec[j+0];
    }
    if (col_ind[index] > row_i) {
      j=col_ind[index]*NVAR;
      block=&matrix[mat_begin];
      prod[0]-= block[0*NVAR+0]*vec[j+0];
    }
    mat_begin += NVAR * NVAR;
 }
  
}
#undef NVAR
