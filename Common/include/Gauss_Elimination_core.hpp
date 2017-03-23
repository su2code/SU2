    su2double d0,d1,d2,d3,d4;
    if (nVar<6) {
#define BLOCK(i,j) block[(i)*nVar+(j)]
    d0 = 1.0000000000000000 / BLOCK(0,0);
    switch(nVar) {
    case 1:
     rhs[0] *= d0;
     break;
/*  NVAR = 2 */
    case 2:
     weight = BLOCK(1,0) * d0;
     BLOCK(1,0) = weight;
     BLOCK(1,1) -= weight*BLOCK(0,1);
     rhs[1]     -= weight*rhs[0];
     rhs[1] = rhs[1]/BLOCK(1,1);
     aux = BLOCK(0,1)*rhs[1];
     rhs[0] = (rhs[0]-aux) * d0;
     break;
    case 3:
    case 4:
    case 5: 
     weight = BLOCK(1,0) * d0;
     BLOCK(1,0)  = weight;
     BLOCK(1,1) -= weight*BLOCK(0,1);
     BLOCK(1,2) -= weight*BLOCK(0,2);
     rhs[1]     -= weight*rhs[0];
     d1 = 1.0000000000000000 / BLOCK(1,1);
     weight = BLOCK(2,0) * d0;
     BLOCK(2,0)  = weight;
     BLOCK(2,1) -= weight*BLOCK(0,1);
     BLOCK(2,2) -= weight*BLOCK(0,2);
     rhs[2]     -= weight*rhs[0];
     weight = BLOCK(2,1)*d1;
     BLOCK(2,1)  = weight;
     BLOCK(2,2) -= weight*BLOCK(1,2);
     d2 = 1.0000000000000000 / BLOCK(2,2);
     rhs[2]     -= weight*rhs[1];
     if(nVar==3) {
      rhs[2]      = rhs[2]*d2;
      aux=BLOCK(1,2)*rhs[2];
      rhs[1] = (rhs[1]-aux) * d1;
      aux=BLOCK(0,1)*rhs[1]+BLOCK(0,2)*rhs[2];
      rhs[0] = (rhs[0]-aux) * d0;
      break;
     }
    //case 4: 
     BLOCK(1,3) -= BLOCK(1,0)*BLOCK(0,3);
     BLOCK(2,3) -= BLOCK(2,0)*BLOCK(0,3);
     BLOCK(2,3) -= BLOCK(2,1)*BLOCK(1,3);
     weight = BLOCK(3,0) * d0;
     BLOCK(3,0)  = weight;
     BLOCK(3,1) -= weight*BLOCK(0,1);
     BLOCK(3,2) -= weight*BLOCK(0,2);
     BLOCK(3,3) -= weight*BLOCK(0,3);
     rhs[3]     -= weight*rhs[0];
     weight = BLOCK(3,1)*d1;
     BLOCK(3,1)  = weight;
     BLOCK(3,2) -= weight*BLOCK(1,2);
     BLOCK(3,3) -= weight*BLOCK(1,3);
     rhs[3]     -= weight*rhs[1];
     weight = BLOCK(3,2)*d2;
     BLOCK(3,2)  = weight;
     BLOCK(3,3) -= weight*BLOCK(2,3);
     d3 = 1.0000000000000000 / BLOCK(3,3);
     rhs[3]     -= weight*rhs[2];
     if(nVar==4) {
     rhs[3]      = rhs[3]*d3;
     aux=BLOCK(2,3)*rhs[3];
     rhs[2] = (rhs[2]-aux) * d2;
     aux=BLOCK(1,2)*rhs[2]+BLOCK(1,3)*rhs[3];
     rhs[1] = (rhs[1]-aux) * d1;
     aux=BLOCK(0,1)*rhs[1]+BLOCK(0,2)*rhs[2]+BLOCK(0,3)*rhs[3];
     rhs[0] = (rhs[0]-aux) * d0;
     break;
     }
    //case 5: 
     BLOCK(1,4) -= BLOCK(1,0)*BLOCK(0,4);
     BLOCK(2,4) -= BLOCK(2,0)*BLOCK(0,4);
     BLOCK(2,4) -= BLOCK(2,1)*BLOCK(1,4);
     BLOCK(3,4) -= BLOCK(3,0)*BLOCK(0,4);
     BLOCK(3,4) -= BLOCK(3,1)*BLOCK(1,4);
     BLOCK(3,4) -= BLOCK(3,2)*BLOCK(2,4);
     weight = BLOCK(4,0) * d0;
     BLOCK(4,0)  = weight;
     BLOCK(4,1) -= weight*BLOCK(0,1);
     BLOCK(4,2) -= weight*BLOCK(0,2);
     BLOCK(4,3) -= weight*BLOCK(0,3);
     BLOCK(4,4) -= weight*BLOCK(0,4);
     rhs[4]     -= weight*rhs[0];
     weight = BLOCK(4,1)*d1;
     BLOCK(4,1)  = weight;
     BLOCK(4,2) -= weight*BLOCK(1,2);
     BLOCK(4,3) -= weight*BLOCK(1,3);
     BLOCK(4,4) -= weight*BLOCK(1,4);
     rhs[4]     -= weight*rhs[1];
     weight = BLOCK(4,2)*d2;
     BLOCK(4,2)  = weight;
     BLOCK(4,3) -= weight*BLOCK(2,3);
     BLOCK(4,4) -= weight*BLOCK(2,4);
     rhs[4]     -= weight*rhs[2];
     weight = BLOCK(4,3)*d3;
     BLOCK(4,3)  = weight;
     BLOCK(4,4) -= weight*BLOCK(3,4);
     d4 = 1.0000000000000000 / BLOCK(4,4);
     rhs[4]     -= weight*rhs[3];
     rhs[4]      = rhs[4]*d4;
     aux=BLOCK(3,4)*rhs[4];
     rhs[3] = (rhs[3]-aux) * d3;
     aux=BLOCK(2,3)*rhs[3]+BLOCK(2,4)*rhs[4];
     rhs[2] = (rhs[2]-aux) * d2;
     aux=BLOCK(1,2)*rhs[2]+BLOCK(1,3)*rhs[3]+BLOCK(1,4)*rhs[4];
     rhs[1] = (rhs[1]-aux) * d1;
     aux=BLOCK(0,1)*rhs[1]+BLOCK(0,2)*rhs[2]+BLOCK(0,3)*rhs[3]+BLOCK(0,4)*rhs[4];
     rhs[0] = (rhs[0]-aux) * d0;
    } // end switch
   } else {
    block[0] = 1.0000000000000000 / block[0];
    for (iVar = 1; iVar < (short)nVar; iVar++) {
      for (jVar = 0; jVar < iVar; jVar++) {
        weight = block[iVar*nVar+jVar] * block[jVar*nVar+jVar];
        block[iVar*nVar+jVar]=0.;
        for (kVar = jVar+1; kVar < (short)nVar; kVar++)
          block[iVar*nVar+kVar] -= weight*block[jVar*nVar+kVar];
        rhs[iVar] -= weight*rhs[jVar];
      }
      block[iVar*nVar+iVar] = 1.0000000000000000 / block[iVar*nVar+iVar];
    }
    
    /*--- Backwards substitution ---*/
    
    rhs[nVar-1] = rhs[nVar-1] * block[nVar*nVar-1];
    for (iVar = (short)nVar-2; iVar >= 0; iVar--) {
      aux = 0.0;
      for (jVar = iVar+1; jVar < (short)nVar; jVar++)
        aux += block[iVar*nVar+jVar]*rhs[jVar];
      rhs[iVar] = (rhs[iVar]-aux) * block[iVar*nVar+iVar];
      if (iVar == 0) break;
    }
   }

#undef BLOCK
