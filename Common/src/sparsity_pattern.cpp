#include "../include/sparsity_pattern.hpp"


CSparsityPattern::CSparsityPattern(){
  nnz = 0;
  row_ptr = NULL;
  col_ind = NULL;
}

CSparsityPattern::~CSparsityPattern(){}

unsigned long CSparsityPattern::GetIndex(const unsigned long iPoint, const unsigned long jPoint){
  
  unsigned long step = 0, index;
  bool found_index = false;
  
  for (index = row_ptr[iPoint]; index < row_ptr[iPoint+1]; index++) {
    step++;
    if (col_ind[index] == jPoint) return row_ptr[iPoint]+step-1;
  }
  
  /*--- Make sure the index has been found ---*/
  
  assert(found_index);
  
  return 0;
}

CDualMeshSparsity::CDualMeshSparsity(CGeometry *geometry, unsigned short fill_in){
  
  unsigned long iPoint, index;
  unsigned short *nNeigh, iNeigh;
  vector<unsigned long> vneighs;
  vector<unsigned long>::iterator it;
  nPoint = geometry->GetnPoint();
  
  /*--- Compute the number of neighbors ---*/
  
  nNeigh = new unsigned short [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nNeigh[iPoint] = (geometry->node[iPoint]->GetnPoint()+1);  // +1 -> to include diagonal element
  }  
  
  row_ptr = new unsigned long [nPoint+1];
  row_ptr[0] = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    row_ptr[iPoint+1] = row_ptr[iPoint] + nNeigh[iPoint];
  nnz = row_ptr[nPoint];
  
  /*--- Create col_ind structure ---*/
  
  col_ind = new unsigned long [nnz];
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    vneighs.clear();
    
    SetNeighbours(geometry, iPoint, 0, fill_in, vneighs);
    
    sort(vneighs.begin(), vneighs.end());
    it = unique(vneighs.begin(), vneighs.end());
    vneighs.resize( it - vneighs.begin() );
    
    index = row_ptr[iPoint];
    for (iNeigh = 0; iNeigh < vneighs.size(); iNeigh++) {
      col_ind[index] = vneighs[iNeigh];
      index++;
    }
  }
  
  delete [] nNeigh;
  
}

CDualMeshSparsity::~CDualMeshSparsity(){
  
  if (row_ptr != NULL) delete [] row_ptr;
  if (col_ind != NULL) delete [] col_ind;
  
}

void CDualMeshSparsity::SetNeighbours(CGeometry *geometry, unsigned long iPoint, unsigned short deep_level, unsigned short fill_level, vector<unsigned long> &vneighs){
  
  unsigned long Point;
  unsigned short iNode;
  
  vneighs.push_back(iPoint);
  for (iNode = 0; iNode < geometry->node[iPoint]->GetnPoint(); iNode++) {
    Point = geometry->node[iPoint]->GetPoint(iNode);
    vneighs.push_back(Point);
    if (deep_level < fill_level) SetNeighbours(geometry, Point, deep_level+1, fill_level, vneighs);
  }
}

CPrimalMeshSparsity::CPrimalMeshSparsity(CGeometry *geometry, unsigned short fill_in){
  
  unsigned long iPoint, nPoint, index, iElem, iNode, Elem;
  unsigned short *nNeigh, iNeigh;
  vector<unsigned long> vneighs;
  vector<unsigned long>::iterator it;
  nPoint = geometry->GetnPoint();
  
  for (iElem = 0; iElem < geometry->node[iPoint]->GetnElem(); iElem++) {
    Elem =  geometry->node[iPoint]->GetElem(iElem);
    for (iNode = 0; iNode < geometry->elem[Elem]->GetnNodes(); iNode++)
      vneighs.push_back(geometry->elem[Elem]->GetNode(iNode));
  }
  vneighs.push_back(iPoint);
  
  sort(vneighs.begin(), vneighs.end());
  it = unique(vneighs.begin(), vneighs.end());
  vneighs.resize(it - vneighs.begin());
  nNeigh[iPoint] = vneighs.size();
  
  /*--- Create row_ptr structure, using the number of neighbors ---*/
  
  row_ptr = new unsigned long [nPoint+1];
  row_ptr[0] = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    row_ptr[iPoint+1] = row_ptr[iPoint] + nNeigh[iPoint];
  nnz = row_ptr[nPoint];
  
  /*--- Create col_ind structure ---*/
  
  col_ind = new unsigned long [nnz];
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    vneighs.clear();
    
    SetNeighbours(geometry, iPoint, 0, fill_in, vneighs);
    
    vneighs.push_back(iPoint);
    sort(vneighs.begin(), vneighs.end());
    it = unique(vneighs.begin(), vneighs.end());
    vneighs.resize( it - vneighs.begin() );
    
    index = row_ptr[iPoint];
    for (iNeigh = 0; iNeigh < vneighs.size(); iNeigh++) {
      col_ind[index] = vneighs[iNeigh];
      index++;
    }
  }
  
  delete [] nNeigh;
}

CPrimalMeshSparsity::~CPrimalMeshSparsity(){
  
  if (row_ptr != NULL) delete [] row_ptr;
  if (col_ind != NULL) delete [] col_ind;

}

void CPrimalMeshSparsity::SetNeighbours(CGeometry *geometry, unsigned long iPoint, unsigned short deep_level, unsigned short fill_level, vector<unsigned long> &vneighs){
  
  unsigned long Point, iElem, Elem;
  unsigned short iNode;
  
  for (iElem = 0; iElem < geometry->node[iPoint]->GetnElem(); iElem++) {
    Elem =  geometry->node[iPoint]->GetElem(iElem);
    for (iNode = 0; iNode < geometry->elem[Elem]->GetnNodes(); iNode++) {
      Point = geometry->elem[Elem]->GetNode(iNode);
      vneighs.push_back(Point);
      if (deep_level < fill_level) SetNeighbours(geometry, Point, deep_level+1, fill_level, vneighs);
    }
  }
}

CILU0Sparsity::CILU0Sparsity(CGeometry *geometry, CConfig *config, CSparsityPattern *sparsity_base){
  
  /*--- For ILU0 the sparsity pattern is the same as the original matrix A, so we can simply copy the pointers ---*/
  
  nnz     = sparsity_base->nnz;
  col_ind = sparsity_base->col_ind;
  row_ptr = sparsity_base->row_ptr;
  
}

CILU0Sparsity::~CILU0Sparsity(){}

