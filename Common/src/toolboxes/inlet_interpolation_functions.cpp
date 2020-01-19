/*!
 * \file inlet_interpolation_functions.cpp
 * \brief Inlet_interpolation_functions
 * \author Aman Baig
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */
 

#include "../../include/toolboxes/inlet_interpolation_functions.hpp"


CInletInterpolation::CInletInterpolation(CGeometry **geometry, CConfig *config,string profile_filename,unsigned short KIND_MARKER, unsigned short iMarker, unsigned short jMarker, unsigned short nDim){

  nColumns = profileReader.GetNumberOfColumnsInProfile(jMarker);
  CMarkerProfileReaderFVM profileReader(geometry[MESH_0], config, profile_filename, KIND_MARKER, 52);
  nRow = profileReader.GetNumberOfRowsInProfile(jMarker);
  Inlet_Data = profileReader.GetDataForProfile(jMarker);
  Inlet_Values.resize(nColumns + nDim);
  InletInterpolatedData.resize((nColumns + nDim)*geometry[MESH_0]->nVertex[iMarker]);

  SetVertex();

  if(config->GetPrintInlet_InterpolatedData() == true)
    PrintInterpolatedData();

}


void CInletInterpolation::Interpolate()
{
  switch(config->GetKindInletInterpolationFunction()){
    case(ONED_LINEAR_SPANWISE):
      LinearInterpolation();
    case(ONED_AKIMASPLINE_SPANWISE):
      AkimaInterpolation(iRow);
  }
}


void CInletInterpolation::SetVertex(){

  for (iVertex = 0; iVertex < geometry[MESH_0]->nVertex[iMarker]; iVertex++){
    iPoint   = geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
    Coord    = geometry[MESH_0]->node[iPoint]->GetCoord();  

    Interp_Radius = sqrt(pow(Coord[0],2)+ pow(Coord[1],2));
    Theta = atan2(Coord[1],Coord[0]);  

    for  (unsigned short iVar=0; iVar < nDim; iVar++)
      Inlet_Values[iVar]=Coord[iVar];

    for (iRow = 0; iRow < profileReader.GetNumberOfRowsInProfile(jMarker); iRow++) {
    /*--- Find the two closest radii for the specific vertex. ---*/
      if (Inlet_Data[nColumns*iRow] <= Interp_Radius && Inlet_Data[nColumns*(iRow+1)] >= Interp_Radius){

      for (index=1; index<nColumns; index++){
      Point_Match = true;

      Interpolate();
      
      if (index > nDim+1)
        Inlet_Values[index+(nDim-1)+1]=interpolated_value;
      else
        Inlet_Values[index+(nDim-1)]=interpolated_value;
      }

      /*--- New interpolated parameters for that Interp_Radius. ---*/
      Parameter1=Inlet_Values[nDim+2];
      Parameter2=Inlet_Values[nDim+3];

      CorrectForInterpolationType();

      Inlet_Values[nDim+2] = unit_r*cos(Theta) - unit_Theta*sin(Theta); //for ix
      Inlet_Values[nDim+3] = unit_r*sin(Theta) + unit_Theta*cos(Theta); //for iy
      Inlet_Values[nDim+4] = sqrt(1-pow(unit_r,2)- pow(unit_Theta,2));  //for iz

      SetInterpolatedData();
      }
    }

    if(Point_Match == false){
      cout<<"No matching radius found for the Point: "<<iPoint<<": "<<"Node: "<<Coord[0]<<" , "<<Coord[1];
    if(nDim==3)
      cout<<" , "<<Coord[2]<<endl;
    else
      cout<<endl;
    }

  }
}

void CInletInterpolation::CorrectForInterpolationType()
{
    switch(config->GetKindInletInterpolationType()){

    case(VR_VTHETA):
      unit_r=Parameter1;
      unit_Theta=Parameter2;

    break;
    
    case(ALPHA_PHI):
      Alpha=Parameter1;
      Phi=Parameter2;
      unit_m = sqrt(1/(1+pow(tan(Alpha),2)));

      unit_Theta = tan(Alpha)*unit_m;
      unit_r=unit_m*sin(Phi);

    break;
  }
}

void CInletInterpolation::SetInterpolatedData()
{
    for  (unsigned short iVar=0; iVar < (nColumns+nDim); iVar++)
      InletInterpolatedData[iVertex*(nColumns+nDim)+iVar]=Inlet_Values[iVar];
}

void CInletInterpolation::LinearInterpolation(){

  slope = (Inlet_Data[index+nColumns*(iRow+1)]-Inlet_Data[index+nColumns*iRow])/(Inlet_Data[nColumns*(iRow+1)]-Inlet_Data[nColumns*iRow]);

  interpolated_value = slope*(Interp_Radius - Inlet_Data[nColumns*iRow]) + Inlet_Data[index+nColumns*iRow];;

}

void CInletInterpolation::PrintInterpolatedData()
{

ofstream myfile;
myfile.precision(16);
myfile.open("Interpolated_Data_"+profileReader.GetTagForProfile(jMarker)+".dat",ios_base::out);


if (myfile.is_open())
{  
  for (unsigned long iVertex = 0; iVertex < geometry[MESH_0]->nVertex[iMarker]; iVertex++) {

            for  (unsigned long iVar=0; iVar < (nColumns+nDim); iVar++)
              {
                myfile<<InletInterpolatedData[iVertex*(nColumns+nDim)+iVar]<<"\t";
              }
            myfile<<endl;
          }
  myfile.close();
}
  else
    cout<<"file cannot be opened"<<endl;
}


void CInletInterpolation::AkimaInterpolation(unsigned long iRow_Akima)
{
    /*--- Finding the cofficients of the third order polynomial for Akima Interpolation ---*/
    dxi = Inlet_Data[nColumns*(iRow_Akima+1)] - Inlet_Data[nColumns*iRow_Akima];
    ai = Inlet_Data[nColumns*iRow_Akima + index];
    bi = Get_Ai_dash(iRow_Akima);
    ci = (3*Get_Pi(iRow_Akima)-2*bi-Get_Ai_dash(iRow_Akima+1))/dxi;
    di = (bi + Get_Ai_dash(iRow_Akima+1) - 2*Get_Pi(iRow_Akima))/pow(dxi,2);
    delta = Interp_Radius - Inlet_Data[nColumns*iRow_Akima];    

    interpolated_value = ai+bi*delta+ci*pow(delta,2)+di*pow(delta,3);
}

su2double CInletInterpolation::Get_Ai_dash(unsigned long iRow_Akima){
    //for Boundary conditions (two first and two last points require special definition, can vary for different codes)
    if(iRow_Akima == 0) {return Get_Pi(iRow_Akima);}
    else if (iRow_Akima == 1) {return (Get_Pi(iRow_Akima-1)+Get_Pi(iRow_Akima))/2;}
    else if (iRow_Akima == nRow-2) {return (Get_Pi(nRow-2)+Get_Pi(nRow-3))/2;}
    else if (iRow_Akima == nRow-1) {return Get_Pi( nRow-2);}
    else{
    if((Get_Wi(iRow_Akima+1)+Get_Wi(iRow_Akima-1)) != 0)
      return (Get_Wi(iRow_Akima+1)*Get_Pi( iRow_Akima-1) + Get_Wi(iRow_Akima-1)*Get_Pi(iRow_Akima))/(Get_Wi(iRow_Akima+1) + Get_Wi(iRow_Akima-1));
    else
      return ((Get_Pi(iRow_Akima-1)-Get_Pi(iRow_Akima))/2);
    }
}


su2double CInletInterpolation::Get_Pi(unsigned long iRow_Akima) {
return (Inlet_Data[nColumns*(iRow_Akima+1)+index] - Inlet_Data[nColumns*iRow_Akima+index])/(Inlet_Data[nColumns*(iRow_Akima+1)] - Inlet_Data[nColumns*iRow_Akima]);}

su2double CInletInterpolation::Get_Wi(unsigned long iRow_Akima) {return fabs(Get_Pi(iRow_Akima) - Get_Pi(iRow_Akima-1));}

CInletInterpolation::~CInletInterpolation(void){};