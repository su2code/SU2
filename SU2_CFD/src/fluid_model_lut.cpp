/*!
 * fluid_model_lut.cpp
 * \brief Source of the look-up table model.
 * \author S. Vitale, A. Rubino
 * \version 4.1.2 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/fluid_model_lut.hpp"

CThermoList::CThermoList(){

	StaticEnergy = 0.0;
	Entropy      = 0.0;
	Enthalpy     = 0.0;
	Density      = 0.0;
	Pressure     = 0.0;
	SoundSpeed2  = 0.0;
	Temperature  = 0.0;
	dPdrho_e     = 0.0;
	dPde_rho     = 0.0;
	dTdrho_e     = 0.0;
	dTde_rho     = 0.0;
	Cp           = 0.0;
	Mu 			 = 0.0;
	dmudrho_T    = 0.0;
	dmudT_rho    = 0.0;
	Kt           = 0.0;
	dktdrho_T    = 0.0;
	dktdT_rho    = 0.0;

}

CThermoList::~CThermoList(){

}

//void CThermoList::CTLprint()
//{
//	cout<<"StaticEnergy:"<<StaticEnergy<<endl;
//	cout<<"Enthalpy    :"<<Enthalpy<<endl;
//	cout<<"Entropy     :"<<Entropy<<endl;
//	cout<<"Density     :"<<Density<<endl;
//	cout<<"Pressure    :"<<Pressure<<endl;
//	cout<<"SoundSpeed2 :"<<SoundSpeed2<<endl;
//	cout<<"Temperature :"<<Temperature<<endl;
//	cout<<"dPdrho_e    :"<<dPdrho_e<<endl;
//	cout<<"dPde_rho    :"<<dPde_rho<<endl;
//	cout<<"dTdrho_e    :"<<dTdrho_e<<endl;
//	cout<<"dTde_rho    :"<<dTde_rho<<endl;
//	cout<<"Cp          :"<<Cp<<endl;
//	cout<<"Mu          :"<<Mu<<endl;
//	cout<<"dmudrho_T   :"<<dmudrho_T<<endl;
//	cout<<"dmudT_rho   :"<<dmudT_rho<<endl;
//	cout<<"Kt          :"<<Kt<<endl;
//	cout<<"dktdrho_T   :"<<dktdrho_T<<endl;
//	cout<<"dktdT_rho   :"<<dktdT_rho<<endl;
//}

CLookUpTable::CLookUpTable() : CFluidModel() {

	ThermoTables = NULL;
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			coeff[i][j] = 0.0;
		}
	}
	HS_tree= NULL;
	iIndex = -1;
	jIndex = -1;
	p_dim  = 0;
	rho_dim  = 0;
}


CLookUpTable::CLookUpTable(CConfig *config ):CFluidModel() {

	ThermoTables = NULL;
//	CThermoList interpolated;
	TableLoadCFX(config->GetLUTFileName());
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			coeff[i][j] = -1.0;
		}
	}
	iIndex = -1;//negative number means it hasn't been preset yet
	jIndex = -1;//same
	//Set the nearest neighbours to -1
	for (int i=0; i<4; i++)
	{
		NN_i[i] = -1;
		NN_j[i] = -1;
	}

	cout<<"p_dim  : "<<p_dim<<endl;
	cout<<"rho_dim: "<<rho_dim<<endl;
	cout<<"Building HS_tree"<<endl;
	su2double* xtemp = new su2double[rho_dim*p_dim];
	su2double* ytemp = new su2double[rho_dim*p_dim];
	int* itemp 		 = new int[rho_dim*p_dim];
	for (int i=0; i<rho_dim; i++)
	{
		for (int j=0; j<p_dim; j++)
		{
			xtemp[p_dim*i+j] = ThermoTables[i][j].Enthalpy;
			ytemp[p_dim*i+j] = ThermoTables[i][j].Entropy;
			itemp[p_dim*i+j] = p_dim*i+j;
		}
	}
	HS_tree = KD_Tree(xtemp, ytemp,itemp, p_dim*rho_dim, 0);
	cout<<"HS_tree built"<<endl;
}


CLookUpTable::~CLookUpTable(void) {
	delete(ThermoTables);
}


void CLookUpTable::free_KD_tree(KD_node* root)
{
	if (root->dim>1)
	{
		free_KD_tree(root->upper);
		free_KD_tree(root->lower);
	}
	delete(root->x_values);
	delete(root->y_values);
	delete(root->i_values);
	delete(root);
}


struct KD_node* CLookUpTable::KD_Tree(su2double* x_values, su2double* y_values, int* i_values, int dim, int depth)
{

	struct KD_node *kdn = new KD_node;//(KD_node*) malloc(sizeof(KD_node));;
	kdn->x_values = x_values;
	kdn->y_values = y_values;
	kdn->i_values = i_values;
	kdn->depth    = depth;
	kdn->dim      = dim;
	if(dim>1)
	{
		su2double temp;
		int itemp = 0;
		int swaps = 0; //for bubblesort
		bool sorted = false;

		if (depth%2 ==0)
		{
			//bubble sort by xvalues
			while (not sorted)
			{
				swaps = 0;

				for (int i=0; i<dim-1; i++)
				{
					if (x_values[i]>x_values[i+1])
					{
						temp   = x_values[i];
						x_values[i] = x_values[i+1];
						x_values[i+1] = temp;

						temp   = y_values[i];
						y_values[i] = y_values[i+1];
						y_values[i+1] = temp;

						itemp   = i_values[i];
						i_values[i] = i_values[i+1];
						i_values[i+1] = itemp;
						//keep a record of the number of swaps performed
						swaps++;
					}
				}
				if (swaps==0) sorted = true;
			}
		}
		else if (depth%2 ==1)
		{
			//bubble sort by yvalues
			while (not sorted)
			{
				swaps = 0;

				for (int i=0; i<dim-1; i++)
				{
					if (y_values[i]>y_values[i+1])
					{
						temp   = y_values[i];
						y_values[i] = y_values[i+1];
						y_values[i+1] = temp;

						temp   = x_values[i];
						x_values[i] = x_values[i+1];
						x_values[i+1] = temp;

						itemp   = i_values[i];
						i_values[i] = i_values[i+1];
						i_values[i+1] = itemp;
						//keep a record of the number of swaps performed
						swaps++;
					}
				}
				if (swaps==0) sorted = true;
			}
		}
		//Create the new upper and lower arrays
		su2double* upperx = new su2double[dim/2];
		su2double* uppery = new su2double[dim/2];
		int*  	   upperi = new int[dim/2];
		su2double* lowerx = new su2double[dim-dim/2];
		su2double* lowery = new su2double[dim-dim/2];
		int*       loweri = new int[dim-dim/2];
		for (int i=dim/2;i<dim;i++)
		{
			upperx[i-dim/2] = x_values[i];
			uppery[i-dim/2] = y_values[i];
			upperi[i-dim/2] = i_values[i];
		}
		for (int i=0;i<dim/2;i++)
		{
			lowerx[i] = x_values[i];
			lowery[i] = y_values[i];
			loweri[i] = i_values[i];
		}

		kdn->upper = KD_Tree( upperx, uppery, upperi, dim/2, depth+1);
		kdn->lower = KD_Tree( lowerx, lowery, loweri, dim-dim/2, depth+1);
	}

	return kdn;
}


su2double CLookUpTable::Dist_KD_Tree (su2double x, su2double y, KD_node *branch)
{
	su2double dist;
	dist = pow((branch->x_values[branch->dim/2]-x)/x,2)\
			+ pow((branch->y_values[branch->dim/2]-y)/y,2);
	return dist;
}


void CLookUpTable::NN_KD_Tree (su2double thermo1, su2double thermo2, KD_node *root, su2double best_dist)
{
	//Recursive nearest neighbor search of the KD tree, without unwinding
	su2double dist;
	dist = Dist_KD_Tree(thermo1, thermo2, root);

	if (dist<=best_dist)
	{
		best_dist = dist;
		iIndex = root->i_values[root->dim/2]/p_dim;
		jIndex = root->i_values[root->dim/2]%p_dim;
	}


	if (root->dim>1)
	{
		if (root->depth%2==0)
		{
			if (root->x_values[root->dim/2]<=thermo1)
			{
				NN_KD_Tree(thermo1, thermo2, root->upper, best_dist);
			}
			else if (root->x_values[root->dim/2]>thermo1)
			{
				NN_KD_Tree(thermo1, thermo2, root->lower, best_dist);
			}
		}
		else if (root->depth%2==1)
		{
			if (root->y_values[root->dim/2]<=thermo2)
			{
				NN_KD_Tree (thermo1, thermo2, root->upper, best_dist);
			}
			else if (root->y_values[root->dim/2]>thermo2)
			{
				NN_KD_Tree (thermo1, thermo2, root->lower, best_dist);
			}
		}
	}

}


void CLookUpTable::NN4_KD_Tree (su2double thermo1, su2double thermo2, KD_node *root, su2double *best_dist)
{
	//Recursive search for the 4 nearest neighbors of the KD tree, with unwinding
	su2double dist;
	dist = Dist_KD_Tree(thermo1, thermo2, root);
	cout<<"Depth "<<root->depth;
	int i=0;
	while (i<4)
	{
		if(dist==best_dist[i]) i=5;
		if (dist<best_dist[i])
		{
			cout<<" i:"<<i;
			for (int j=3;j>i;j--)
			{
				cout<<" j: "<<j;
				best_dist[j] = best_dist[j-1];
				NN_i[j] = NN_i[j-1];
				NN_j[j] = NN_j[j-1];
			}
			best_dist[i] = dist;
			NN_i[i] = root->i_values[root->dim/2]/p_dim;
			NN_j[i] = root->i_values[root->dim/2]%p_dim;
			i = 4;
		}
		i++;

	}
	cout<<endl;
	cout<<"best_dist: ";
	for (int i=0; i<4;i++) cout<<best_dist[i]<<" , ";
	cout<<endl;
	if ((root->dim>1))
	{
		if (root->depth%2==0)
		{
			if (root->x_values[root->dim/2]<=thermo1)
			{
				NN4_KD_Tree(thermo1, thermo2, root->upper, best_dist);
				if (dist<best_dist[3])
				{
					cout<<"Here"<<endl;
					NN4_KD_Tree (thermo1, thermo2, root->lower, best_dist);
				}
			}
			else if (root->x_values[root->dim/2]>thermo1)
			{
				NN4_KD_Tree(thermo1, thermo2, root->lower, best_dist);
				if (dist<best_dist[3])
				{
					cout<<"Here"<<endl;
					NN4_KD_Tree (thermo1, thermo2, root->upper, best_dist);
				}
			}
		}
		else if (root->depth%2==1)
		{
			if (root->y_values[root->dim/2]<=thermo2)
			{
				NN4_KD_Tree (thermo1, thermo2, root->upper, best_dist);
				if (dist<best_dist[3])
				{
					cout<<"Here"<<endl;
					NN4_KD_Tree (thermo1, thermo2, root->lower, best_dist);
				}
			}
			else if (root->y_values[root->dim/2]>thermo2)
			{
				NN4_KD_Tree (thermo1, thermo2, root->lower, best_dist);
				if (dist<best_dist[3])
				{
					cout<<"Here"<<endl;
					NN4_KD_Tree (thermo1, thermo2, root->upper, best_dist);
				}
			}
		}
	}
}




void CLookUpTable::SetTDState_rhoe (su2double rho, su2double e ) {

	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((rho>Density_limits[1]) or (rho<Density_limits[0]))
	{
		cerr<<"RHOE Input Density out of bounds\n";
	}
	if ((e>StaticEnergy_limits[1]) or (e<StaticEnergy_limits[0]))
	{
		cerr<<"RHOE Input StaticEnergy out of bounds\n";
	}
	cout<<endl<<"rho desired : "<<rho<<endl;
	cout<<"e desired   : "<<e<<endl;

	su2double RunVal;
	unsigned int CFL    = 2;
	unsigned int LowerI;
	unsigned int UpperI;
	unsigned int LowerJ;
	unsigned int UpperJ;
	//Restart search from previously used index if it exists, else go to middle
	//		if (jIndex<0) LowerJ = 0;
	//		if (jIndex>=0) LowerJ = jIndex;
	//		if (jIndex<ceil(p_dim/2))
	//		{
	//			UpperJ = ceil(p_dim/2);
	//		}
	//		else UpperJ=p_dim-1;
	UpperJ = p_dim -1;
	LowerJ = 0;

	//Determine the I index: rho is equispaced (no restart)
	LowerI = floor((rho-Density_limits[0])/(Density_limits[1]-Density_limits[0])*(rho_dim-1));
	UpperI = LowerI + 1;

	su2double grad, x00, y00, y10, x10;

	while(UpperJ-LowerJ>1)
	{
		//Check current value
		y00 = ThermoTables[LowerI][LowerJ].StaticEnergy;
		y10 = ThermoTables[UpperI][LowerJ].StaticEnergy;
		x00 = ThermoTables[LowerI][LowerJ].Density;
		x10 = ThermoTables[UpperI][LowerJ].Density;
		//BUG FIXED: interpolate the 'e' value along the line between the two rho values
		//also applied to all other searches.
		RunVal = y00 + (y10-y00)/(x10-x00)*(rho-x00);
		grad = ThermoTables[LowerI][UpperJ].StaticEnergy-y00;
		if (grad*RunVal<=grad*e)
		{
			LowerJ = LowerJ + ceil((UpperJ-LowerJ)/CFL);
		}
		else if (grad*RunVal>grad*e)
		{
			UpperJ  = LowerJ;
			LowerJ = LowerJ/CFL;
		}
		cout<<"Here"<<endl;
		cout<<LowerJ<<"  "<<UpperJ<<endl;
	}

	iIndex = LowerI;
	jIndex = LowerJ;
//	cout<<"i "<<LowerI<<"  "<<UpperI<<endl;
//	cout<<"j "<<LowerJ<<"  "<<UpperJ<<endl;


//	cout<<"Closest fit box :"<<endl;
//	cout<<"Point i j :"<<endl;
//	ThermoTables[iIndex][jIndex].CTLprint();
//	cout<<"Point i+1 j :"<<endl;
//	ThermoTables[iIndex+1][jIndex].CTLprint();
//	cout<<"Point i j+1 :"<<endl;
//	ThermoTables[iIndex][jIndex+1].CTLprint();
//	cout<<"Point i+1 j+1 :"<<endl;
//	ThermoTables[iIndex+1][jIndex+1].CTLprint();
	//Now use the closest fit box to interpolate


	su2double x, y;
	x = rho - ThermoTables[iIndex][jIndex].Density;
	y = e - ThermoTables[iIndex][jIndex].StaticEnergy;
	//Set the nearest neigbours to the adjacent i and j vertexes
	NN_i[0]=iIndex; NN_i[1]=iIndex+1;NN_i[2]=iIndex  ; NN_i[3]=iIndex+1;
	NN_j[0]=jIndex; NN_j[1]=jIndex  ;NN_j[2]=jIndex+1; NN_j[3]=jIndex+1;
	//Determine interpolation coefficients
	Interp2D_ArbitrarySkewCoeff(x,y,"RHOE");
	cout<<"Interpolation matrix inverse \n";
	for (int j=0; j<3; j++)
	{
//		cout<<setw(15)<<coeff[j][0]<<"   "<<coeff[j][1]<<"   "<<coeff[j][2]<<endl;
	}

	StaticEnergy      = e;
	Density           = rho ;
	Entropy           = Interp2D_lin(x, y, "Entropy" );
	Pressure          = Interp2D_lin(x, y, "Pressure" );
	SoundSpeed2       = Interp2D_lin(x, y, "SoundSpeed2" );
	Temperature       = Interp2D_lin(x, y, "Temperature" );
	dPdrho_e          = Interp2D_lin(x, y, "dPdrho_e" );
	dPde_rho          = Interp2D_lin(x, y, "dPde_rho" );
	dTdrho_e          = Interp2D_lin(x, y, "dTdrho_e" );
	dTde_rho          = Interp2D_lin(x, y, "dTde_rho" );
	Cp                = Interp2D_lin(x, y, "Cp" );
	Mu                = Interp2D_lin(x, y, "Mu" );
	dmudrho_T         = Interp2D_lin(x, y, "dmudrho_T" );
	dmudT_rho         = Interp2D_lin(x, y, "dmudT_rho" );
	Kt                = Interp2D_lin(x, y, "Kt" );
	dktdrho_T         = Interp2D_lin(x, y, "dktdrho_T" );
	dktdT_rho         = Interp2D_lin(x, y, "dktdT_rho" );

//	cout<<"Interpolated fit:"<<endl;
//	CTLprint ();
	if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
	{
		cerr<<"RHOE Interpolated Density out of bounds\n";
	}
	if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
	{
		cerr<<"RHOE Interpolated Pressure out of bounds\n";
	}
}

void CLookUpTable::SetTDState_PT (su2double P, su2double T ) {
	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((P>Pressure_limits[1]) or (P<Pressure_limits[0]))
	{
		cerr<<"PT Input Pressure out of bounds\n";
	}
	if ((T>Temperature_limits[1]) or (T<Temperature_limits[0]))
	{
		cerr<<"PT Input Temperature out of bounds\n";
	}
	cout<<endl<<"P desired : "<<P<<endl;
	cout<<"T desired   : "<<T<<endl;

	su2double RunVal;
	unsigned int CFL = 2;
	unsigned int LowerI;
	unsigned int UpperI;
	unsigned int LowerJ;
	unsigned int UpperJ;
	//Restart search from previously used index if it exists, else go to middle
	//		if (iIndex<0) LowerI = 0;
	//		if (iIndex>=0) LowerI = iIndex;
	//		if  (iIndex<ceil(rho_dim/2))
	//		{
	//			UpperI = ceil(rho_dim/2); //probably can be made more efficient (smaller square)
	//		}
	//		else UpperI = rho_dim-1;
	LowerI = 0;
	UpperI = rho_dim-1;


	//Determine the J index: P is equispaced (no restart)
	LowerJ = floor((P-Pressure_limits[0])/(Pressure_limits[1]-Pressure_limits[0])*(p_dim-1));
	UpperJ = LowerJ + 1;

	//Determine the I index (for T)
	su2double grad, x00, y00, y01, x01;
	while(UpperI-LowerI>1)
	{
		y00 = ThermoTables[LowerI][LowerJ].Pressure;
		y01 = ThermoTables[LowerI][UpperJ].Pressure;
		x00 = ThermoTables[LowerI][LowerJ].Temperature;
		x01 = ThermoTables[LowerI][UpperJ].Temperature;
		grad = ThermoTables[UpperI][LowerJ].Temperature-x00;
		RunVal = x00 + (x01-x00)/(y01-y00)*(P-y00);
		if (grad*RunVal<=grad*T)
		{
			LowerI = LowerI + ceil((UpperI-LowerI)/CFL);
		}
		else if (grad*RunVal>grad*T)
		{
			UpperI = LowerI;
			LowerI = LowerI/CFL;
		}
		cout<<LowerI<<"  "<<UpperI<<endl;
	}

	iIndex = LowerI;
	jIndex = LowerJ;
	cout<<"i "<<LowerI<<"  "<<UpperI<<endl;
	cout<<"j "<<LowerJ<<"  "<<UpperJ<<endl;

//	cout<<"Closest fit box :"<<endl;
//	cout<<"Point i j :"<<endl;
//	ThermoTables[iIndex][jIndex].CTLprint();
//	cout<<"Point i+1 j :"<<endl;
//	ThermoTables[iIndex+1][jIndex].CTLprint();
//	cout<<"Point i j+1 :"<<endl;
//	ThermoTables[iIndex][jIndex+1].CTLprint();
//	cout<<"Point i+1 j+1 :"<<endl;
//	ThermoTables[iIndex+1][jIndex+1].CTLprint();
//	//Now use the closest fit box to interpolate

	su2double x, y;
	x = T - ThermoTables[iIndex][jIndex].Temperature;
	y = P - ThermoTables[iIndex][jIndex].Pressure;
	//Set the nearest neigbours to the adjacent i and j vertexes
	NN_i[0]=iIndex; NN_i[1]=iIndex+1;NN_i[2]=iIndex  ; NN_i[3]=iIndex+1;
	NN_j[0]=jIndex; NN_j[1]=jIndex  ;NN_j[2]=jIndex+1; NN_j[3]=jIndex+1;
	//Determine interpolation coefficients
	Interp2D_ArbitrarySkewCoeff(x,y,"PT");
	cout<<"Interpolation matrix inverse \n";
//	for (int j=0; j<3; j++)
//	{
//		cout<<setw(15)<<coeff[j][0]<<"   "<<coeff[j][1]<<"   "<<coeff[j][2]<<endl;
//	}

	Temperature       = T;
	Pressure          = P ;
	StaticEnergy      = Interp2D_lin(x, y, "StaticEnergy" );
	Entropy           = Interp2D_lin(x, y, "Entropy" );
	Density           = Interp2D_lin(x, y, "Density" );
	SoundSpeed2       = Interp2D_lin(x, y, "SoundSpeed2" );
	dPdrho_e          = Interp2D_lin(x, y, "dPdrho_e" );
	dPde_rho          = Interp2D_lin(x, y, "dPde_rho" );
	dTdrho_e          = Interp2D_lin(x, y, "dTdrho_e" );
	dTde_rho          = Interp2D_lin(x, y, "dTde_rho" );
	Cp                = Interp2D_lin(x, y, "Cp" );
	Mu                = Interp2D_lin(x, y, "Mu" );
	dmudrho_T         = Interp2D_lin(x, y, "dmudrho_T" );
	dmudT_rho         = Interp2D_lin(x, y, "dmudT_rho" );
	Kt                = Interp2D_lin(x, y, "Kt" );
	dktdrho_T         = Interp2D_lin(x, y, "dktdrho_T" );
	dktdT_rho         = Interp2D_lin(x, y, "dktdT_rho" );
	//Intermediate variables only needed for StandAlone version
//	su2double Density d.Density;
//	su2double Pressure = Pressure;
//	cout<<"Interpolated fit:"<<endl;
//	CTLprint ();
	if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
	{
		cerr<<"PT Interpolated Density out of bounds\n";
	}
	if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
	{
		cerr<<"PT Interpolated Pressure out of bounds\n";
	}
}


void CLookUpTable::SetTDState_Prho (su2double P, su2double rho ) {
	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((P>Pressure_limits[1]) or (P<Pressure_limits[0]))
	{
		cerr<<"PRHO Input Pressure out of bounds\n";
	}
	if ((rho>Density_limits[1]) or (rho<Density_limits[0]))
	{
		cerr<<"PRHO Input Density out of bounds\n";
	}
	cout<<endl<<"rho desired : "<<rho<<endl;
	cout<<"P desired   : "<<P<<endl;

	unsigned int LowerI;
	unsigned int UpperI;
	unsigned int LowerJ;
	unsigned int UpperJ;

	//Determine the I index: RHO is equispaced
	LowerI = floor((rho-Density_limits[0])/(Density_limits[1]-Density_limits[0])*(rho_dim-1));
	UpperI = LowerI + 1;

	//Determine the J index: P is equispaced (no restart)
	LowerJ = floor((P-Pressure_limits[0])/(Pressure_limits[1]-Pressure_limits[0])*(p_dim-1));
	UpperJ = LowerJ + 1;
	cout<<"i "<<LowerI<<"  "<<UpperI<<endl;
	cout<<"j "<<LowerJ<<"  "<<UpperJ<<endl;

	iIndex = LowerI;
	jIndex = LowerJ;

//	cout<<"Closest fit box :"<<endl;
//	cout<<"Point i j :"<<endl;
//	ThermoTables[iIndex][jIndex].CTLprint();
//	cout<<"Point i+1 j :"<<endl;
//	ThermoTables[iIndex+1][jIndex].CTLprint();
//	cout<<"Point i j+1 :"<<endl;
//	ThermoTables[iIndex][jIndex+1].CTLprint();
//	cout<<"Point i+1 j+1 :"<<endl;
//	ThermoTables[iIndex+1][jIndex+1].CTLprint();
	//Now use the closest fit box to interpolate


	su2double x, y;
	x = rho - ThermoTables[iIndex][jIndex].Density;
	y = P - ThermoTables[iIndex][jIndex].Pressure;
	//Set the nearest neigbours to the adjacent i and j vertexes
	NN_i[0]=iIndex; NN_i[1]=iIndex+1;NN_i[2]=iIndex  ; NN_i[3]=iIndex+1;
	NN_j[0]=jIndex; NN_j[1]=jIndex  ;NN_j[2]=jIndex+1; NN_j[3]=jIndex+1;
	//Determine interpolation coefficients
	Interp2D_ArbitrarySkewCoeff(x,y,"PRHO");
	cout<<"Interpolation matrix inverse \n";
//	for (int j=0; j<3; j++)
//	{
//		cout<<setw(15)<<coeff[j][0]<<"   "<<coeff[j][1]<<"   "<<coeff[j][2]<<endl;
//	}

	Pressure           = P;
	Density           = rho ;
	StaticEnergy      = Interp2D_lin(x, y, "StaticEnergy" );
	Entropy           = Interp2D_lin(x, y, "Entropy" );
	SoundSpeed2       = Interp2D_lin(x, y, "SoundSpeed2" );
	Temperature       = Interp2D_lin(x, y, "Temperature" );
	dPdrho_e          = Interp2D_lin(x, y, "dPdrho_e" );
	dPde_rho          = Interp2D_lin(x, y, "dPde_rho" );
	dTdrho_e          = Interp2D_lin(x, y, "dTdrho_e" );
	dTde_rho          = Interp2D_lin(x, y, "dTde_rho" );
	Cp                = Interp2D_lin(x, y, "Cp" );
	Mu                = Interp2D_lin(x, y, "Mu" );
	dmudrho_T         = Interp2D_lin(x, y, "dmudrho_T" );
	dmudT_rho         = Interp2D_lin(x, y, "dmudT_rho" );
	Kt                = Interp2D_lin(x, y, "Kt" );
	dktdrho_T         = Interp2D_lin(x, y, "dktdrho_T" );
	dktdT_rho         = Interp2D_lin(x, y, "dktdT_rho" );
//	//Intermediate variables only needed for StandAlone version
//	su2double Density = Density;
//	su2double Pressure = Pressure;
//	cout<<"Interpolated fit:"<<endl;
//	CTLprint ();
	if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
	{
		cerr<<"PRHO Interpolated Density out of bounds\n";
	}
	if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
	{
		cerr<<"PRHO Interpolated Pressure out of bounds\n";
	}

}

void CLookUpTable::SetEnergy_Prho (su2double P, su2double rho ) {


}

void CLookUpTable::SetTDState_hs (su2double h, su2double s ) {
	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((h>Enthalpy_limits[1]) or (h<Enthalpy_limits[0]))
	{
		cerr<<"HS Input Enthalpy out of bounds\n";
	}
	if ((s>Entropy_limits[1]) or (s<Entropy_limits[0]))
	{
		cerr<<"HS Input Entropy out of bounds\n";
	}
	cout<<endl<<"h desired : "<<h<<endl;
	cout<<"s desired   : "<<s<<endl;
	iIndex = HS_tree->i_values[HS_tree->dim/2]/p_dim;
	jIndex = HS_tree->i_values[HS_tree->dim/2]%p_dim;

	for (int i=0; i<4; i++)
	{
		NN_i[i] = -1;
		NN_j[i] = -1;
	}

	su2double best_dist[4];
	for (int i=0;i<4;i++)
	{
		best_dist[i]=1E10;
	}

	cout<<"Search the HS_tree"<<endl;
	NN4_KD_Tree(h,s,HS_tree,best_dist);
	cout<<"HS_tree searched"<<endl;

	cout<<"NNi"<<endl;
	for (int j=0;j<4;j++)
	{

		cout<<NN_i[j]<<", ";
	}
	cout<<endl;
	cout<<"NNj"<<endl;
	for (int j=0;j<4;j++)
	{

		cout<<NN_j[j]<<", ";
	}
	cout<<endl;

//	cout<<"Closest fit box :"<<endl;
//	cout<<"Point i j :"<<endl;
//	ThermoTables[NN_i[0]][NN_j[0]].CTLprint();
//	cout<<"Point i+1 j :"<<endl;
//	ThermoTables[NN_i[1]][NN_j[1]].CTLprint();
//	cout<<"Point i j+1 :"<<endl;
//	ThermoTables[NN_i[2]][NN_j[2]].CTLprint();
//	cout<<"Point i+1 j+1 :"<<endl;
//	ThermoTables[NN_i[3]][NN_j[3]].CTLprint();
	//Now use the closest fit box to interpolat
	su2double x,y;
	Entropy           = s;
//	Enthalpy          = h;
	StaticEnergy      = Interp2D_Inv_Dist("StaticEnergy", best_dist);
	Density           = Interp2D_Inv_Dist("Density", best_dist);
	Pressure          = Interp2D_Inv_Dist("Pressure", best_dist);
	SoundSpeed2       = Interp2D_Inv_Dist("SoundSpeed2", best_dist);
	Temperature       = Interp2D_Inv_Dist("Temperature", best_dist);
	dPdrho_e          = Interp2D_Inv_Dist("dPdrho_e", best_dist);
	dPde_rho          = Interp2D_Inv_Dist("dPde_rho", best_dist);
	dTdrho_e          = Interp2D_Inv_Dist("dTdrho_e", best_dist);
	dTde_rho          = Interp2D_Inv_Dist("dTde_rho", best_dist);
	Cp                = Interp2D_Inv_Dist("Cp", best_dist);
	Mu                = Interp2D_Inv_Dist("Mu", best_dist);
	dmudrho_T         = Interp2D_Inv_Dist("dmudrho_T", best_dist);
	dmudT_rho         = Interp2D_Inv_Dist("dmudT_rho", best_dist);
	Kt                = Interp2D_Inv_Dist("Kt", best_dist);
	dktdrho_T         = Interp2D_Inv_Dist("dktdrho_T", best_dist);
	dktdT_rho         = Interp2D_Inv_Dist("dktdT_rho", best_dist);
	//Intermediate variables only needed for StandAlone version
	su2double Density = Density;
	su2double Pressure = Pressure;
	cout<<"Interpolated fit:"<<endl;
//	CTLprint ();
	if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
	{
		cerr<<"HS Interpolated Density out of bounds\n";
	}
	if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
	{
		cerr<<"HS Interpolated Pressure out of bounds\n";
	}
}

void CLookUpTable::SetTDState_Ps (su2double P, su2double s )
{
	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((P>Pressure_limits[1]) or (P<Pressure_limits[0]))
	{
		cerr<<"PS Input Pressure out of bounds\n";
	}
	if ((s>Entropy_limits[1]) or (s<Entropy_limits[0]))
	{
		cerr<<"PS Input Entropy  out of bounds\n";
	}
	cout<<endl<<"P desired : "<<P<<endl;
	cout<<"s desired   : "<<s<<endl;

	su2double RunVal;
	unsigned int CFL = 2;
	unsigned int LowerI;
	unsigned int UpperI;
	unsigned int LowerJ;
	unsigned int UpperJ;
	//Restart search from previously used index if it exists, else go to middle
	//		if (iIndex<0) LowerI = 0;
	//		if (iIndex>=0) LowerI = iIndex;
	//		if  (iIndex<ceil(rho_dim/2))
	//		{
	//			UpperI = ceil(rho_dim/2); //probably can be made more efficient (smaller square)
	//		}
	//		else UpperI = rho_dim-1;
	LowerI = 0;
	UpperI = rho_dim-1;
	cout<<LowerI<<"  "<<UpperI<<endl;

	//Determine the I index: RHO is equispaced (no restart)
	LowerJ = floor((P-Pressure_limits[0])/(Pressure_limits[1]-Pressure_limits[0])*(p_dim-1));
	UpperJ = LowerJ + 1;

	//Determine the J index (for s)
	su2double grad, x00, y00, y01, x01;
	while(UpperI-LowerI>1)
	{
		y00 = ThermoTables[LowerI][LowerJ].Pressure;
		y01 = ThermoTables[LowerI][UpperJ].Pressure;
		x00 = ThermoTables[LowerI][LowerJ].Entropy;
		x01 = ThermoTables[LowerI][UpperJ].Entropy;
		grad = ThermoTables[UpperI][LowerJ].Entropy - x00;
		RunVal = x00 + (x01-x00)/(y01-y00)*(P-y00);
		if (grad*RunVal<=grad*s)
		{
			LowerI = LowerI + ceil((UpperI-LowerI)/CFL);
		}
		else if (grad*RunVal>grad*s)
		{
			UpperI = LowerI;
			LowerI = LowerI/2;
		}
		cout<<LowerI<<"  "<<UpperI<<endl;
	}


	iIndex = LowerI;
	jIndex = LowerJ;
	cout<<"i "<<LowerI<<"  "<<UpperI<<endl;
	cout<<"j "<<LowerJ<<"  "<<UpperJ<<endl;


//	cout<<"Closest fit box :"<<endl;
//	cout<<"Point i j :"<<endl;
//	ThermoTables[iIndex][jIndex].CTLprint();
//	cout<<"Point i+1 j :"<<endl;
//	ThermoTables[iIndex+1][jIndex].CTLprint();
//	cout<<"Point i j+1 :"<<endl;
//	ThermoTables[iIndex][jIndex+1].CTLprint();
//	cout<<"Point i+1 j+1 :"<<endl;
//	ThermoTables[iIndex+1][jIndex+1].CTLprint();
	//Now use the closest fit box to interpolate


	su2double x, y;
	y = P - ThermoTables[iIndex][jIndex].Pressure;
	x = s - ThermoTables[iIndex][jIndex].Entropy ;
	//Set the nearest neigbours to the adjacent i and j vertexes
	NN_i[0]=iIndex; NN_i[1]=iIndex+1;NN_i[2]=iIndex  ; NN_i[3]=iIndex+1;
	NN_j[0]=jIndex; NN_j[1]=jIndex  ;NN_j[2]=jIndex+1; NN_j[3]=jIndex+1;
	//Determine interpolation coefficients
	Interp2D_ArbitrarySkewCoeff(x,y,"PS");
	cout<<"Interpolation matrix inverse \n";
	for (int j=0; j<3; j++)
	{
//		cout<<setw(15)<<coeff[j][0]<<"   "<<coeff[j][1]<<"   "<<coeff[j][2]<<endl;
	}

	Entropy           = s;
	Pressure          = P ;
	StaticEnergy      = Interp2D_lin(x, y, "StaticEnergy" );
	Density           = Interp2D_lin(x, y, "Density" );
	SoundSpeed2       = Interp2D_lin(x, y, "SoundSpeed2" );
	Temperature       = Interp2D_lin(x, y, "Temperature" );
	dPdrho_e          = Interp2D_lin(x, y, "dPdrho_e" );
	dPde_rho          = Interp2D_lin(x, y, "dPde_rho" );
	dTdrho_e          = Interp2D_lin(x, y, "dTdrho_e" );
	dTde_rho          = Interp2D_lin(x, y, "dTde_rho" );
	Cp                = Interp2D_lin(x, y, "Cp" );
	Mu                = Interp2D_lin(x, y, "Mu" );
	dmudrho_T         = Interp2D_lin(x, y, "dmudrho_T" );
	dmudT_rho         = Interp2D_lin(x, y, "dmudT_rho" );
	Kt                = Interp2D_lin(x, y, "Kt" );
	dktdrho_T         = Interp2D_lin(x, y, "dktdrho_T" );
	dktdT_rho         = Interp2D_lin(x, y, "dktdT_rho" );
	//Intermediate variables only needed for StandAlone version
	su2double Density = Density;
	su2double Pressure = Pressure;
	cout<<"Interpolated fit:"<<endl;
//	CTLprint ();
	if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
	{
		cerr<<"PS Interpolated Density out of bounds\n";
	}
	if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
	{
		cerr<<"PS Interpolated Pressure out of bounds\n";
	}
}

void CLookUpTable::SetTDState_rhoT (su2double rho, su2double T ) {
	//Check if inputs are in total range (necessary but not sufficient condition)
	if ((rho>Density_limits[1]) or (rho<Density_limits[0]))
	{
		cerr<<"RHOT Input Density out of bounds\n";
	}
	if ((T>Temperature_limits[1]) or (T<Temperature_limits[0]))
	{
		cerr<<"RHOT Input Temperature out of bounds\n";
	}
	cout<<endl<<"rho desired : "<<rho<<endl;
	cout<<"T desired   : "<<T<<endl;

	su2double RunVal;
	unsigned int CFL= 2;
	unsigned int LowerI;
	unsigned int UpperI;
	unsigned int LowerJ;
	unsigned int UpperJ;
	//Restart search from previously used index if it exists, else go to middle
	//		if (jIndex<0) LowerJ = 0;
	//		if (jIndex>=0) LowerJ = jIndex;
	//		if (jIndex<ceil(p_dim/2))
	//		{
	//			UpperJ = ceil(p_dim/2);
	//		}
	//		else UpperJ=p_dim-1;
	LowerJ = 0;
	UpperJ = p_dim-1;

	//Determine the I index: RHO is equispaced (no restart)
	LowerI = floor((rho-Density_limits[0])/(Density_limits[1]-Density_limits[0])*(rho_dim-1));
	UpperI = LowerI + 1;

	//Determine the I index (for T)
	su2double grad, x00, y00, y10, x10;

	while(UpperJ-LowerJ>1)
	{
		//Check current value
		y00 = ThermoTables[LowerI][LowerJ].Temperature;
		y10 = ThermoTables[UpperI][LowerJ].Temperature;
		x00 = ThermoTables[LowerI][LowerJ].Density;
		x10 = ThermoTables[UpperI][LowerJ].Density;
		//BUG FIXED: interpolate the 'e' value along the line between the two rho values
		//also applied to all other searches.
		RunVal = y00 + (y10-y00)/(x10-x00)*(rho-x00);
		grad = ThermoTables[LowerI][UpperJ].Temperature-y00;
		if (grad*RunVal<grad*T)
		{
			LowerJ = LowerJ + ceil((UpperJ-LowerJ)/CFL);
		}
		else if (grad*RunVal>grad*T)
		{
			UpperJ  = LowerJ;
			LowerJ = LowerJ/CFL;
		}
		cout<<LowerJ<<"  "<<UpperJ<<endl;
	}

	iIndex = LowerI;
	jIndex = LowerJ;
	cout<<"i "<<LowerI<<"  "<<UpperI<<endl;
	cout<<"j "<<LowerJ<<"  "<<UpperJ<<endl;

//
//	cout<<"Closest fit box :"<<endl;
//	cout<<"Point i j :"<<endl;
//	ThermoTables[iIndex][jIndex].CTLprint();
//	cout<<"Point i+1 j :"<<endl;
//	ThermoTables[iIndex+1][jIndex].CTLprint();
//	cout<<"Point i j+1 :"<<endl;
//	ThermoTables[iIndex][jIndex+1].CTLprint();
//	cout<<"Point i+1 j+1 :"<<endl;
//	ThermoTables[iIndex+1][jIndex+1].CTLprint();
	//Now use the closest fit box to interpolate


	su2double x, y;
	x = rho - ThermoTables[iIndex][jIndex].Density;
	y = T - ThermoTables[iIndex][jIndex].Temperature;
	//Set the nearest neigbours to the adjacent i and j vertexes
	NN_i[0]=iIndex; NN_i[1]=iIndex+1;NN_i[2]=iIndex  ; NN_i[3]=iIndex+1;
	NN_j[0]=jIndex; NN_j[1]=jIndex  ;NN_j[2]=jIndex+1; NN_j[3]=jIndex+1;
	//Determine interpolation coefficients
	Interp2D_ArbitrarySkewCoeff(x,y,"RHOT");
	cout<<"Interpolation matrix inverse \n";
//	for (int j=0; j<3; j++)
//	{
//		cout<<setw(15)<<coeff[j][0]<<"   "<<coeff[j][1]<<"   "<<coeff[j][2]<<endl;
//	}

	Temperature       = T;
	Density           = rho ;
	StaticEnergy      = Interp2D_lin(x, y, "StaticEnergy" );
	Entropy           = Interp2D_lin(x, y, "Entropy" );
	Pressure          = Interp2D_lin(x, y, "Pressure" );
	SoundSpeed2       = Interp2D_lin(x, y, "SoundSpeed2" );
	dPdrho_e          = Interp2D_lin(x, y, "dPdrho_e" );
	dPde_rho          = Interp2D_lin(x, y, "dPde_rho" );
	dTdrho_e          = Interp2D_lin(x, y, "dTdrho_e" );
	dTde_rho          = Interp2D_lin(x, y, "dTde_rho" );
	Cp                = Interp2D_lin(x, y, "Cp" );
	Mu                = Interp2D_lin(x, y, "Mu" );
	dmudrho_T         = Interp2D_lin(x, y, "dmudrho_T" );
	dmudT_rho         = Interp2D_lin(x, y, "dmudT_rho" );
	Kt                = Interp2D_lin(x, y, "Kt" );
	dktdrho_T         = Interp2D_lin(x, y, "dktdrho_T" );
	dktdT_rho         = Interp2D_lin(x, y, "dktdT_rho" );
	//Intermediate variables only needed for StandAlone version
	su2double Density = Density;
	su2double Pressure = Pressure;
	cout<<"Interpolated fit:"<<endl;
	if ((Density>Density_limits[1]) or (Density<Density_limits[0]))
	{
		cerr<<"RHOT Interpolated Density out of bounds\n";
	}
	if ((Pressure>Pressure_limits[1]) or (Pressure<Pressure_limits[0]))
	{
		cerr<<"RHOT Interpolated Pressure out of bounds\n";
	}
}


void CLookUpTable::Interp2D_ArbitrarySkewCoeff(su2double x, su2double y, std::string grid_var)
{
	//Distances in along x and y axis are taken relative to i,j point (x00, y00).
	//This reduces the interpolation to a 3by3 system rather than 4by4
	//x and y are not strictrly necessary for the calculation of the coefficients.
	//However, they do allow for checking whether the point of interest is contained in
	//the quad under consideration.
	su2double x00, y00, dx10, dx01, dx11, dy10, dy01, dy11;
	//Interpolation LHM
	su2double A[3][3];
	//Helper variable for Gaussian elimination
	su2double c;
	//Load in the coordinates of the qudrilateral (values relative to i,j)
	if(grid_var=="RHOE")
	{
		x00  = ThermoTables[NN_i[0]][NN_j[0]].Density     ;
		y00  = ThermoTables[NN_i[0]][NN_j[0]].StaticEnergy;
		dx01 = ThermoTables[NN_i[2]][NN_j[2]].Density      -x00;
		dy01 = ThermoTables[NN_i[2]][NN_j[2]].StaticEnergy -y00;
		dx10 = ThermoTables[NN_i[1]][NN_j[1]].Density      -x00;
		dy10 = ThermoTables[NN_i[1]][NN_j[1]].StaticEnergy -y00;
		dx11 = ThermoTables[NN_i[3]][NN_j[3]].Density      -x00;
		dy11 = ThermoTables[NN_i[3]][NN_j[3]].StaticEnergy -y00;
	}
	else if(grid_var=="PT")
	{
		y00  = ThermoTables[NN_i[0]][NN_j[0]].Pressure   ;
		x00  = ThermoTables[NN_i[0]][NN_j[0]].Temperature;
		dy01 = ThermoTables[NN_i[2]][NN_j[2]].Pressure    -y00;
		dx01 = ThermoTables[NN_i[2]][NN_j[2]].Temperature -x00;
		dy10 = ThermoTables[NN_i[1]][NN_j[1]].Pressure    -y00;
		dx10 = ThermoTables[NN_i[1]][NN_j[1]].Temperature -x00;
		dy11 = ThermoTables[NN_i[3]][NN_j[3]].Pressure    -y00;
		dx11 = ThermoTables[NN_i[3]][NN_j[3]].Temperature -x00;
	}
	else if(grid_var=="PRHO")
	{
		y00  = ThermoTables[NN_i[0]][NN_j[0]].Pressure;
		x00  = ThermoTables[NN_i[0]][NN_j[0]].Density ;
		dy01 = ThermoTables[NN_i[2]][NN_j[2]].Pressure -y00;
		dx01 = ThermoTables[NN_i[2]][NN_j[2]].Density  -x00;
		dy10 = ThermoTables[NN_i[1]][NN_j[1]].Pressure -y00;
		dx10 = ThermoTables[NN_i[1]][NN_j[1]].Density  -x00;
		dy11 = ThermoTables[NN_i[3]][NN_j[3]].Pressure -y00;
		dx11 = ThermoTables[NN_i[3]][NN_j[3]].Density  -x00;
	}
	else if(grid_var=="RHOT")
	{
		x00  = ThermoTables[NN_i[0]][NN_j[0]].Density    ;
		y00  = ThermoTables[NN_i[0]][NN_j[0]].Temperature;
		dx01 = ThermoTables[NN_i[2]][NN_j[2]].Density     -x00;
		dy01 = ThermoTables[NN_i[2]][NN_j[2]].Temperature -y00;
		dx10 = ThermoTables[NN_i[1]][NN_j[1]].Density     -x00;
		dy10 = ThermoTables[NN_i[1]][NN_j[1]].Temperature -y00;
		dx11 = ThermoTables[NN_i[3]][NN_j[3]].Density     -x00;
		dy11 = ThermoTables[NN_i[3]][NN_j[3]].Temperature -y00;
	}
	else if(grid_var=="PS")
	{
		y00  = ThermoTables[NN_i[0]][NN_j[0]].Pressure;
		x00  = ThermoTables[NN_i[0]][NN_j[0]].Entropy ;
		dy01 = ThermoTables[NN_i[2]][NN_j[2]].Pressure -y00;
		dx01 = ThermoTables[NN_i[2]][NN_j[2]].Entropy  -x00;
		dy10 = ThermoTables[NN_i[1]][NN_j[1]].Pressure -y00;
		dx10 = ThermoTables[NN_i[1]][NN_j[1]].Entropy  -x00;
		dy11 = ThermoTables[NN_i[3]][NN_j[3]].Pressure -y00;
		dx11 = ThermoTables[NN_i[3]][NN_j[3]].Entropy  -x00;
	}
	else if(grid_var=="HS")
	{
		x00  = ThermoTables[NN_i[0]][NN_j[0]].Enthalpy;
		y00  = ThermoTables[NN_i[0]][NN_j[0]].Entropy     ;
		dx01 = ThermoTables[NN_i[2]][NN_j[2]].Enthalpy  -x00;
		dy01 = ThermoTables[NN_i[2]][NN_j[2]].Entropy   -y00;
		dx10 = ThermoTables[NN_i[1]][NN_j[1]].Enthalpy  -x00;
		dy10 = ThermoTables[NN_i[1]][NN_j[1]].Entropy   -y00;
		dx11 = ThermoTables[NN_i[3]][NN_j[3]].Enthalpy 	-x00;
		dy11 = ThermoTables[NN_i[3]][NN_j[3]].Entropy   -y00;
	}
	//Check if x, y is indeed in the quad
	//Some extra logic is needed as the both monotonically increasing and monotonically decreasing functions
	//have to be anticipated
	bool BOTTOM, TOP, LEFT, RIGHT;
	BOTTOM = (y*dx10)<(x*dy10);
	TOP = ((y-dy01)*(dx11-dx01))>((dy11-dy01)*(x-dx01));
	RIGHT = ((x-dx10)*(dy11-dy10))>((dx11-dx10)*(y-dy10));
	LEFT = (x*dy01)<(dx01*y);
	//Check BOTTOM quad boundary
	if(BOTTOM and !TOP)
	{
		//added table limit detection
		if (jIndex==0)
		{
			cerr<<grid_var+" interpolation point lies below the LUT\n";
		}
		else
		{
			cerr<<grid_var+" interpolation point lies below bottom boundary of selected quad\n";
		}
	}
	//Check RIGHT quad boundary
	if(RIGHT and !LEFT)
	{
		//added table limit detection
		if (iIndex==(rho_dim-1))
		{
			cerr<<grid_var+" interpolation point lies right of the LUT\n";
		}
		else
		{
			cerr<<grid_var+" interpolation point lies to the right of the boundary of selected quad\n";
		}
	}
	//Check TOP quad boundary
	if(TOP and !BOTTOM)
	{
		//added table limit detection
		if (jIndex==p_dim-1)
		{
			cerr<<grid_var+" interpolation point lies above the LUT\n";
		}
		else
		{
			cerr<<grid_var+" interpolation point lies above the boundary of selected quad\n";
		}
	}
	//Check LEFT quad boundary

	if(LEFT and !RIGHT)
	{
		//added table limit detection
		if (iIndex==0)
		{
			cerr<<grid_var+" interpolation point lies left of the LUT\n";
		}
		else
		{
			cerr<<grid_var+" interpolation point lies to the left of the boundary of selected quad\n";
		}
	}
	//Setup the LHM matrix for the interpolation
	A[0][0] = dx10;
	A[0][1] = dy10;
	A[0][2] = dx10*dy10;
	A[1][0] = dx01;
	A[1][1] = dy01;
	A[1][2] = dx01*dy01;
	A[2][0] = dx11;
	A[2][1] = dy11;
	A[2][2] = dx11*dy11;
	cout<<"Interpolation LHM matrix \n"<<"[";

//	for (int j=0; j<3; j++)
//	{
//		cout<<setw(15)<<"["<<A[j][0]<<" ,  "<<A[j][1]<<"  , "<<A[j][2]<<"]"<<endl;
//	}
//	cout<<"]\n";

	//Store the inverse of the LHM matrix as coeff
	coeff[0][0] = 1;
	coeff[0][1] = 0;
	coeff[0][2] = 0;
	coeff[1][0] = 0;
	coeff[1][1] = 1;
	coeff[1][2] = 0;//solved interpolation bug
	coeff[2][0] = 0;
	coeff[2][1] = 0;
	coeff[2][2] = 1;

	//Compute inverse of LHM using Gaussian elimination
	//Reduced Echelon form of the LHM
	if (A[0][0] != 0)
	{
		c = A[1][0]/A[0][0];
		coeff[1][0] = coeff[1][0] -coeff[0][0]*c;
		for (int i=0; i<3; i++)
		{
			A[1][i] = A[1][i] -A[0][i]*c;
		}
		c =A[2][0]/A[0][0];
		coeff[2][0] = coeff[2][0] -coeff[0][0]*c;

		for (int i=0; i<3; i++)
		{
			A[2][i] = A[2][i] -A[0][i]*c;
		}
	}

	if (A[1][1] != 0)
	{
		c = A[2][1]/A[1][1];
		for (int i=0; i<2; i++)
		{
			coeff[2][i] = coeff[2][i] -coeff[1][i]*c;
		}
		for (int i=0; i<3; i++)
			A[2][i] = A[2][i] -A[1][i]*c;
	}
	//Reduced reduced Echelon form of LHM
	if (A[2][2] != 0)
	{

		for (int i=0; i<3; i++)
		{
			coeff[1][i] = coeff[1][i] -coeff[2][i]*A[1][2]/A[2][2];
		}
		A[1][2] = 0;
		for (int i=0; i<3; i++)
		{
			coeff[0][i] = coeff[0][i] -coeff[2][i]*A[0][2]/A[2][2];
		}
		A[0][2] = 0;
	}
	if (A[1][1] != 0)
	{
		for (int i=0; i<3; i++)
			coeff[0][i] = coeff[0][i] -coeff[1][i]*A[0][1]/A[1][1];
	}
	A[0][1] = 0;

	//Normalize the RR Echelon form
	if (A[0][0] != 0)
	{
		for (int i=0; i<3; i++)
		{
			coeff[0][i] = coeff[0][i]/A[0][0];
		}
		if (A[1][1] != 0)
		{
			for (int i=0; i<3; i++)
			{
				coeff[1][i] = coeff[1][i]/A[1][1];
			}
		}
		if (A[2][2] != 0)
		{
			for (int i=0; i<3; i++)
			{
				coeff[2][i] = coeff[2][i]/A[2][2];
			}

		}
	}
	return;
}
su2double CLookUpTable::Interp2D_Inv_Dist(std::string interpolant_var, su2double* dist)
{
	su2double interp_result=0;
	//The function values to interpolate from
	su2double F[4];
	//For each case the values are filled differently
	if(interpolant_var=="StaticEnergy")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].StaticEnergy;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].StaticEnergy;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].StaticEnergy;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].StaticEnergy;
	}
	else if(interpolant_var=="Entropy")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].Entropy;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].Entropy;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].Entropy;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].Entropy;
	}
	else if(interpolant_var=="Density")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].Density;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].Density;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].Density;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].Density;
	}
	else if(interpolant_var=="Pressure")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].Pressure;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].Pressure;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].Pressure;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].Pressure;
	}
	else if(interpolant_var=="SoundSpeed2")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].SoundSpeed2;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].SoundSpeed2;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].SoundSpeed2;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].SoundSpeed2;
	}
	else if(interpolant_var=="Temperature")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].Temperature;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].Temperature;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].Temperature;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].Temperature;
	}
	else if(interpolant_var=="dPdrho_e")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].dPdrho_e;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].dPdrho_e;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].dPdrho_e;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].dPdrho_e;
	}
	else if(interpolant_var=="dPde_rho")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].dPde_rho;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].dPde_rho;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].dPde_rho;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].dPde_rho;
	}
	else if(interpolant_var=="dTdrho_e")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].dTdrho_e;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].dTdrho_e;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].dTdrho_e;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].dTdrho_e;
	}
	else if(interpolant_var=="dTde_rho")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].dTde_rho;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].dTde_rho;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].dTde_rho;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].dTde_rho;
	}
	else if(interpolant_var=="Cp")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].Cp;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].Cp;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].Cp;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].Cp;
	}
	else if(interpolant_var=="Mu")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].Mu;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].Mu;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].Mu;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].Mu;
	}
	else if(interpolant_var=="dmudrho_T")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].dmudrho_T;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].dmudrho_T;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].dmudrho_T;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].dmudrho_T;
	}
	else if(interpolant_var=="dmudT_rho")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].dmudT_rho;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].dmudT_rho;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].dmudT_rho;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].dmudT_rho;
	}
	else if(interpolant_var=="Kt")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].Kt;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].Kt;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].Kt;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].Kt;
	}
	else if(interpolant_var=="dktdrho_T")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].dktdrho_T;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].dktdrho_T;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].dktdrho_T;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].dktdrho_T;
	}
	else if(interpolant_var=="dktdT_rho")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].dktdT_rho;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].dktdT_rho;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].dktdT_rho;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].dktdT_rho;
	}
	else if(interpolant_var=="Enthalpy")
	{
		F[0] = ThermoTables[NN_i[0]][NN_j[0]].Enthalpy;
		F[1] = ThermoTables[NN_i[1]][NN_j[1]].Enthalpy;
		F[2] = ThermoTables[NN_i[2]][NN_j[2]].Enthalpy;
		F[3] = ThermoTables[NN_i[3]][NN_j[3]].Enthalpy;
	}

	su2double dist_sum = 0;
	for (int i=0; i<4;i++)
	{
		interp_result += (1/dist[i])*F[i];
		dist_sum += 1/dist[i];
	}

	interp_result = interp_result/dist_sum;

	return interp_result;
}


su2double CLookUpTable::Interp2D_lin(su2double x, su2double y, string interpolant_var)
{
	//F is the RHS part of the interpolation equation
	su2double F[3];
	//The solution vector for the interpolation equation
	su2double C[3];
	//The function values
	su2double f00, f10, f01, f11;
	//For each case the values are filled differently
	if(interpolant_var=="StaticEnergy")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].StaticEnergy;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].StaticEnergy;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].StaticEnergy;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].StaticEnergy;
	}
	else if(interpolant_var=="Entropy")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].Entropy;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].Entropy;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].Entropy;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].Entropy;
	}
	else if(interpolant_var=="Density")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].Density;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].Density;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].Density;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].Density;
	}
	else if(interpolant_var=="Pressure")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].Pressure;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].Pressure;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].Pressure;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].Pressure;
	}
	else if(interpolant_var=="SoundSpeed2")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].SoundSpeed2;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].SoundSpeed2;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].SoundSpeed2;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].SoundSpeed2;
	}
	else if(interpolant_var=="Temperature")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].Temperature;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].Temperature;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].Temperature;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].Temperature;
	}
	else if(interpolant_var=="dPdrho_e")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].dPdrho_e;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].dPdrho_e;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].dPdrho_e;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].dPdrho_e;
	}
	else if(interpolant_var=="dPde_rho")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].dPde_rho;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].dPde_rho;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].dPde_rho;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].dPde_rho;
	}
	else if(interpolant_var=="dTdrho_e")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].dTdrho_e;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].dTdrho_e;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].dTdrho_e;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].dTdrho_e;
	}
	else if(interpolant_var=="dTde_rho")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].dTde_rho;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].dTde_rho;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].dTde_rho;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].dTde_rho;
	}
	else if(interpolant_var=="Cp")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].Cp;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].Cp;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].Cp;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].Cp;
	}
	else if(interpolant_var=="Mu")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].Mu;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].Mu;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].Mu;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].Mu;
	}
	else if(interpolant_var=="dmudrho_T")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].dmudrho_T;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].dmudrho_T;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].dmudrho_T;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].dmudrho_T;
	}
	else if(interpolant_var=="dmudT_rho")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].dmudT_rho;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].dmudT_rho;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].dmudT_rho;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].dmudT_rho;
	}
	else if(interpolant_var=="Kt")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].Kt;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].Kt;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].Kt;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].Kt;
	}
	else if(interpolant_var=="dktdrho_T")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].dktdrho_T;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].dktdrho_T;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].dktdrho_T;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].dktdrho_T;
	}
	else if(interpolant_var=="dktdT_rho")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].dktdT_rho;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].dktdT_rho;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].dktdT_rho;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].dktdT_rho;
	}
	else if(interpolant_var=="Enthalpy")
	{
		f00 = ThermoTables[NN_i[0]][NN_j[0]].Enthalpy;
		f10 = ThermoTables[NN_i[1]][NN_j[1]].Enthalpy;
		f01 = ThermoTables[NN_i[2]][NN_j[2]].Enthalpy;
		f11 = ThermoTables[NN_i[3]][NN_j[3]].Enthalpy;
	}

	//Using offset relative to i,j point yields a 3by3 system rather than 4by4
	F[0] = f10 - f00;
	C[0] = 0;
	F[1] = f01 - f00;
	C[1] = 0;
	F[2] = f11 - f00;
	C[2] = 0;
	for (int i = 0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			C[i] = C[i] + F[i]*coeff[i][j];
		}
	}

	return f00 + C[0]*x + C[1]*y + C[2]*x*y;
}


void CLookUpTable::LUTprint(void)
{
//	for (int i=0; i<rho_dim; i++)
//	{
//		for (int j=0; j<p_dim; j++)
//		{
////			ThermoTables[i][j].CTLprint();
//		}
//	}
}

void CLookUpTable::RecordState(char* file)
{
//	fstream fs;
//	fs.open(file, fstream::app);
//	assert(fs.is_open());
//	fs << Temperature<<", ";
//	fs << Density<<", ";
//	fs << Enthalpy<<", ";
//	fs << StaticEnergy<<", ";
//	fs << Entropy<<", ";
//	fs << Pressure<<", ";
//	fs << SoundSpeed2<<", ";
//	fs << dPdrho_e<<", ";
//	fs << dPde_rho<<", ";
//	fs << dTdrho_e<<", ";
//	fs << dTde_rho<<", ";
//	fs << Cp<<", ";
//	fs << Mu<<", ";
//	fs << dmudrho_T<<", ";
//	fs << dmudT_rho<<", ";
//	fs << Kt<<", ";
//	fs << dktdrho_T<<", ";
//	fs << dktdT_rho<<", ";
//	fs << "\n";
//	fs.close();

}

void CLookUpTable::TableDump(char* filename)
{
//	for (int i=0; i<rho_dim;i++)
//	{
//		for (int j=0; j<p_dim;j++)
//		{
//			iIndex = i;
//			jIndex = j;
//			Temperature       = ThermoTables[iIndex][jIndex].Temperature;
//			Density           = ThermoTables[iIndex][jIndex].Density;
//			Enthalpy          = ThermoTables[iIndex][jIndex].Enthalpy;
//			StaticEnergy      = ThermoTables[iIndex][jIndex].StaticEnergy;
//			Entropy           = ThermoTables[iIndex][jIndex].Entropy;
//			Pressure          = ThermoTables[iIndex][jIndex].Pressure;
//			SoundSpeed2       = ThermoTables[iIndex][jIndex].SoundSpeed2;
//			dPdrho_e          = ThermoTables[iIndex][jIndex].dPdrho_e;
//			dPde_rho          = ThermoTables[iIndex][jIndex].dPde_rho;
//			dTdrho_e          = ThermoTables[iIndex][jIndex].dTdrho_e;
//			dTde_rho          = ThermoTables[iIndex][jIndex].dTde_rho;
//			Cp                = ThermoTables[iIndex][jIndex].Cp;
//			Mu                = ThermoTables[iIndex][jIndex].Mu;
//			dmudrho_T         = ThermoTables[iIndex][jIndex].dmudrho_T;
//			dmudT_rho         = ThermoTables[iIndex][jIndex].dmudT_rho;
//			Kt                = ThermoTables[iIndex][jIndex].Kt;
//			dktdrho_T         = ThermoTables[iIndex][jIndex].dktdrho_T;
//			dktdT_rho         = ThermoTables[iIndex][jIndex].dktdT_rho;
//			RecordState(filename);
//		}
//	}

}

void CLookUpTable::TableLoadCFX(string filename){
	int N_PARAM = 0;
	int set_x = 0;
	int set_y = 0;
	int var_steps = 0;
	int var_scanned=0;

	string line;
	string value;


	ifstream table (filename.c_str());
	assert(table.is_open());
	cout<<"Looking for number of parameters"<<endl;
	while ( getline(table,line) )
	{
		unsigned int found;
		found = line.find("$$PARAM");
		if (found<10)
		{
			getline(table,line);
			istringstream in(line);
			in>>N_PARAM;
			N_PARAM++;
			cout<<"Number of parameters "<<N_PARAM<<endl;
		}
		for (int var=var_scanned; var<N_PARAM+1; var++)
		{
			string svar = static_cast<ostringstream*>( &(ostringstream() << var) )->str();
			found = line.find("$TABLE_"+svar);
			if (found<10)
			{
				var_scanned = var;
				cout<<found<<' '<<line<<endl;
				getline(table,line);
				istringstream in(line);
				int x,y;
				in>>x>>y;
				if (var==1)
				{
					ThermoTables = new CThermoList*[x];
					for (int i=0; i<x; i++)
					{
						ThermoTables[i] = new CThermoList[y];
					}
					set_x = x;
					set_y = y;
					cout<<"Tables have been allocated"<<var<<endl;

					//Fill in the densities
					Density_limits[0] = 10E15;
					Density_limits[1] = 0;
					var_steps = 10;
					su2double* vD = new su2double[set_x];

					for (int k =0; k<ceil(float(set_x)/10.0);k++)
					{
						getline(table,line);
						cout<<line<<endl;
						istringstream inD(line);
						if ((set_x-k*10)<10) var_steps = (set_x-k*10);
						for (int i =0; i<var_steps; i++)
						{
							inD>>vD[10*k+i];
							if (vD[10*k+i]>Density_limits[1])
							{
								Density_limits[1]=vD[10*k+i];
							}
							if (vD[10*k+i]<Density_limits[0])
							{
								Density_limits[0]=vD[10*k+i];
							}
						}
					}
					for(int i =0; i<set_x; i++)
					{
						for (int j =0; j<set_y; j++)
						{
							ThermoTables[i][j].Density = vD[i];
						}
					}
					delete vD;


					//Fill in the pressures
					su2double* vP = new su2double[set_y];
					var_steps = 10; //solved pressure reading bug
					Pressure_limits[0] = 10E15; //lower limit
					Pressure_limits[1] = 0; //upper limit
					//Each line contains at most 10 pressure values
					for (int k =0; k<ceil(float(set_y)/10.0);k++)
					{

						getline(table,line);
						cout<<line<<endl;
						istringstream inP(line);
						//Check if line contains less than 10 values
						if ((set_y-k*10)<10) var_steps = (set_y-k*10);
						for (int j =0; j<var_steps; j++)
						{
							inP>>vP[10*k+j];
							if (vP[10*k+j]>Pressure_limits[1])
							{
								Pressure_limits[1]=vP[10*k+j];
							}
							if (vP[10*k+j]<Pressure_limits[0])
							{
								Pressure_limits[0]=vP[10*k+j];
							}
						}
					}
					for(int i =0; i<set_x; i++)
					{
						for (int j =0; j<set_y; j++)
						{
							ThermoTables[i][j].Pressure = vP[j];
						}
					}
					delete vP;
					cout<<"Tables have been filled with D and P values "<<var<<endl;

				}
				// Check that additional tables all adhere to the same x,y dimensions, otherwise throw an error
				else if (x != set_x && y!=set_y)
				{
					cerr<<"The encountered dimensions of the CFX table are not the same throughout. They should be; for this to work.\n";

				}
				//Go through each one of the variables of interest
				if(var==16)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					StaticEnergy_limits[0] = 10E20;//lower limit
					StaticEnergy_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].StaticEnergy = inp[i%10];
							if (inp[i%10]>StaticEnergy_limits[1])
							{
								StaticEnergy_limits[1]= inp[i%10];
							}
							if (inp[i%10]<StaticEnergy_limits[0])
							{
								StaticEnergy_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with Static Energy values "<<var<<endl;
				}
				if(var==1)
				{
					//Fixed a bug: lines already skipped for var==1
					//for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					//for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure
					cout<<"set_x "<<set_x<<endl;
					cout<<"set_y "<<set_y<<endl;
					Enthalpy_limits[0] = 10E20;//lower limit
					Enthalpy_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Enthalpy = inp[i%10];
							if (inp[i%10]>Enthalpy_limits[1])
							{
								Enthalpy_limits[1]= inp[i%10];
							}
							if (inp[i%10]<Enthalpy_limits[0])
							{
								Enthalpy_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with speed of specific enthalpy values "<<var<<endl;
				}
				if(var==2)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					SoundSpeed2_limits[0] = 10E20;//lower limit
					SoundSpeed2_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].SoundSpeed2 = inp[i%10];
							if (inp[i%10]>SoundSpeed2_limits[1])
							{
								SoundSpeed2_limits[1]= inp[i%10];
							}
							if (inp[i%10]<SoundSpeed2_limits[0])
							{
								SoundSpeed2_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with speed of sound values "<<var<<endl;
				}
				if(var==5)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					Cp_limits[0] = 10E20;//lower limit
					Cp_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Cp = inp[i%10];
							if (inp[i%10]>Cp_limits[1])
							{
								Cp_limits[1]= inp[i%10];
							}
							if (inp[i%10]<Cp_limits[0])
							{
								Cp_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with isobaric heat capacity values "<<var<<endl;
				}
				if(var==7)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					Entropy_limits[0] = 10E20;//lower limit
					Entropy_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Entropy = inp[i%10];
							if (inp[i%10]>Entropy_limits[1])
							{
								Entropy_limits[1]= inp[i%10];
							}
							if (inp[i%10]<Entropy_limits[0])
							{
								Entropy_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with entropy values "<<var<<endl;

				}
				if(var==8)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					Mu_limits[0] = 10E20;//lower limit
					Mu_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Mu = inp[i%10];
							if (inp[i%10]>Mu_limits[1])
							{
								Mu_limits[1]= inp[i%10];
							}
							if (inp[i%10]<Mu_limits[0])
							{
								Mu_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with viscosity values "<<var<<endl;

				}
				if(var==9)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					Kt_limits[0] = 10E20;//lower limit
					Kt_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Kt = inp[i%10];
							if (inp[i%10]>Kt_limits[1])
							{
								Kt_limits[1]= inp[i%10];
							}
							if (inp[i%10]<Kt_limits[0])
							{
								Kt_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with thermal conductivity values "<<var<<endl;

				}
				if(var==10)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					dPdrho_e_limits[0] = 10E20;//lower limit
					dPdrho_e_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].dPdrho_e = inp[i%10];
							if (inp[i%10]>dPdrho_e_limits[1])
							{
								dPdrho_e_limits[1]= inp[i%10];
							}
							if (inp[i%10]<dPdrho_e_limits[0])
							{
								dPdrho_e_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with specific dPdrho_e values "<<var<<endl;

				}
				if(var==11)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					dPde_rho_limits[0] = 10E20;//lower limit
					dPde_rho_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].dPde_rho = inp[i%10];
							if (inp[i%10]>dPde_rho_limits[1])
							{
								dPde_rho_limits[1]= inp[i%10];
							}
							if (inp[i%10]<dPde_rho_limits[0])
							{
								dPde_rho_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with specific dPde_rho values "<<var<<endl;
				}
				if(var==12)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure

					dTdrho_e_limits[0] = 10E20;//lower limit
					dTdrho_e_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].dTdrho_e = inp[i%10];
							if (inp[i%10]>dTdrho_e_limits[1])
							{
								dTdrho_e_limits[1]= inp[i%10];
							}
							if (inp[i%10]<dTdrho_e_limits[0])
							{
								dTdrho_e_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with specific dTdrho_e values "<<var<<endl;
				}
				if(var==13)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					dTde_rho_limits[0] = 10E20;//lower limit
					dTde_rho_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].dTde_rho = inp[i%10];
							if (inp[i%10]>dTde_rho_limits[1])
							{
								dTde_rho_limits[1]= inp[i%10];
							}
							if (inp[i%10]<dTde_rho_limits[0])
							{
								dTde_rho_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with specific dTde_rho values "<<var<<endl;
				}
				if(var==15)
				{
					for (int k =0; k<ceil(float(set_x)/10.0);k++) getline(table,line); //skip density
					for (int k =0; k<ceil(float(set_y)/10.0);k++) getline(table,line); //skip pressure


					Temperature_limits[0] = 10E20;//lower limit
					Temperature_limits[1] = 0;//upper limit

					su2double inp[10];

					for (int j =0; j<set_y; j++)
					{
						for(int i =0; i<set_x; i++)
						{
							if ((j*set_x+i)%10==0)
							{
								getline(table,line);
								cout<<line<<endl;
								istringstream in(line);
								for (int z = 0; z<10; z++)
								{
									in>>inp[z];
								}
							}
							ThermoTables[i][j].Temperature = inp[i%10];
							if (inp[i%10]>Temperature_limits[1])
							{
								Temperature_limits[1]= inp[i%10];
							}
							if (inp[i%10]<Temperature_limits[0])
							{
								Temperature_limits[0]=inp[i%10];
							}
						}

					}
					cout<<"Tables have been filled with Temperature values "<<var<<endl;
				}

			}
		}
	}
	rho_dim = set_x;
	p_dim = set_y;
	table.close();
}










