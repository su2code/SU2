/*!
 * \file numerics_structure.hpp
 * \brief Headers of the main subroutines for the dumerical definition of the problem.
 *        The subroutines and functions are in the <i>numerics_structure.cpp</i>,
 *        <i>numerics_convective.cpp</i>, <i>numerics_viscous.cpp</i>, and
 *        <i>numerics_source.cpp</i> files.
 * \author B. Tracey
 * \version 4.0.0 "Cardinal"
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

#pragma once

#include "../../Common/include/mpi_structure.hpp"

#ifdef HAVE_JSONCPP
  #include <json/json.h>
#endif
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <math.h>
#include <sstream>

#include "../include/numerics_machine_learning_turbulent.hpp"

using namespace std;

class CScaler{
public:
	CScaler();
  virtual ~CScaler();
	virtual void Scale(su2double *) = 0;
	virtual void Unscale(su2double *) = 0;
};

class CNormalScaler: public CScaler{
private:
	su2double * mu;
	su2double * sigma;
	int dim;

public:
	CNormalScaler();
	CNormalScaler(int, su2double*, su2double*);
#ifdef HAVE_JSONCPP
  CNormalScaler(Json::Value);
#endif
	~CNormalScaler();
	void Scale(su2double *);
	void Unscale(su2double *);
};

class CMulInputScaler : public CScaler{
public:
  su2double MulScale;
  CScaler* InnerScaler;
  CMulInputScaler();
#ifdef HAVE_JSONCPP
  CMulInputScaler(Json::Value);
#endif
  ~CMulInputScaler();
  void Scale(su2double *);
	void Unscale(su2double *);
};

class CMulOutputScaler : public CScaler{
public:
  su2double MulScale;
  CMulOutputScaler();
#ifdef HAVE_JSONCPP
  CMulOutputScaler(Json::Value);
#endif
  ~CMulOutputScaler();
  void Scale(su2double *);
	void Unscale(su2double *);
  
};

class CActivator{
public:
	CActivator();
	virtual ~CActivator();
	virtual su2double Activate(su2double combination) {cout<< "IN BASE ACTIVATOR THIS IS BAD" << endl; return 0;};
};

class CTanhActivator : public CActivator{
public:
	CTanhActivator();
#ifdef HAVE_JSONCPP
  CTanhActivator(Json::Value);
#endif
	~CTanhActivator();
	su2double Activate(su2double combination);
};

class CLinearActivator : public CActivator{
public:
	CLinearActivator();
#ifdef HAVE_JSONCPP
  CLinearActivator(Json::Value);
#endif
	~CLinearActivator();
	su2double Activate(su2double combination);
};

class CNeuron{
public:
  CNeuron();
  ~CNeuron();
  virtual su2double Activate(su2double combination) {cout << "In base neuron. Bad";return 0;};
  virtual su2double Combine(su2double * parameters, int nParameters, su2double *inputs, int nInputs) {cout << "In base neuron. Bad";return 0;};
};

class CSumNeuron : public CNeuron{
private:
		CActivator *activator;
public:
	CSumNeuron();
	CSumNeuron(CActivator*); // activator, parameterStart, nParameters
#ifdef HAVE_JSONCPP
  CSumNeuron(Json::Value);
#endif
	~CSumNeuron();
	//Activator* GetActivator(void);
	su2double Combine(su2double * parameters, int nParameters, su2double * inputs, int nInputs);
  su2double Activate(su2double combination);
};


class CPredictor{
protected:
  int inputDim;
  int outputDim;
public:
  CPredictor();
  virtual  ~CPredictor();
  virtual void Predict(su2double *, su2double *) {cout << "In base Predict, this is bad";};
  int InputDim();
  int OutputDim();
};

class CScalePredictor{
public:
  CScalePredictor();
  CScalePredictor(string filename);
  ~CScalePredictor();
  CPredictor *Pred;
  CScaler *InputScaler;
  CScaler *OutputScaler;
public:
  void Predict(su2double *inputs, su2double *outputs);
  // Need to add predict method
};

class CMulPredictor : public CPredictor{
public:
  CMulPredictor();
#ifdef HAVE_JSONCPP
  CMulPredictor(Json::Value);
#endif
  ~CMulPredictor();
  CPredictor* Inner;
  void Predict(su2double *, su2double *);
  
};

class CNeurNet : public CPredictor {
private:
  int maxNeurons; // Number of neurons in the largest layer
  int nLayers;
  CNeuron ***neurons; // Array of arrays to pointers to neuron
  int* nNeuronsInLayer; //one list for each layer
  int** nParameters; // Number of parameters for the neuron
//  int inputDim;
  
  void processLayer(su2double *, int, CNeuron **, su2double **, int, int * , su2double *);
  
  //----
//	int outputDim;
public:
	CNeurNet();
#ifdef HAVE_JSONCPP
  CNeurNet(Json::Value);
#endif
	~CNeurNet();
//	int InputDim();
//	int OutputDim();
	void Predict(su2double *, su2double *);
    su2double*** parameters; // Array of parameters for each neuron
};

class CSANondimInputs{
private:
  int nDim;
public:
  CSANondimInputs(int);
  ~CSANondimInputs();
  void Set(SpalartAllmarasInputs*);
  void NondimensionalizeSource(int, su2double*);
  void DimensionalizeSource(int, su2double*);
  su2double Chi;
  su2double OmegaNondim;
  su2double OmegaBar;
  su2double SourceNondim;
  su2double NuGradNondim;
  su2double * DNuHatDXBar;
  su2double NuHatGradNorm;
  su2double NuHatGradNormBar;
};

#include "numerics_machine_learning.inl"
