/*!
 * \file numerics_structure.hpp
 * \brief Headers of the main subroutines for the dumerical definition of the problem.
 *        The subroutines and functions are in the <i>numerics_structure.cpp</i>,
 *        <i>numerics_convective.cpp</i>, <i>numerics_viscous.cpp</i>, and
 *        <i>numerics_source.cpp</i> files.
 * \author B. Tracey
 * \version 3.2.5 "eagle"
 *
 * Copyright (C) 2012-2014 SU2 <https://github.com/su2code>.
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

#ifdef HAVE_MPI
  #include "mpi.h"
#endif
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
	virtual void Scale(double *) = 0;
	virtual void Unscale(double *) = 0;
};

class CNormalScaler: public CScaler{
private:
	double * mu;
	double * sigma;
	int dim;

public:
	CNormalScaler();
	CNormalScaler(int,double*,double*);
#ifdef HAVE_JSONCPP
  CNormalScaler(Json::Value);
#endif
	~CNormalScaler();
	void Scale(double *);
	void Unscale(double *);
};

class CMulInputScaler : public CScaler{
public:
  double MulScale;
  CScaler* InnerScaler;
  CMulInputScaler();
#ifdef HAVE_JSONCPP
  CMulInputScaler(Json::Value);
#endif
  ~CMulInputScaler();
  void Scale(double *);
	void Unscale(double *);
};

class CMulOutputScaler : public CScaler{
public:
  double MulScale;
  CMulOutputScaler();
#ifdef HAVE_JSONCPP
  CMulOutputScaler(Json::Value);
#endif
  ~CMulOutputScaler();
  void Scale(double *);
	void Unscale(double *);
  
};

class CActivator{
public:
	CActivator();
	~CActivator();
	virtual double Activate(double combination){cout<< "IN BASE ACTIVATOR THIS IS BAD" <<endl; return 0;};
};

class CTanhActivator : public CActivator{
public:
	CTanhActivator();
#ifdef HAVE_JSONCPP
  CTanhActivator(Json::Value);
#endif
	~CTanhActivator();
	double Activate(double combination);
};

class CLinearActivator : public CActivator{
public:
	CLinearActivator();
#ifdef HAVE_JSONCPP
  CLinearActivator(Json::Value);
#endif
	~CLinearActivator();
	double Activate(double combination);
};

class CNeuron{
public:
  CNeuron();
  ~CNeuron();
  virtual double Activate(double combination){cout << "In base neuron. Bad";return 0;};
  virtual double Combine(double * parameters, int nParameters, double *inputs, int nInputs){cout << "In base neuron. Bad";return 0;};
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
	double Combine(double * parameters, int nParameters, double * inputs, int nInputs);
  double Activate(double combination);
};


class CPredictor{
protected:
  int inputDim;
  int outputDim;
public:
  CPredictor();
  virtual  ~CPredictor();
  virtual void Predict(double *, double *){cout << "In base Predict, this is bad";};
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
  void Predict(double *inputs, double *outputs);
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
  void Predict(double *, double *);
  
};

class CNeurNet : public CPredictor {
private:
  int maxNeurons; // Number of neurons in the largest layer
  int nLayers;
  CNeuron ***neurons; // Array of arrays to pointers to neuron
  int* nNeuronsInLayer; //one list for each layer
  int** nParameters; // Number of parameters for the neuron
//  int inputDim;
  
  void processLayer(double *, int,CNeuron **, double **, int, int * ,double *);
  
  //----
//	int outputDim;
	int totalNumParameters;
public:
	CNeurNet();
#ifdef HAVE_JSONCPP
  CNeurNet(Json::Value);
#endif
	~CNeurNet();
//	int InputDim();
//	int OutputDim();
	void Predict(double *, double *);
    double*** parameters; // Array of parameters for each neuron
};

class CSANondimInputs{
private:
  int nDim;
public:
  CSANondimInputs(int);
  ~CSANondimInputs();
  void Set(SpalartAllmarasInputs*);
  void NondimensionalizeSource(int,double*);
  void DimensionalizeSource(int,double*);
  double Chi;
  double OmegaNondim;
  double OmegaBar;
  double SourceNondim;
  double NuGradNondim;
  double * DNuHatDXBar;
  double NuHatGradNorm;
  double NuHatGradNormBar;
};

#include "numerics_machine_learning.inl"
