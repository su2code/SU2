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
