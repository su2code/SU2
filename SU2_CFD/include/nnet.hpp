#pragma once

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>

using namespace std;

class CScaler{
public:
	CScaler();
  virtual ~CScaler();
	virtual void Scale(double *) = 0;
	virtual void Unscale(double *) = 0;
	//virtual void LoadJSON(int) = 0;
};

class CNormalScaler: public CScaler{
private:
	double * mu;
	double * sigma;
	int nInputs;

public:
	CNormalScaler();
	CNormalScaler(int,double*,double*);
	~CNormalScaler();
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
	~CTanhActivator();
	double Activate(double combination);
};

class CLinearActivator : public CActivator{
public:
	CLinearActivator();
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
	~CSumNeuron();
	//Activator* GetActivator(void);
	double Combine(double * parameters, int nParameters, double * inputs, int nInputs);
  double Activate(double combination);
};


class CPredictor{
public:
  CPredictor();
  ~CPredictor();
  virtual void Predict(double *, double *){cout << "In base Predict, this is bad";};
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

class CNeurNet : public CPredictor{
private:
  int maxNeurons; // Number of neurons in the largest layer
  int nLayers;
  CNeuron ***neurons; // Array of arrays to pointers to neuron
  double*** parameters; // Array of parameters for each neuron
  int* nNeuronsInLayer; //one list for each layer
  int** nParameters; // Number of parameters for the neuron
  int inputDim;
  
  void processLayer(double *, int,CNeuron **, double **, int, int * ,double *);
  
  //----
	int outputDim;
	int totalNumParameters;
public:
	CNeurNet();
    //CNeurNet(string filename); // Should have custom constructor with nInputs, etc.
	~CNeurNet();
//	int InputDim();
//	int OutputDim();
	void Predict(double *, double *);
};
