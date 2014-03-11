#pragma once

#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <sstream>

#ifndef NO_JSONCPP
#include <json/json.h>
#endif


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
#ifndef NO_JSONCPP
  CNormalScaler(Json::Value);
#endif
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
#ifndef NO_JSONCPP
  CTanhActivator(Json::Value);
#endif
	~CTanhActivator();
	double Activate(double combination);
};

class CLinearActivator : public CActivator{
public:
	CLinearActivator();
#ifndef NO_JSONCPP
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
#ifndef NO_JSONCPP
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
  ~CPredictor();
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

class CNeurNet : public CPredictor{
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
#ifndef NO_JSONCPP
  CNeurNet(Json::Value);
#endif
	~CNeurNet();
//	int InputDim();
//	int OutputDim();
	void Predict(double *, double *);
    double*** parameters; // Array of parameters for each neuron
};