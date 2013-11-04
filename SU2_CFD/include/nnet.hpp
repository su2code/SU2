#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>

using namespace std;

class Scaler{
public:
	Scaler();
	~Scaler();
	virtual void Scale(double *) = 0;
	virtual void Unscale(double *) = 0;
	//virtual void LoadJSON(int) = 0;
};

class NormalScaler: public Scaler{
private:
	double * mu;
	double * sigma;
	int nInputs;

public:
	NormalScaler();
	NormalScaler(int,double*,double*);
	~NormalScaler();
	NormalScaler(double *, double *);
	void Scale(double *);
	void Unscale(double *);
};

class Activator{
public:
	Activator();
	~Activator();
	virtual double Activate(double sum){cout<< "IN BASE ACTIVATOR THIS IS BAD" <<endl; return 0;};
};

class TanhActivator : public Activator{
public:
	TanhActivator();
	~TanhActivator();
	double Activate(double sum);
};

class LinearActivator : public Activator{
public:
	LinearActivator();
	~LinearActivator();
	double Activate(double sum);
};

class Neuron{
private:
	
	int parameterStart;
	int nParameters;
public:
	Activator *activator;
	Neuron();
	Neuron(Activator*, int, int); // activator, parameterStart, nParameters
	~Neuron();
	//Activator* GetActivator(void);
	double Combine(double * input, double * parameters);
	int ParameterStart();
	int ParameterEnd();
};

class CNeurNet {
private:
	int nNodes;
	int nInputs;
	int nOutputs;
	int nParameters;
	int nLayers;
	int *nNeuronsPerLayer;
	double *parameters;
	int **parameterIdx;
	int **nParametersPerNeuron;
	Neuron ***neurons; // Array of arrays to pointers to neuron
	Scaler* LoadScaler(ifstream&,int);
	int LoadInteger(ifstream&, string);
public:
	Scaler *inputScaler;
	Scaler *outputScaler;
	CNeurNet();
	CNeurNet(string, string);
	~CNeurNet();
	int NumInputs();
	int NumOutputs();
	void Predict(double *, double *);
	bool CheckPredictions(string);
};
