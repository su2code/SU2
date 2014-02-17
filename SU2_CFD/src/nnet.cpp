#include "../include/nnet.hpp"
using namespace std;

CScaler::CScaler(){}
CScaler::~CScaler(){}

CNormalScaler::CNormalScaler(){}
CNormalScaler::CNormalScaler(int nInputs, double *mu, double *sigma){
	this->nInputs = nInputs;
	this->mu = new double[nInputs];
	this->sigma = new double[nInputs];
	for (int i=0; i < nInputs; i++){
		this->mu[i] = mu[i];
		this->sigma[i] = sigma[i];
	}
	return;
}
CNormalScaler::~CNormalScaler(){
	delete [] mu;
	delete [] sigma;
}

void CNormalScaler::Scale(double * inputs){
	for (int i=0; i<nInputs; i++){
		inputs[i] = (inputs[i]-mu[i])/sigma[i];
	}
	return;
}
void CNormalScaler::Unscale(double * inputs){
	for (int i=0; i<nInputs; i++){
		inputs[i] = inputs[i]*sigma[i] + mu[i];
	}
	return;
}

CActivator::CActivator(){}
CActivator::~CActivator(){}

CTanhActivator::CTanhActivator(){}
CTanhActivator::~CTanhActivator(){}
double CTanhActivator::Activate(double sum){
	double val =  1.7159 * tanh(2.0/3.0 *sum);
	return val;
}

CLinearActivator::CLinearActivator(){}
CLinearActivator::~CLinearActivator(){}
double CLinearActivator::Activate(double sum){
	return sum;
}

CNeuron::CNeuron(){}
CNeuron::~CNeuron(){}

CSumNeuron::CSumNeuron(){}
CSumNeuron::CSumNeuron(CActivator* activator){
	this->activator = activator;
}
CSumNeuron::~CSumNeuron(){
	delete this->activator;
}

double CSumNeuron::Combine(double *parameters, int nParameters, double *inputs, int nInputs){
  if (nParameters != nInputs +1){
    cout << "parameter size mismatch" <<endl;
  }
  double combination = 0;
	for (int i = 0; i < nInputs; i++){
		combination += inputs[i] * parameters[i];
	}
	// Add in the bias term
  combination+= parameters[nParameters -1];
	return combination;
}

double CSumNeuron::Activate(double combination){
  return this->activator->Activate(combination);
}

CPredictor::CPredictor(){}
CPredictor::~CPredictor(){}

CScalePredictor::CScalePredictor(){}
CScalePredictor::CScalePredictor(string filename){
  cout << "filename is " << filename << endl;
  return;
}
CScalePredictor::~CScalePredictor(){
  delete this->Pred;
  delete this->InputScaler;
  delete this->OutputScaler;
  return;
}

CNeurNet::CNeurNet(){}
CNeurNet::~CNeurNet(){}

void CNeurNet::processLayer(double * input, int nInput, CNeuron **neurons, double **parameters, int nNeurons, int * nParameters, double *output){
  for (int i = 0; i < nNeurons; i++){
    double combination = neurons[i]->Combine(parameters[i], nParameters[i], input, nInput);
    output[i] = neurons[i]->Activate(combination);
  }
  return;
}

void CNeurNet::Predict(double * input, double * output){
  double * prevTmpOutput = new double[this->maxNeurons];
  double *tmpOutput = new double[this->maxNeurons];
  
  int nLayers = this->nLayers;
  if (nLayers == 1){
    this->processLayer(input, this->inputDim, this->neurons[0], this->parameters[0], this->nNeuronsInLayer[0], this->nParameters[0], output);
    return;
  }
  
  // First layer uses the real input as the input
  this->processLayer(input, this->inputDim, this->neurons[0], this->parameters[0], this->nNeuronsInLayer[0], this->nParameters[0], tmpOutput);
  
  // Middle layers use the previous output as input
  for (int i= 1; i < nLayers -1; i++){
    double *tmp;
    tmp = prevTmpOutput;
    prevTmpOutput = tmpOutput;
    tmpOutput = tmp;
    
    int inputDim = this->nNeuronsInLayer[i-1];
    processLayer(prevTmpOutput, inputDim, this->neurons[i], this->parameters[i], this->nNeuronsInLayer[i], this->nParameters[i], tmpOutput);
  }
  int layer = nLayers -1;
  int inputDim = this->nNeuronsInLayer[nLayers-1];
  // Last layer has the actual output
  processLayer(tmpOutput, inputDim, this->neurons[layer], this->parameters[layer], this->nNeuronsInLayer[layer],this->nParameters[layer], output);
  return;
}

/*
void CNeurNet::Predict(double * input, double * output){
	this->inputScaler->Scale(input);
  
	// Predict over all the layers
  
	// Make data to store the layers
	//double ** neuronCombinations = new double*[this->nLayers];
	double ** neuronOutputs = new double*[this->nLayers];
	for (int i = 0; i < this->nLayers; i++){
		neuronOutputs[i] = new double[this -> nNeuronsPerLayer[i]];
	}
  
	// Make the prediction
	// The first layer has the input as an input
	int firstLayer = 0;
	for (int i= 0; i < nNeuronsPerLayer[firstLayer]; i++){
		double combination = this->neurons[firstLayer][i]->Combine(input, this->parameters);
		double output = this->neurons[firstLayer][i]->activator->Activate(combination);
		neuronOutputs[firstLayer][i] = output;
	}
  
	// For all the other layers, the input is the output of the previous layer
	for (int j = 1; j < this->nLayers; j++){
		for (int i=0; i < nNeuronsPerLayer[j]; i++){
			//cout << "Layer " << j << " neuron " << i << endl;
			double combination = this->neurons[j][i]->Combine(neuronOutputs[j-1], this->parameters);
			//cout << "combination " << combination <<endl;
			double output = this->neurons[j][i]->activator->Activate(combination);
			//cout << "output " << output <<endl;
			neuronOutputs[j][i] = output;
		}
	}
  
	// Copy the last layer
	int lastLayer = this->nLayers - 1;
	for (int i = 0; i < this->nOutputs; i++){
		output[i] = neuronOutputs[lastLayer][i];
	}
  
	this->inputScaler->Unscale(input);
	this->outputScaler->Unscale(output);
  
	for (int i=0; i < this->nLayers; i++){
		//delete neuronSums[i];
		delete [] neuronOutputs[i];
	}
	delete [] neuronOutputs;
	return;
}


CNeurNet::CNeurNet(){}

CNeurNet::CNeurNet(string filename, string predFilename){
	// Open the file
	ifstream f;
	string text_line;
	f.open(filename.c_str(), ios::in);

	// Construct the neural net
	// this algorithm is slow, but should be 
	// good enough for now

	// If you are a future person reading this, be careful, this
	// code not at all robust to changes in the JSON format

	int position;
	this->nInputs = this->LoadInteger(f, "\"NumInputs\":");
	f.clear();
	f.seekg(0,ios::beg);
	this->nOutputs = this->LoadInteger(f, "\"NumOutputs\":");
	f.clear();
	f.seekg(0,ios::beg);
	this->nParameters = this->LoadInteger(f, "\"TotalNumParameters\":");
	f.clear();
	f.seekg(0,ios::beg);
	this->nLayers = this->LoadInteger(f, "\"NumLayers\":");
	f.clear();
	f.seekg(0,ios::beg);
	//cout << "Integers read"<<endl;

	// Read in the number of neurons per layer
	this->nNeuronsPerLayer = new int[this->nLayers];
	bool nNeuronsRead = false;
	while(getline(f,text_line)){
		position = text_line.find("\"NumNeuronsPerLayer\":");
		if (position == string::npos){
			continue;
		}
		nNeuronsRead = true;
		for (int i = 0; i < this->nLayers; i++){
			getline(f,text_line);
			stringstream ss(text_line);
			int neur = 0;
			string tmp;
			ss >> neur >> tmp;
			//cout << "neuron in layer: " << neur <<endl;
			this->nNeuronsPerLayer[i] = neur;
		}
	}
	if (!nNeuronsRead){
		cout << "Neurons not read"<<endl;
		throw 10;
	}
	f.clear();
	f.seekg(0, ios::beg);

	// Now that we have those, create memory for the parameters
	this->parameters = new double[nParameters];
	this->parameterIdx = new int*[this->nLayers];
	this->nParametersPerNeuron = new int*[this->nLayers];
	for (int i=0; i < this->nLayers; i++){
		this->parameterIdx[i] = new int[this->nNeuronsPerLayer[i]];
		this->nParametersPerNeuron[i] = new int[this->nNeuronsPerLayer[i]];
	}

	//cout << "before while loop"<<endl;
	// Read in the rest of the data
	bool inputScalerFound = false;
	bool outputScalerFound = false;
	while (getline(f,text_line)){
		position = text_line.find("InputScaler");
		if (position != string::npos){
			inputScalerFound = true;
			this->inputScaler = LoadScaler(f, this->nInputs);
		}
		position = text_line.find("OutputScaler");
		if (position != string::npos){
			outputScalerFound = true;
			this->outputScaler = LoadScaler(f, this->nOutputs);
		}
		position = text_line.find("\"Parameters\":");
		if (position != string::npos){
			for (int i=0; i< this->nParameters; i++){
				getline(f,text_line);
				stringstream ss(text_line);
				string tmp;
				double param;
				ss >> param >> tmp;
				this->parameters[i] = param;
			}
		}
		position = text_line.find("\"NumParametersPerNeuron\":");
		if (position != string::npos){
			for (int i = 0; i < this->nLayers; i++){
				// Get the line from the first brace
				getline(f,text_line);
				for (int j = 0; j < this->nNeuronsPerLayer[i]; j++){
					getline(f,text_line);
					stringstream ss(text_line);
					int neur = -1 ;
					ss >> neur;
					if (neur == -1){
						cout << "Misread nParam: text is: " <<text_line<<endl;
						throw 10;
					}
					this->nParametersPerNeuron[i][j] = neur;
				}
				// Get the ending brace
				getline(f,text_line);
			}
		}
		position = text_line.find("\"ParameterIndex\":");
		if (position != string::npos){
			for (int i = 0; i < this->nLayers; i++){
				// Get the line from the first brace
				getline(f,text_line);
				for (int j = 0; j < this->nNeuronsPerLayer[i]; j++){
					getline(f,text_line);
					stringstream ss(text_line);
					int neur = -1 ;
					ss >> neur;
					if (neur == -1){
						cout << "Misread param idx: text is: " <<text_line<<endl;
						throw 10;
					}
					this->parameterIdx[i][j] = neur;
				}
				// Get the ending brace
				getline(f,text_line);
			}
		}
	}

	// Create all neurons. Right now they must all be Tanh except the last
	// layer which is linear. Should check this later
	this->neurons = new Neuron**[this->nLayers];
	for (int i = 0; i < this->nLayers; i++){
		this->neurons[i] = new Neuron*[this->nNeuronsPerLayer[i]];
	}
	for (int i = 0; i < this->nLayers - 1; i++){
		for (int j = 0; j < this->nNeuronsPerLayer[i]; j++){
			Activator *activator = new TanhActivator;
			this->neurons[i][j] = new Neuron(activator, this->parameterIdx[i][j], this->nParametersPerNeuron[i][j]);
		}
	}
	int lastLayer = this->nLayers -1;
	for (int j = 0; j < this->nNeuronsPerLayer[lastLayer]; j++){
		Activator * activator = new LinearActivator;
		this -> neurons[lastLayer][j] = new Neuron(activator, this->parameterIdx[lastLayer][j], this->nParametersPerNeuron[lastLayer][j]);
	}

	// Check predictions
	this->CheckPredictions(f);
  
  f.close();
	return;
}

CNeurNet::~CNeurNet(){
//	delete this->inputScaler;
//	delete this->outputScaler;
	delete this->parameters;
	for (int i=0; i < this->nLayers; i++){
		delete [] this->parameterIdx[i];
		delete [] this->nParametersPerNeuron[i];
		for (int j = 0; j < this -> nNeuronsPerLayer[i]; j++){
			delete [] this -> neurons[i][j];
		}
	}
	for (int i= 0; i < this->nLayers; i++){
		delete [] this->neurons[i];
	}
	delete [] this->nNeuronsPerLayer;
	delete [] this->parameterIdx;
	delete [] this->nParametersPerNeuron;
	delete [] this->neurons;
}

int CNeurNet::LoadInteger(ifstream& f, string name){
	int position;
	int i;
	string text_line;
	bool stringFound = false;
	while(getline(f,text_line)){
		position = text_line.find(name);
		if (position != string::npos){
			stringstream ss(text_line);
			string tmpstr;
			string comma;
			ss >> name >> i >> comma;
			stringFound = true;
		}
	}
	if (!stringFound){
		cout << "String "<< name <<" not found"<<endl;
    exit(1);
	}
	return i;
}

Scaler* CNeurNet::LoadScaler(ifstream& f, int nVals){
	string text_line;
	getline(f, text_line);
	int position;
	// Next line should contain the go package path
	position = text_line.find("Type");
	if (position == string::npos){
		cout <<"Type not on line after InputScaler"<<endl;
		throw 10;
	}
	position = text_line.find("github.com/btracey/nnet/scale/");
	if (position == string::npos){
		cout << "Scaler must come from scale"<<endl;
		throw 10;
	}
	// This could be replaced with some sort of string parsing and switch statement
	position = text_line.find("Normal");
	if (position != string::npos){
		// Net line should have Value
		getline(f, text_line);
		// Next line should be the start of Mu
		getline(f, text_line);
		position = text_line.find("Mu");
		if (position == string::npos){
			cout << "Mu not found reading normal scaler"<< endl;
			throw 10;
		}
	
		// If it is there, should be nInputs lines of doubles
		double *mus = new double[nVals];
		for (int i=0; i< nVals; i++){
			// Get the line
			getline(f, text_line);
			stringstream ss(text_line);
			double mu;
			string tmp2;
			ss >> mu >> tmp2;
			mus[i] = mu;
		}
		getline(f, text_line);
		// Now, should read in the sigmas
		double *sigmas = new double[nVals];
		getline(f,text_line);
		position = text_line.find("Sigma");
		if (position == string::npos){
			throw 10;
		}
		for (int i = 0; i<nVals; i++){
			getline(f,text_line);
			stringstream ss(text_line);
			double sigma;
			string tmp;
			ss >> sigma >> tmp;
			sigmas[i] = sigma;
		}

		Scaler * s = new NormalScaler(nVals,mus,sigmas);

		delete [] mus;
		delete [] sigmas;

		return s;
	}
	throw "Scaler type not implemented";
}

void CNeurNet::Predict(double * input, double * output){
	this->inputScaler->Scale(input);

	// Predict over all the layers

	// Make data to store the layers
	//double ** neuronCombinations = new double*[this->nLayers];
	double ** neuronOutputs = new double*[this->nLayers];
	for (int i = 0; i < this->nLayers; i++){
		neuronOutputs[i] = new double[this -> nNeuronsPerLayer[i]];
	}

	// Make the prediction
	// The first layer has the input as an input
	int firstLayer = 0;
	for (int i= 0; i < nNeuronsPerLayer[firstLayer]; i++){
		double combination = this->neurons[firstLayer][i]->Combine(input, this->parameters);
		double output = this->neurons[firstLayer][i]->activator->Activate(combination);
		neuronOutputs[firstLayer][i] = output;
	}

	// For all the other layers, the input is the output of the previous layer
	for (int j = 1; j < this->nLayers; j++){
		for (int i=0; i < nNeuronsPerLayer[j]; i++){
			//cout << "Layer " << j << " neuron " << i << endl;
			double combination = this->neurons[j][i]->Combine(neuronOutputs[j-1], this->parameters);
			//cout << "combination " << combination <<endl;
			double output = this->neurons[j][i]->activator->Activate(combination);
			//cout << "output " << output <<endl;
			neuronOutputs[j][i] = output;
		}
	}

	// Copy the last layer 
	int lastLayer = this->nLayers - 1;
	for (int i = 0; i < this->nOutputs; i++){
		output[i] = neuronOutputs[lastLayer][i];
	}

	this->inputScaler->Unscale(input);
	this->outputScaler->Unscale(output);

	for (int i=0; i < this->nLayers; i++){
		//delete neuronSums[i];
		delete [] neuronOutputs[i];
	}
	delete [] neuronOutputs;
	return;
}

int CNeurNet::NumInputs(){
	return nInputs;
}

int CNeurNet::NumOutputs(){
	return nOutputs;
}

bool CNeurNet::CheckPredictions(ifstream& f){
  int position;
  string text_line;
  bool lineFound = false;
  // Rewind to top
  f.clear();
	f.seekg(0, ios::beg);
  
  while (getline(f,text_line)){
		position = text_line.find("\"PredictionCheck\":");
		if (position != string::npos){
      lineFound = true;
      break;
		}
  }
  if (!lineFound){
    cout << "PredictionCheck line not found" << endl;
    exit(10);
  }
  
	// Check that the predictions from this net match the predictions 
	// from the file
	// Should be a "["
	double * input = new double[this->nInputs];
	double * output = new double[this->nOutputs];
	double * predOutput = new double[this->nOutputs];
	string tmp;
	double val = 0;

	while (true){
		// Read starting { or ending ]
		getline(f,text_line);
		position = text_line.find("]");
		if (position != string::npos){
			break;
		}
		// Read inputs line
		getline(f, text_line);
		position = text_line.find("\"Input\":");
		if (position == string::npos){
			cout <<"No Input line"<<endl;
			cout << "Line is: " << text_line<<endl;
			exit(10);
		}
		for (int i = 0; i < this->nInputs; i++){
			getline(f,text_line);
			stringstream ss(text_line);
			// Read in double
			ss >> val;
			input[i] = val;
		}
		// Read ending braces
		getline(f,text_line);
		// Read outputs line
		getline(f,text_line);
		position = text_line.find("\"Output\":");
		if (position == string::npos){
			cout <<"No Output line"<<endl;
			cout << "Line is: " << text_line<<endl;
			throw 10;
		}
		for (int i=0; i < this->nOutputs; i++){
			getline(f,text_line);
			stringstream ss(text_line);
			ss >> val;
			output[i] = val;
		}
		// Get the ]
		getline(f,text_line);
		// get the }
		getline(f, text_line);

		// Check that the prediction matches
		this->Predict(input,predOutput);
		for (int i = 0 ; i < this->nOutputs; i++){
			if (abs(output[i] - predOutput[i]) > 10e-14){
				cout << "predictions don't match"<< endl;
				cout << "real output " << output[i] << endl;
			cout << "pred output " << predOutput[i] << endl;
			cout << "diff = "<< output[i] - predOutput[i]<<endl;
			exit(10);
			}
		}
	}
	delete [] input;
	delete [] output;
	delete [] predOutput;
	return true;
}
*/