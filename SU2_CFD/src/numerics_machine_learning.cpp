/*!
 * \file numerics_template.cpp
 * \brief This file contains all the convective term discretization.
 * \author B. Tracey
 * \version 3.2.8.1 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
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

#include "../include/numerics_machine_learning.hpp"
using namespace std;

// This is up here because otherwise it's not in scope of the functions below
#ifdef HAVE_JSONCPP
CPredictor* parse_predictor(Json::Value json) {
  string type = json["Type"].asString();
  Json::Value value = json["Value"];
  if (type.compare("github.com/reggo/reggo/supervised/nnet/Net*")==0) {
    CPredictor* predictor = new CNeurNet(value);
    return predictor;
  }
  if (type.compare("github.com/btracey/ransuq/mlalg/MulPredictor")==0) {
    CPredictor* predictor = new CMulPredictor(value);
    return predictor;
  }
  cout << "No Match for predictor type: " << type << endl;
  return NULL;
}

CScaler* parse_cscaler(Json::Value json) {
  
  string type = json["Type"].asString();
  Json::Value value = json["Value"];
  if (type.compare("github.com/reggo/reggo/scale/Normal*") == 0) {
    // We matched the normal scaler. Now, allocate a new one
    CScaler * scaler = new CNormalScaler(value);
    return scaler;
  } else if (type.compare("github.com/btracey/ransuq/mlalg/MulInputScaler*") == 0) {
    // Allocate a MulScaler
    CScaler * scaler = new CMulInputScaler(value);
    return scaler;
  } else if (type.compare("github.com/btracey/ransuq/mlalg/MulOutputScaler*") == 0) {
    // Allocate a MulScaler
    CScaler * scaler = new CMulOutputScaler(value);
    return scaler;
  } else{
    cout << "NoMatch for scaler type: "<< type << endl;
    exit(EXIT_FAILURE);
  }
  cout << "Shouldnt' be here" << endl;
  exit(EXIT_FAILURE);
  
  return NULL;
}
#endif



CScaler::CScaler() {}
CScaler::~CScaler() {}

CNormalScaler::CNormalScaler() {}
CNormalScaler::CNormalScaler(int dim, double *mu, double *sigma) {
	this->dim = dim;
	this->mu = new double[dim];
	this->sigma = new double[dim];
	for (int i=0; i < dim; i++) {
		this->mu[i] = mu[i];
		this->sigma[i] = sigma[i];
	}
	return;
}
CNormalScaler::~CNormalScaler() {
	delete [] mu;
	delete [] sigma;
}

#ifdef HAVE_JSONCPP
CNormalScaler::CNormalScaler(Json::Value json) {
  int nDim = json["Dim"].asInt();
  Json::Value muVal = json["Mu"];
  
  int muSize = muVal.size();
  if (muSize != nDim) {
    cout << "musize and Dim mismatch" << endl;
  }
  double *mu = new double[muSize];
  for (int i = 0; i < nDim; i++) {
    mu[i] = muVal[i].asDouble();
  }

  Json::Value sigmaVal = json["Sigma"];
  int sigmaSize = sigmaVal.size();
  if (sigmaSize != nDim) {
    cout << "sigmasize and and Dim mismatch" << endl;
  }
  
  double * sigma = new double[sigmaSize];
  for (int i = 0; i < nDim; i++) {
    sigma[i] = sigmaVal[i].asDouble();
  }

  this->dim = nDim;
  this->mu = mu;
  this->sigma = sigma;
  
  return;
}
#endif


void CNormalScaler::Scale(double * inputs) {
	for (int i=0; i<dim; i++) {
		inputs[i] = (inputs[i]-mu[i])/sigma[i];
	}
	return;
}
void CNormalScaler::Unscale(double * inputs) {
	for (int i=0; i<dim; i++) {
		inputs[i] = inputs[i]*sigma[i] + mu[i];
	}
	return;
}

CMulInputScaler::CMulInputScaler() {}

#ifdef HAVE_JSONCPP
CMulInputScaler::CMulInputScaler(Json::Value json) {
  // Get the scaler subfield
  this->InnerScaler = parse_cscaler(json["Scaler"]);

  // Get the constant
  this->MulScale = json["MulOutputScaler"]["MulScale"].asDouble();
}
#endif

CMulInputScaler::~CMulInputScaler() {}
void CMulInputScaler::Scale(double * inputs) {
  inputs[0] /= this->MulScale;
  double * second = &inputs[1];
  this->InnerScaler->Scale(second);
  return;
}

void CMulInputScaler::Unscale(double * inputs) {
  inputs[0] *= this->MulScale;
  double * second = &inputs[1];
  this->InnerScaler->Unscale(second);
  return;
}

CMulOutputScaler::CMulOutputScaler() {}
CMulOutputScaler::~CMulOutputScaler() {}
void CMulOutputScaler::Scale(double * outputs) {
  outputs[0] /= this->MulScale;
  return;
}

#ifdef HAVE_JSONCPP
CMulOutputScaler::CMulOutputScaler(Json::Value json) {
  this->MulScale = json["MulScale"].asDouble();
}
#endif

void CMulOutputScaler::Unscale(double * outputs) {
  outputs[0] *= this->MulScale;
  return;
}

CActivator::CActivator() {}
CActivator::~CActivator() {}

CTanhActivator::CTanhActivator() {}
#ifdef HAVE_JSONCPP
CTanhActivator::CTanhActivator(Json::Value json) {
}
#endif
CTanhActivator::~CTanhActivator() {}
double CTanhActivator::Activate(double sum) {
  //cout << "In CTanhActivator" << endl;
	double val =  1.7159 * tanh(2.0/3.0 *sum);
  //cout << "Done CTanhActivator" << endl;
	return val;
}

CLinearActivator::CLinearActivator() {}
#ifdef HAVE_JSONCPP
CLinearActivator::CLinearActivator(Json::Value) {
}
#endif

CLinearActivator::~CLinearActivator() {}

double CLinearActivator::Activate(double sum) {
//  cout << "In CLinear Activate" << endl;
	return sum;
}

CNeuron::CNeuron() {}
CNeuron::~CNeuron() {}

CSumNeuron::CSumNeuron() {}
CSumNeuron::CSumNeuron(CActivator* activator) {
	this->activator = activator;
}
#ifdef HAVE_JSONCPP
CSumNeuron::CSumNeuron(Json::Value json) {
  string type = json["Type"].asString();
  if (type.compare("github.com/reggo/reggo/supervised/nnet/Tanh") == 0) {
    this->activator = new CTanhActivator(json["Value"]);
  } else if (type.compare("github.com/reggo/reggo/supervised/nnet/Linear") == 0) {
    this->activator = new CLinearActivator(json["Value"]);
  } else{
    cout << "Unknown activator type: " << type << endl;
    exit(EXIT_FAILURE);
  }
}
#endif

CSumNeuron::~CSumNeuron() {
	delete this->activator;
}

double CSumNeuron::Combine(double *parameters, int nParameters, double *inputs, int nInputs) {
  double combination = 0;
	for (int i = 0; i < nInputs; i++) {
		combination += inputs[i] * parameters[i];
	}
 
	// Add in the bias term
  combination+= parameters[nParameters -1];
//  cout << "Done CSumNeuronCombine" << endl;
	return combination;
}

double CSumNeuron::Activate(double combination) {
  return this->activator->Activate(combination);
}

CPredictor::CPredictor() {}
CPredictor::~CPredictor() {}

int CPredictor::InputDim() {
  return this->inputDim;
}

int CPredictor::OutputDim() {
  return this->outputDim;
}

CMulPredictor::CMulPredictor() {}
#ifdef HAVE_JSONCPP
CMulPredictor::CMulPredictor(Json::Value json) {
  this->Inner = parse_predictor(json["Inner"]);
  this->inputDim = this->Inner->InputDim() + 1;
  this->outputDim = this->Inner->OutputDim();
}
#endif

CMulPredictor::~CMulPredictor() {}

void CMulPredictor::Predict(double * input, double * output) {
  double * secondInput = &input[1];
  this->Inner->Predict(secondInput, output);
  for (int i = 0; i < this->OutputDim(); i++) {
    output[i] *= input[0];
  }
  return;
}

CNeurNet::CNeurNet() {}

#ifdef HAVE_JSONCPP
CNeurNet::CNeurNet(Json::Value json) {
  this-> inputDim = json["InputDim"].asInt();
  this-> outputDim = json["OutputDim"].asInt();
  this-> totalNumParameters = json["TotalNumParameters"].asInt();

  // Get the number of neurons
  Json::Value layers = json["Neurons"];
  Json::Value parameters = json["Parameters"];
  int nLayers = layers.size();
  int nParameterLayers = parameters.size();
  if (nLayers != nParameterLayers) {
    cout << "neurons and parameters disagree on number of layers" << endl;
  }
  this->nLayers = nLayers;
  
  // Allocate memory for neuron and parameters by layer
  this->nNeuronsInLayer = new int[nLayers];
  this->nParameters = new int *[nLayers];
  this->neurons = new CNeuron **[nLayers];
  this->parameters = new double**[nLayers];
  
  // Per layer, get the number of neurons in the layer and then read in the neurons
  for (int i = 0; i < nLayers; i++) {
    Json::Value layer = layers[i];
    Json::Value parameterLayer = parameters[i];
    int neuronsInLayer = layer.size();
    int parameterNeurons = parameterLayer.size();
    if (neuronsInLayer != parameterNeurons) {
      cout << "neurons and parameters disagree on the number of neurons in layer i" << endl;
    }
    this->nNeuronsInLayer[i] = neuronsInLayer;
    if (i == nLayers-1) {
      if (neuronsInLayer != this->outputDim) {
        cout << "Size of initial layer is not equal to input dimension" << endl;
      }
    }
    this->neurons[i] = new CNeuron *[neuronsInLayer];
    this->nParameters[i] = new int[neuronsInLayer];
    this->parameters[i] = new double *[neuronsInLayer];
    
    // Loop over all the neurons in the layer and add the parameters and neuron itself
    for (int j = 0; j < neuronsInLayer; j++) {
      Json::Value neuron = layer[j];
      
      // get the parameters
      int nParametersInNeuron = parameterLayer[j].size();
      this->nParameters[i][j] = nParametersInNeuron;
      this->parameters[i][j] = new double [nParametersInNeuron];
      for (int k = 0; k < nParametersInNeuron; k++) {
        
        double val = parameterLayer[j][k].asDouble();
        this->parameters[i][j][k] = val;
        /*
        long *ptrFloatAsInt =(long *)( &val);
        cout.precision(16);
        cout << hex;
        cout << "params: i = " << i << "\tj = " << j << "\tk = " << k << "\tparam = " << *ptrFloatAsInt << endl;
         */
      }
      
      // get the neurons
      string type = neuron["Type"].asString();
      if (type.compare("github.com/reggo/reggo/supervised/nnet/SumNeuron") == 0) {
        this->neurons[i][j] = new CSumNeuron(neuron["Value"]);
      } else{
        cout << "neuron type unknown: " << type << endl;
      }
    }
    // TODO: Should add in extra checking for num parameters and such
  }
  
  // Find the maximum number of neurons
  this->maxNeurons = 0;
  for (int i = 0; i < this->nLayers; i++) {
    if (this->maxNeurons < this->nNeuronsInLayer[i]) {
      this->maxNeurons = this->nNeuronsInLayer[i];
    }
  }
}
#endif

CNeurNet::~CNeurNet() {
  
  for (int i = 0; i < this->nLayers; i++) {
    for (int j = 0; j < this->nNeuronsInLayer[i]; j++) {
      delete [] this-> parameters[i][j];
      delete [] this->neurons[i][j];
    }
    delete [] this->nParameters[i];
    delete [] this->parameters[i];
    delete [] this->neurons[i];
  }
  delete [] this->nNeuronsInLayer;
  delete [] this->neurons;
  delete [] this->parameters;
  delete [] this->nParameters;
}

void CNeurNet::processLayer(double * input, int nInput, CNeuron **neurons, double **parameters, int nNeurons, int * nParameters, double *output) {
  for (int i = 0; i < nNeurons; i++) {
    double combination = neurons[i]->Combine(parameters[i], nParameters[i], input, nInput);
    output[i] = neurons[i]->Activate(combination);
  }
  return;
}

void CNeurNet::Predict(double * input, double * output) {
  double * prevTmpOutput = new double[this->maxNeurons];
  double *tmpOutput = new double[this->maxNeurons];
  
  int nLayers = this->nLayers;
  
  /*
  cout << "Input is:" << endl;
  for (int i=0; i < this->InputDim(); i++) {
    cout << input[i] << endl;
  }
  cout << "Done input" << endl;
   */
  
  if (nLayers == 1) {
    this->processLayer(input, this->inputDim, this->neurons[0], this->parameters[0], this->nNeuronsInLayer[0], this->nParameters[0], output);
    return;
  }
  
  // First layer uses the real input as the input
  this->processLayer(input, this->inputDim, this->neurons[0], this->parameters[0], this->nNeuronsInLayer[0], this->nParameters[0], tmpOutput);
  
  /*
  cout << "First layer tmp output" << endl;
  for (int i = 0; i < this->nNeuronsInLayer[0]; i++) {
    cout << tmpOutput[i] << endl;
    double tmp = tmpOutput[i];
    long *ptrFloatAsInt =(long *)( &tmp);
    cout << hex;
    cout << *ptrFloatAsInt << endl;
    cout << dec;
  }
  cout << "done tmp output" << endl;
   */
  
  // Middle layers use the previous output as input
  for (int i= 1; i < nLayers -1; i++) {
    double *tmp;
    tmp = prevTmpOutput;
    prevTmpOutput = tmpOutput;
    tmpOutput = tmp;
    
    int inputDim = this->nNeuronsInLayer[i-1];
    processLayer(prevTmpOutput, inputDim, this->neurons[i], this->parameters[i], this->nNeuronsInLayer[i], this->nParameters[i], tmpOutput);
  }
  int layer = nLayers -1;
  int inputDim = this->nNeuronsInLayer[nLayers-2];
  
  // Last layer has the actual output
  processLayer(tmpOutput, inputDim, this->neurons[layer], this->parameters[layer], this->nNeuronsInLayer[layer],this->nParameters[layer], output);
  
  // Clean up garbage
  delete [] prevTmpOutput;
  delete [] tmpOutput;
  return;
}



// get_file_contents gets all of the file contents and returns them as a string
string get_file_contents(string filename) {
  
  const char * charfile = filename.c_str();
  
  ifstream in(charfile, ios::in | ios::binary);
  string contents;
  if (in)
  {
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());
    in.close();
    return(contents);
  }
  cout << "Predictor filename " << filename << " not found" << endl;
  exit(EXIT_FAILURE);
}

// TODO: Separate filename from parse script. (make a function of a Node)
CScalePredictor::CScalePredictor() {}
#ifdef HAVE_JSONCPP
CScalePredictor::CScalePredictor(string filename) {
  
  if (filename.compare("none")==0) {
    return;
  }
  
  string contents = get_file_contents(filename);
  
  Json::Value root;
  Json::Reader reader;
  bool parsingSuccessful = reader.parse(contents, root);
  if (!parsingSuccessful) {
    std::cout << "Failed to parse \n" << reader.getFormatedErrorMessages()<< endl;
  }
  
  // Get the input scaler
  this->InputScaler = parse_cscaler(root["InputScaler"]);
  this->OutputScaler = parse_cscaler(root["OutputScaler"]);
  this->Pred = parse_predictor(root["Predictor"]);
  
  // Check the predictions
  int nTestInputs = root["TestInputs"].size();
  int nTestOutputs = root["TestOutputs"].size();
  if (nTestInputs != nTestOutputs) {
    cout << "Number of test inputs and number of test outputs doesn't match" << endl;
  }
  Json::Value testInputs = root["TestInputs"];
  Json::Value testOutputs = root["TestOutputs"];
  
  for (int i = 0; i < nTestInputs; i++) {
    int nInputs = testInputs[i].size();
    int nOutputs = testOutputs[i].size();
    double *input = new double[nInputs];
    double *output = new double[nOutputs];
    for (int j = 0; j < nInputs; j++) {
      input[j] = testInputs[i][j].asDouble();
    }
    for (int j=0; j < nOutputs; j++) {
      output[j] = testOutputs[i][j].asDouble();
    }
    double *predOutput = new double[nOutputs];
    
    this->Predict(input, predOutput);
    bool mismatch = 0;
    for (int j = 0; j < nOutputs; j++) {
      double max = abs(output[j]);
      if (predOutput[j] > max) {
        max = abs(predOutput[j]);
      }
      if (max < 1.0) {
        max = 1.0;
      }
      if (abs(output[j] - predOutput[j])/(max) > 1e-12) {
        mismatch = 1;
      }
      cout.precision(16);
      if (mismatch) {
        cout << "Prediction mismatch" << endl;
        for (int j = 0; j < nOutputs; j++) {
          cout << "j = " <<  " true: " << output[j] << " pred: " << predOutput[j] << " rel error: " << abs(output[j] - predOutput[j])/(max) <<  endl;
        }
        throw string("mismatch");
      }
    }
//    if (mismatch) {
//      throw(-1);
//    }

    delete [] predOutput;
    delete [] input;
    delete [] output;
  }
  return;
}
#else
  CScalePredictor::CScalePredictor(string filename) {
    cout << "Must have JsonCpp installed" << endl;
    exit(EXIT_FAILURE);
  }
#endif
CScalePredictor::~CScalePredictor() {
  delete this->Pred;
  delete this->InputScaler;
  delete this->OutputScaler;
  return;
}
void CScalePredictor::Predict(double *input, double *output) {
  // Scale the input
  this->InputScaler->Scale(input);
  
  // Call the predict method
  this->Pred->Predict(input, output);
  
  // Unscale
  this->InputScaler->Unscale(input);
	this->OutputScaler->Unscale(output);
}

CSANondimInputs::CSANondimInputs(int nDim) {
  this->nDim = nDim;
  this->DNuHatDXBar = new double[nDim];
};
CSANondimInputs::~CSANondimInputs() {
  delete [] this->DNuHatDXBar;
};

void CSANondimInputs::Set(SpalartAllmarasInputs* sainputs) {
  double kin_visc = sainputs->Laminar_Viscosity / sainputs->Density;
  this->Chi = sainputs->Turbulent_Kinematic_Viscosity / kin_visc;
  double nuSum = sainputs->Turbulent_Kinematic_Viscosity + kin_visc;
  double dist = sainputs->dist;
  if (dist < 1e-10) {
    dist = 1e-10;
  }
  
  double distSq = dist* dist;
  
  this->OmegaNondim = 1.0 / ( distSq / (nuSum) );
  this->OmegaBar = sainputs->Omega / this->OmegaNondim;
  this->SourceNondim = 1.0 /( distSq / (nuSum * nuSum) );
  this->NuGradNondim = 1.0 /( dist / nuSum );
  
  if ((this->NuGradNondim == this->NuGradNondim) &&
      ((this->NuGradNondim - this->NuGradNondim) != 0.0)) {
    cout << this->NuGradNondim;
    cout << "Inf NuGrad" << endl;
    cout << "dist = " << dist << endl;
    cout << "distSq = " << distSq << endl;
    cout << "nuSum = " << nuSum << endl;
    exit(EXIT_FAILURE);
  }
  
  double * turbViscGrad = sainputs->GetTurbKinViscGradient();
  this->NuHatGradNorm = 0;
  for (int i = 0; i < this->nDim; i++) {
    this->DNuHatDXBar[i] = turbViscGrad[i] / this->NuGradNondim;
    this->NuHatGradNorm += turbViscGrad[i] * turbViscGrad[i];
  }
  this->NuHatGradNormBar = this->NuHatGradNorm / this->SourceNondim;
};
void CSANondimInputs::NondimensionalizeSource(int nResidual, double * source) {
  for (int i = 0; i < nResidual; i++) {
    source[i] /= this->SourceNondim;
  }
}

void CSANondimInputs::DimensionalizeSource(int nResidual, double * source) {
  
  for (int i = 0; i < nResidual; i++) {
    //    cout << "pre" << source[i] << endl;
    source[i] *= this->SourceNondim;
    //    cout << "post" << source[i] << endl;
  }
}
