#include <vector>
#include <cstdlib>
#include <cmath>
#include <utility>

using namespace std;

class ActivationFunction{
    public:
        double sigmoid(double x){
            return (1/(1+exp(-x)));
        }
        double tanh(double x){
            return (exp(x)-exp(-x))/(exp(x)+exp(-x));
        }
        double relu(double x){
            if(x<0){
                return 0.0;
            }
            else{
                return x;
            }
        }
        double sigmoidDerivative(double x){
            return sigmoid(x)*(1-sigmoid(x));
        }
};

class Neuron{
    public:
        int bias;
        vector<double> weights;
        vector<double> lastInputs;
        int lastSum;
        Neuron(int numInputs){
            for (int i=0; i<numInputs; i++){
                weights.push_back(randomWeight());
            }
            bias = randomWeight();
        }
        
        double Activate(vector<double> inputs){
            double sum = bias;
            for(int i=0;i< weights.size(); i++){
                sum += weights.at(i) * inputs[i];
            }
            ActivationFunction funct;
            return funct.sigmoid(sum);
        } 
        void saveInputs(vector<double> inputs){
            lastInputs = inputs;
            lastSum = Activate(inputs);
        }

    private:
        double randomWeight() {
            return ((double) rand() / RAND_MAX) * 2 - 1; // Range: [-1, 1]
        }
    
};

// Base class for a layer in the neural network
class Layer{
    public:
        vector<Neuron> layer;

        Layer(int layerSize, int prevLayerSize){
            if (prevLayerSize==0){
                prevLayerSize = 1;
            }
            for(int i; i<layerSize; i++){
                layer.push_back(Neuron(prevLayerSize));
            }
        }
        virtual vector<double> feedForward(vector<double> inputs, bool training){
            vector<double> Activations;
            for(int i=0; i<inputs.size();i++){
                Activations.push_back(layer[i].Activate(inputs));
                if(training == true){
                    layer[i].saveInputs(inputs);
                }
            }
            return Activations;
        }
};

// Convolutional Layer
class ConvolutionalLayer : public Layer {

    public:
        int stride;
        vector<vector<double>> kernel;
        pair<int, int> kernalSize;

        ConvolutionalLayer(int numFilters, int filterSize, int inputSize, int stride, pair<int, int> kernalSize) : Layer(0,0) ,stride(stride), kernalSize(kernalSize.first, kernalSize.second)  { 
            for(int i; i<kernalSize.first; i++){
                for(int ii; ii<kernalSize.second; ii++){
                    kernel[i][ii] = ((double)rand() /RAND_MAX)*2.0 -1.0;
                }
            }
        }

        vector<double> feedForward(vector<double> inputs, bool training) override {
            
        }

};

