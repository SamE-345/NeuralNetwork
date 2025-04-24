#include <vector>
#include <cstdlib>
#include <cmath>
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
};

class Neuron{
    public:
        int bias;
        vector<double> weights;

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

    private:
        double randomWeight() {
            return ((double) rand() / RAND_MAX) * 2 - 1; // Range: [-1, 1]
        }
    
};

class Layer{
    public:
        static vector<Neuron> layer;

        Layer(int layerSize, int prevLayerSize){
            if (prevLayerSize==0){
                prevLayerSize = 1;
            }
            for(int i; i<layerSize; i++){
                layer.push_back(Neuron(prevLayerSize));
            }
        }
        vector<double> feedForward(vector<double> inputs){
            vector<double> Activations;
            for(int i=0; i< sizeof(inputs);i++){
                Activations.push_back(layer[i].Activate(inputs));
            }
        }
};


class NeuralNetwork{
    public:
        vector<Layer> Layers;
        NeuralNetwork(int layerSizes[]){
            int x = 0;
            for(int i=0; i<sizeof(layerSizes); i++){
                Layers.push_back(Layer(layerSizes[i],x));
                x=layerSizes[i];
            }
        }
        void trainModel(vector<double, double> trainingData){
            //Put labels in column 0, rest of data in other columns.

        }
        vector<double> predict(vector<double> input){
            return activateLayer(0, input);
        }
    private:
        vector<double> activateLayer(int i, vector<double> inputs){
            if (i>= Layers.size()){
                return inputs;
            }
            else{
                activateLayer(i+1, Layers[i].feedForward(inputs));
            }
        }
        void backpropagate(vector<double> labels, vector<double> inputs){

        }
        
};