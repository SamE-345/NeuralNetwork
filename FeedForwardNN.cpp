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
        vector<double> feedForward(vector<double> inputs, bool training){
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


class FeedForwardNN{
    public:
        vector<Layer> Layers;
        FeedForwardNN(vector<int> layerSizes){
            int x = 0;
            for(int i=0; i<layerSizes.size(); i++){
                Layers.push_back(Layer(layerSizes.at(i),x));
                x=layerSizes.at(i);
            }
        }
        void trainModel(vector<vector<double>> trainingData, int epochs, vector<vector<double>> labels, const double lr = 0.02){
            for (int e = 0; e < epochs;e++) {
                for (int i = 0; i < trainingData.size(); ++i) {
                    backpropagate(trainingData[i], labels[i], lr);
                }
            }

        }
        vector<double> predict(vector<double> input){
            return activateLayer(0, input);
        }
        void backpropagate(vector<double> labels, vector<double> inputs, const int learningRate){
            ActivationFunction funct;
            vector<double> L0 = inputs;
            vector<vector<double>> activations = {L0}; //The activation for the layer1 neurons are the inputs
            for (int i = 0; i < Layers.size(); i++) {
                L0 = Layers[i].feedForward(L0, true); 
                activations.push_back(L0);
            }
            vector<vector<double>> DeltasList (Layers.size()); //List of deltas for each layer
            vector<double> delta; // List of deltas for one layer
            int end = Layers.size();
            for(int i=0; i< Layers[end].layer.size(); i++){ // Calculates error for output layer
                double error = funct.sigmoid(Layers[end].layer[i].lastSum) - labels[i] ;
                delta.push_back(error * funct.sigmoidDerivative(Layers[end].layer[i].lastSum));
            }   
            DeltasList[end] = delta; // Final layer has deltas calculated above


            for(int iii = end-1; iii>= 0; iii--){ // Calculate error in hidden layers
                vector<double> hiddenDelta;
                for(int iv=0; iv<Layers[iii].layer.size(); iv++){ //Neuron
                    double error = 0.0;
                    for(int v=0; v<Layers[iii+1].layer.size(); v++){ //Weights
                        error += DeltasList[iii+1][v] * Layers[iii+1].layer[v].weights[iv];
                    }
                    hiddenDelta.push_back(error * funct.sigmoidDerivative(Layers[iii].layer[iv].lastSum));

                }
                DeltasList[iii] = hiddenDelta;
            }

            //Update weights next
            for(int i=0; i<Layers.size(); i++){ //Iterate all layers
                for(int ii=0; ii<Layers[i].layer.size(); i++){ //Iterate all neurons
                    for(int iii=0; iii<Layers[i].layer[ii].weights.size(); iii++){ //Iterate all weights
                        Layers[i].layer[ii].weights[iii] -= DeltasList[i][ii] * learningRate * Layers[i].layer[ii].lastInputs[iii];
                    }
                    Layers[i].layer[ii].bias -= learningRate * DeltasList[i][ii];
                }
            }
        }


    private:
        vector<double> activateLayer(int i, vector<double> inputs){
            if (i>= Layers.size()){
                return inputs;
            }
            else{
                return activateLayer(i+1, Layers[i].feedForward(inputs, false));
            }
        }
        double meanSqrError(vector<double> labels, vector<double> inputs){
            double loss = 0;
            for(int i=0; i<inputs.size();i++){
                loss += pow((inputs[i]-labels[i]),2);
            }
            return loss/inputs.size();
        }
};