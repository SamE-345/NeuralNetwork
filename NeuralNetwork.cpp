#include <vector>
#include <cstdlib>
#include <cmath>
using namespace std;

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
        double Sigmoid(double x){
            return (1/(1+exp(-x)));
        }
        double Activate(vector<double> inputs){
            double sum = bias;
            for(int i=0;i< weights.size(); i++ ){
                sum += weights.at(i) * inputs[i];
            }
            return sum;
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
            //
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
        
};