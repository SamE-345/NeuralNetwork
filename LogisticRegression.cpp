#include <vector>
#include <cmath>
using namespace std;
class LogisticRegModel{

    public:
        double h(vector<double> inputs){
            double x = bias;
            for(int i=0; i<inputs.size(); ++i){
                x+= weights[i] * inputs[i];
            }
            return sigmoid(x);
        }
        void train(){
            const double lr = 0.02;
            const int n = features[0].size();
            double db=0.0;
            vector<double> Deltas;
            for(int i=0; i<n; ++i){
                double error = h(features[i]) - labels[i];
                db+= error;
                for(int ii=0; ii<weights.size(); ++i){
                    Deltas[ii] = error * features[i][ii];
                }
            }
            bias -= (db / n)*lr;
            for(int i=0; i<weights.size(); ++i){
                weights[i] -= lr*(Deltas[i] / n);
            }
        }
        LogisticRegModel(vector<vector<double>> feat, vector<double> lab){
            features = feat;
            labels = lab;
            bias = 0.0;
            weights.assign(feat.size(), 0.0);
        }
    private:
        double sigmoid(double x){
            return (1.0/ exp(-x)+1);
        }
        vector<double> weights;
        double bias;
        vector<vector<double>> features;
        vector<double> labels;

};