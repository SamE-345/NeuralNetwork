#include <vector>
using namespace std;
class SingleVarLinearReg{

    private:
        vector<double> features;
        vector<double> labels;
        double w = 0.0;
        double b = 0.0;
        int m;
        double J(){
            double sum=0;
            for(int i=0; i<features.size(); i++){
                sum += (H(features[i]) - labels[i])*(H(features[i]) - labels[i]);
            }
            return sum*(1/2.0*m);
        }
        double dJ(){
            double sum=0;
            for(int i=0; i<features.size(); i++){
                sum += (H(features[i]) - labels[i])*features[i];
            }
            return sum;
        }
    public:
        SingleVarLinearReg(vector<double> feat,vector<double> lab ){
            features = feat;
            labels = lab;
            m = lab.size();
        }
        double H(double x){
            return (w*x +b);
        }
        void train(int epochs, double lr){
            for(int epoch=0; epoch<epochs; ++epoch){
                double dw = 0.0; double db = 0.0;
                for(int i=0; i<m; ++i){
                    double predi = H(features[i]);
                    double error = predi - labels[i];
                    dw += error * features[i];
                    db += error;
                }
                dw /= m;
                db /= m;
                w -= lr* dw;
                b -= lr*db;
            }
        }
};
class MultiVarLinearReg{
    private:
        vector<double> weights;
        vector<vector<double>> features;
        vector<double> labels;
        double bias = 0.0;
        double J(int index){
            return (h(features[index]) - labels[index]) * (h(features[index]) - labels[index]);
        }
        vector<double> dJ(int index){
            vector<double> delta(weights.size());
            double error = h(features[index]) - labels[index];
            for(int i=0; i< weights.size(); ++i){
                delta[i] = error * 2 * features[index][i];
            }
            bias += error;
            return delta;
        }

    public:
        MultiVarLinearReg(vector<vector<double>> feat, vector<double> lab){
            features = feat;
            for(int i=0; i<features[0].size(); ++i){
                weights.push_back(0.0);
            }
            features = feat;
            labels = lab;
        }
        double h(vector<double> theta){
            double sum =0;
            for(int i=0; i<weights.size(); ++i){
                sum += weights[i] * theta[i];
            }
            return sum;
        }
        void train(){
            const double lr = 0.02;
            const int n = features[0].size();
            vector<double> Deltas;
            for(int i=0; i<n; ++i){
                Deltas = dJ(i);
                for(int ii=0; ii<weights.size();++ii){
                    weights[ii] -= lr*Deltas[ii];
                }
            }

        }
};