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
            for(int i=0; i<sizeof(features); i++){
                sum += (H(features[i]) - labels[i])*(H(features[i]) - labels[i]);
            }
            return sum*(1/2*m);
        }
        double dJ(){
            double sum=0;
            for(int i=0; i<sizeof(features); i++){
                sum += (H(features[i]) - labels[i])*features[i];
            }
            return sum;
        }
    public:
        SingleVarLinearReg(vector<double> feat,vector<double> lab ){
            features = feat;
            labels = lab;
            m = sizeof(lab);

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
    //TODO
};