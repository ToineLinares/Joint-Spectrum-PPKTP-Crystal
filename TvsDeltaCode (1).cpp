#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int main() //The goal is to analyze the combinations of lambda pump, delta lambda and temperature of the crystal at which the signal and idler photons have parallel wave vectors.
    
{
    //Define constants relevant to the experiment and the PPKTP crystal
    double A_0 = 3.42445;
    double AA = 2.12725;
    double BB = 1.18431;
    double CC = 5.14852e-2;
    double DD = 0.6603;
    double EE = 100.00507;
    double FF = 9.68956e-3;
    double t_0 = 24.09;
    double lambdas[101];
    double T[1001];
    double deltas[1000];
    double w_p;
    double w_s;
    double w_i;
    double lam_s;
    double lam_i;
    double temps[1000];
    double pi = 3.14159;
    double c = 2.998e14;
    double periodo;
    double n_s;
    double n_i;
    double n_p;
    double n_1s;
    double n_1i;
    double n_2s;
    double n_2i;
    double n_1p;
    double n_2p;
    double dnz;
    double nzs;
    double nzi;
    double nzp;
    double theta_s;
    int count;
    ofstream outfile_1;
    ofstream outfile_2;
    ofstream outfile_3;
    ofstream outfile_4;
    double DELTA;
    for(int i = 0; i <= 100; i++)
    {
        lambdas[i] = 0.00002*i + 0.404; //generates a list with all possible pump wavelengths between 404 and 406 nm
    }
    for(int i = 0; i <= 1000; i++)
    {
        T[i] = 0.1*i; //generates a list with the temperatures at which I wish to evaluate the behavior of the angle of emission of the signal photon
    }

    outfile_1.open("LDTpsi.txt");

    for (int j = 0; j <= 100; j++) //First loop of a triple loop: Iterates over all possible values of lambda pump.
    {
        w_p = (2.0*pi*c)/(lambdas[j]); 
        for (int s = -500; s < 500; s++)  //Generates a list with all possible differences in frequency that may give a collinear configuration
        {
            deltas[s+500] = s*0.0008*w_p;
        }
        for (int r = 0; r < 1000; r++)//Second loop of a triple loop: iterates over all possible values of difference in signal and idler frequency.
        {
            count = 0;
            DELTA = 0.4*w_p;
            
            for (int k = 0; k <= 1000; k++) //Third loop of a triple loop: iterates over all possible values of temperature at which the crystal may be set.
            {
                
                w_s = (w_p/2) + deltas[r];
                w_i = (w_p/2) - deltas[r];
                lam_s = 2.0*pi*c/w_s;
                lam_i = 2.0*pi*c/w_i;
                periodo = A_0*(1 + 0.00000077*(T[k] - 25.0) + 0.0000000071*( (T[k] - 25.0)*(T[k] - 25.0) ));
                
                //Calculation of the index of refraction according to Sellmeier equations
                
                n_1s = 0.0000099587 + (0.0000099228/lam_s) + ((-0.0000089603)/((lam_s)*(lam_s))) + ((0.0000041010)/((lam_s)*(lam_s)*(lam_s)));
                n_1i = 0.0000099587 + (0.0000099228/lam_i) + ((-0.0000089603)/((lam_i)*(lam_i))) + ((0.0000041010)/((lam_i)*(lam_i)*(lam_i)));
                n_1p = 0.0000099587 + (0.0000099228/lambdas[j]) + ((-0.0000089603)/((lambdas[j])*(lambdas[j]))) + ((0.0000041010)/((lambdas[j])*(lambdas[j])*(lambdas[j])));
                
                n_2s = (0.00000010459/lam_s) + ((-0.000000098136)/((lam_s)*(lam_s))) + ((0.000000031481)/((lam_s)*(lam_s)*(lam_s))) - 0.000000011882;
                n_2i = (0.00000010459/lam_i) + ((-0.000000098136)/((lam_i)*(lam_i))) + ((0.000000031481)/((lam_i)*(lam_i)*(lam_i))) - 0.000000011882; 
                n_2p = (0.00000010459/lambdas[j]) + ((-0.000000098136)/((lambdas[j])*(lambdas[j]))) + ((0.000000031481)/((lambdas[j])*(lambdas[j])*(lambdas[j]))) - 0.000000011882; 
                
                nzs = sqrt( AA + ((BB)/(1 - ( (CC)/(lam_s*lam_s) )  )) +  ((DD)/(1 - ( (EE)/(lam_s*lam_s) )  )) - FF*lam_s*lam_s ) + (n_1s*(T[k] - t_0)) + (n_2s*(T[k] - t_0)*(T[k] - t_0));
                
                nzi = sqrt( AA + ((BB)/(1 - ( (CC)/(lam_i*lam_i) )  )) +  ((DD)/(1 - ( (EE)/(lam_i*lam_i) )  )) - FF*lam_i*lam_i ) + (n_1i*(T[k] - t_0)) + (n_2i*(T[k] - t_0)*(T[k] - t_0));
                
                nzp = sqrt( AA + ((BB)/(1 - ( (CC)/(lambdas[j]*lambdas[j]) )  )) +  ((DD)/(1 - ( (EE)/(lambdas[j]*lambdas[j]) )  )) - FF*lambdas[j]*lambdas[j] ) + (n_1p*(T[k] - t_0)) + (n_2p*(T[k] - t_0)*(T[k] - t_0));
                
                theta_s = sqrt( ((2.0*nzi*w_i)/(nzs*w_s))*(1 + (   ((2.0*pi*c/periodo) - (nzp*w_p))/(nzs*w_s + nzi*w_i)) )); //Calculate angle of emission of signal photon
                if(theta_s <= 0.001) //Evaluate is theta_s is approximately zero
                {
                    outfile_1 << lambdas[j] << " " << deltas[r] << " " << T[k] << " " << w_p << " " << w_s << " " << w_i << endl; //Write the combination in a text file to analyze in python
                }
                
            }
        }
        
    }
    outfile_1.close();

    
}            
