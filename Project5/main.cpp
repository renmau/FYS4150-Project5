#include <iostream>
#include <random>
#include <fstream>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;


// Writing money to file to make histogram:
void WriteFile(vec cash, int agents){
    ofstream myfile;
    myfile.open("e_gamma_1_alpha_1_lambda_0_1e5trans_1e3sim.txt");

    for (int i=0; i<agents;i++){
        myfile <<cash(i)<<endl;
    }
    myfile.close();
}


void Transactions_MC(mat transactions, vec& cash, int agents, int m_0, double lambda,double alpha, double gamma){
    int trans_number = 1e5; // at least 1e7

    //double prob;
    int i; int j;

    //RNG:
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> intRNG(0, agents-1);
    uniform_real_distribution<double> doubleRNG(0,1);
    ofstream myfile;
    myfile.open("test2.txt");


    for(int k=0; k<trans_number;k++){

        // Picking two random agents
        i = intRNG(gen);
        j = intRNG(gen);

        double avg_m = sum(cash)/agents;

        if (i!=j){

            int c_ij = transactions(i,j);
            double max_trans = transactions.max();
            double prob = 2*pow(fabs(cash(i) - cash(j))/avg_m, -alpha)*(pow((c_ij+1)/(max_trans+1), gamma));
            double test_nr = doubleRNG(gen);


            if (test_nr < prob){
                double eps = doubleRNG(gen);
                double sum_ij = cash(i) + cash(j);

                cash(i) = lambda*cash(i) + eps*(1.-lambda)*sum_ij;
                cash(j) = lambda*cash(j) + (1.-eps)*(1.-lambda)*sum_ij;

                transactions(i,j) += 1;
                transactions(j,i) += 1;
            }

            //double eps = doubleRNG(gen);

            // Just plain transaction:
            //double sum_ij = cash(i) + cash(j);
            //cash(i) = eps*sum_ij;
            //cash(j) = (1-eps)*sum_ij;

            // transaction + saving:
            //double delta_m = (1.0-lambda)*(eps*cash(j)-(1.0-eps)*cash(i));
            //cash(i) += delta_m;
            //cash(j) -= delta_m;
        }
    }

    cash = sort(cash);
    myfile.close();
}



int main()
{

    int agents = 500;      // Number of trading agents
    double m_0 = 100.0;      // initial money per agent
    double lambda = 0.;   // saving
    double alpha = 1.0;
    double gamma = 1.0;

    int simulations = 1e3;

    vec final_distribution = zeros<vec>(agents);


    // MC simulations:
    for(int k=0; k<simulations;k++){

        cout <<k<<endl;

        mat transactions = zeros<mat>(agents,agents);
        vec cash = zeros<vec>(agents)+m_0;
        Transactions_MC(transactions, cash, agents, m_0, lambda, alpha, gamma);

        final_distribution += cash;
    }


    final_distribution /= (simulations);
    cout << accu(final_distribution) << endl;

    WriteFile(final_distribution, agents);
}
