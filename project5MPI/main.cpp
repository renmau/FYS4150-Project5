#include <iostream>
#include <random>
#include <fstream>
#include <algorithm>
#include <mpi.h>
#include <armadillo>

using namespace std;
using namespace arma;

double sum_money(double* money,int nr_agents){
    double tot_money = 0;
    for (int i = 0; i < nr_agents; i++){
        tot_money += money[i];
    }
    return tot_money;
}

// Writing money to file to make histogram:
void WriteFile(string outfilename, double* cash, int agents){
    ofstream myfile;
    myfile.open(outfilename);

    for (int i=0; i<agents;i++){
        myfile <<cash[i]<<endl;
    }
    myfile.close();
}


void Transactions_MC(mat transactions, vec& cash, int agents, int m_0, double lambda, double alpha, double gamma){
    int trans_number = 1e5; // at least 1e7

    //double prob;
    int i; int j;

    //RNG:
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> intRNG(0, agents-1);
    uniform_real_distribution<double> doubleRNG(0,1);
    ofstream myfile;



    for(int k=0; k<trans_number;k++){

        // Picking two random agents
        i = intRNG(gen);
        j = intRNG(gen);

        double avg_m = sum(cash)/agents;
        //double max_trans = 0.;
        if (i!=j){
            int c_ij = transactions(i,j);
            double max_trans = (transactions).max();
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
        }

    }

    cash=sort(cash);
}



int main(int argn, char*argv[])
{
    MPI_Init(&argn, &argv);
    int numprocs, my_rank;
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    cout << "processor " << my_rank << "of " << numprocs << endl;

    /*int agents = 500;      // Number of trading agents
    double m_0 = 100.0;      // initial money per agent
    double lambda = 0.;   // saving
    double alpha = 1.0;
    double gamma = 1.0;*/

    string outfilename = argv[1];
    int agents = atoi(argv[2]);      // Number of trading agents
    double m_0 = 100.0;              // initial money per agent
    double lambda = atoi(argv[3]);   // saving
    double alpha = atoi(argv[4]);
    double gamma = atoi(argv[5]);
    int simulations = 1e3;


    vec final_distribution = zeros<vec>(agents);

    double* final_distribution_arr = new double[agents]; //convert arma::vec to c++ array after MC cycles are done
    double* total_distribution = new double[agents];
    for(int i=0; i< agents; i++) total_distribution[i] = 0.0;

    // MC simulations:
    for(int k=0; k<simulations;k++){

        cout <<k<<endl;
        mat transactions = zeros<mat>(agents,agents);
        vec cash = zeros<vec>(agents)+m_0;
        Transactions_MC(transactions, cash, agents, m_0, lambda, alpha, gamma);
        final_distribution += cash;
    }
    // Convert from armadillo vector to c++ array
    //for(int i=0; i< agents; i++) final_distribution_arr[i] = final_distribution(i);
    for(int i = 0; i < agents; i++){
        final_distribution_arr[i] = final_distribution(i);
        MPI_Reduce(&final_distribution_arr[i],&total_distribution[i] 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     }
    int my_rank=0;
    int numprocs = 1;
    if (my_rank==0){
        for (int i =0; i<agents;i++){
            total_distribution[i] /=((double) (numprocs*simulations));
            //final_distribution[i] /=((double) (simulations));
        }
        double tot_money = sum_money(total_distribution,agents);
        cout << "SUM: " << tot_money << endl;
        WriteFile(outfilename,total_distribution, agents);
    }
    MPI_Finalize();
}
