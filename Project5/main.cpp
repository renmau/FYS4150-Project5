#include <iostream>
#include <random>
#include <fstream>
#include <algorithm>

using namespace std;

// Array containing money of all agents:
double* CashArray(int agents, int m_0){

    double* cash = new double[agents];

    for (int i=0; i<agents; i++){
        cash[i] = m_0;
    }

    return cash;
}

// Matrix containing all transactions:
int** TransactionMatrix(int agents){
    int** transactions = new int*[agents];

    for (int i=0; i<agents;i++){
        transactions[i] = new int[agents];
    }

    for (int i=0; i<agents;i++){
        for (int j=0; j<agents;j++){
            transactions[i][j] = 0;
        }
    }

    return transactions;
}

// total sum of all money:
void sum_cash(double* cash, int agents){
    double sum_cash = 0;

    for (int i=0; i<agents;i++){
        sum_cash += cash[i];
    }

    cout << 'Sum of all money:'<<sum_cash<<endl;
}

// Writing money to file to make histogram:
void WriteFile(double* cash, int agents){
    ofstream myfile;
    myfile.open('test.txt');

    for (int i=0; i<agents;i++){
        myfile <<cash[i]<<endl;
    }
    myfile.close();
}


double* Transactions_MC(int** transactions, double* cash, int agents, int m_0){
    int trans_number = 1e4; // at least 1e7
    int simcounter = 0;
    int transcounter = 0;

    double prob;
    int i; int j;

    //RNG:
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> intRNG(0, agents-1);
    uniform_real_distribution<double> doubleRNG(0,1);
    ofstream myfile;
    myfile.open("test2.txt");

    while (transcounter < trans_number){
        // Picking two random agents
        i = intRNG(gen);
        j = intRNG(gen);
        if (i!=j){

            int c_ij = transactions[i][j];
            //prob = pow(fabs(cash[i] - cash[j]), -alpha)*(pow(c_i_j+1, gamma));
            double test_nr = doubleRNG(gen);

            /*
            if (test_nr < prob){
                double eps = doubleRNG(gen);
                double sum_ij = cash[i] + cash[j];

                cash[i] = lambda*cash[i] + eps*(1.-lambda)*sum_ij;
                cash[j] = lambda*cash[j] + (1.-eps)*(1.-lambda)*sum_ij;
                transcounter += 1;
                transactions[i][j] += 1;
                sum_cash(cash, agents);
            }*/

            double epsilon = doubleRNG(gen);
            double sum_ij = cash[i] + cash[j];

            cash[i] = eps*sum_ij;
            cash[j] = (1-eps)*sum_ij;


        }

        /*
        double sum_variance = 0;
        for (int i=0; i<agents;i++){
            sum_variance += (cash-[i] - m_0)*(cash[i]-m_0);

        }
        double new_variance = sum_variance/agents;
        myfile << new_variance << endl;
        cout << sum_variance << endl;
        if (new_variance > m_0*m_0/2){
            cout << i << endl;
            cout << old_variance << "   " << new_variance << endl;
            break;

        }
        if (transcounter > 2){
            double old_variance = sum_variance/agents;
        }*/

    }


    simcounter += 1;
    sort(cash, cash + agents);
    myfile.close();
    return cash;
}



int main()
{
    int agents = 500;       // Number of trading agents
    double m_0 = 2;         // initial money

    int simulations = 1e4;
    int simulations_counter = 0;

    double* final_distribution = CashArray(agents, m_0);


    // transactions:
    while (simulations_counter < simulations){

        int** transactions = TransactionMatrix(agents);

        double* cash = CashArray(agents, m_0);
        sum_cash(cash,agents);

        double* new_cash = Transactions_MC(transactions, cash, agents, m_0);

        for (int i=0; i<agents; i++){
            final_distribution[i] += new_cash[i];
        }

        simulations_counter += 1;
    }

    for (int i=0; i<agents; i++){
        final_distribution[i] /= simulations;
    }

    WriteFile(final_distribution, agents);
}
