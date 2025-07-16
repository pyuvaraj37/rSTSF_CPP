//Use this to define autogregressive transformation 
#include "../include/autoreg.hpp"


/* In Python..
#𝓒.𝓙𝓪𝓷𝓮: Autoregressive Representation Conversion 
def ar_coefs(X):
    X_transform = []                                #Transformed Matrix Declaration         
    lags = int(12*(X.shape[1]/100.)**(1/4.))        #Generates lag 
    for i in range(X.shape[0]):                     #Goes through each element in the vecotr 
        coefs,_ = burg(X[i,:],order=lags)           #Coefficients 
        X_transform.append(coefs)                   #Adds to the transformed matrix 
    return np.array(X_transform)                    #returns the transformed matrix 
*/

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <vector>
using namespace std;

/** 
 * A test harness for Cedrick Collomb's Burg algorithm variant.
 *
 * Taken from Cedrick Collomb. "Burg's method, algorithm, and recursion",
 * November 2009 available at http://www.emptyloop.com/technotes/.
 */

/**
 * Returns in vector coefficients calculated using Burg algorithm applied to
 * the input source data x
 */


//Burg for AR computation 
#include <vector>
#include <cmath>
using namespace std;

vector<double> BurgAlgorithm(const vector<double>& x, int m)
{
    size_t N = x.size() - 1;

    vector<double> Ak(m + 1, 0.0);
    Ak[0] = 1.0;

    vector<double> f(x);
    vector<double> b(x);

    double Dk = 0.0;
    for (size_t j = 0; j <= N; j++)
    {
        Dk += 2.0 * f[j] * f[j];
    }
    Dk -= f[0] * f[0] + b[N] * b[N];

    for (size_t k = 0; k < m; k++)
    {
        double mu = 0.0;
        for (size_t n = 0; n <= N - k - 1; n++)
        {
            mu += f[n + k + 1] * b[n];
        }
        mu *= -2.0 / Dk;

        for (size_t n = 0; n <= (k + 1) / 2; n++)
        {
            double t1 = Ak[n] + mu * Ak[k + 1 - n];
            double t2 = Ak[k + 1 - n] + mu * Ak[n];
            Ak[n] = t1;
            Ak[k + 1 - n] = t2;
        }

        for (size_t n = 0; n <= N - k - 1; n++)
        {
            double t1 = f[n + k + 1] + mu * b[n];
            double t2 = b[n] + mu * f[n + k + 1];
            f[n + k + 1] = t1;
            b[n] = t2;
        }

        Dk = (1.0 - mu * mu) * Dk - f[k + 1] * f[k + 1] - b[N - k - 1] * b[N - k - 1];
    }

    // Return AR coefficients excluding Ak[0] (which is always 1)
    return vector<double>(Ak.begin() + 1, Ak.end());
}


/* AR Transformation Function Declaration 
@param      Matrix "X" to be transformed
@returns    Corresponding AR Representation 
*/
vector<vector<double>> ar_coeffs(const vector<vector<double>> &X){
    vector<vector<double>> X_ar; 
    size_t num_columns = X[0].size();
    int lags = static_cast<int>(12.0 * pow(static_cast<double>(num_columns) / 100.0, 0.25));
    cout << "C++ lags: " << lags << endl; 
    for (const auto& row : X)
    {
        vector<double> coeffs = BurgAlgorithm(row, lags);
        X_ar.push_back(coeffs);
    }
    return X_ar;
}