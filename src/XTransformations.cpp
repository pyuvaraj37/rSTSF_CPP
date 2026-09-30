//Use this to define autogregressive transformation 
#include "../include/XTransformations.hpp"


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
#include <complex>
#include <cmath>
#include <fftw3.h>
#include <stdexcept>

using namespace std;


//Burg for AR computation 
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

//AR Transformation
vector<vector<double>> ar_coeffs(const vector<vector<double>> &X){
    //Declare X_ar 
    vector<vector<double>> X_ar; 
    size_t num_columns = X[0].size();

    //Calculate Lags
    int lags = static_cast<int>(12.0 * pow(static_cast<double>(num_columns) / 100.0, 0.25));

    //Compute Burg on each row of X 
    for (const auto& row : X)
    {
        vector<double> coeffs = BurgAlgorithm(row, lags);

        //Flip sign convention to match aeon's AR representation
        for (auto& c : coeffs) {
            c = -c;
        }

        X_ar.push_back(coeffs);
    }

  

    return X_ar;
}


//Periodogram Transformation - TODO
vector<vector<double>> periodogram(const vector<vector<double>> &X){
if (X.empty()) return {};

    int n_samples = static_cast<int>(X.size());
    int n_feats   = static_cast<int>(X[0].size());
    int half      = n_feats / 2;

    std::vector<std::vector<double>> per_X(n_samples, std::vector<double>(half));

    // FFTW input/output buffers for a single row (real input -> complex output
    // is possible, but pyfftw.builders.fft assumes complex input/output by
    // default, so we replicate that exactly here).
    fftw_complex* in  = fftw_alloc_complex(n_feats);
    fftw_complex* out = fftw_alloc_complex(n_feats);

    // Plan once, reuse for every row (much faster than replanning each time)
    fftw_plan plan = fftw_plan_dft_1d(n_feats, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

for (int i = 0; i < n_samples; ++i) {
        // Load row into complex input (imaginary part = 0, since X is real)
        for (int j = 0; j < n_feats; ++j) {
            in[j][0] = X[i][j]; // real part
            in[j][1] = 0.0;     // imag part
        }

        fftw_execute(plan);

        // Magnitude of first half of the spectrum
        for (int j = 0; j < half; ++j) {
            double re = out[j][0];
            double im = out[j][1];
            per_X[i][j] = std::sqrt(re * re + im * im);
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return per_X;

}




//Difference Transformation
vector<vector<double>> difference(const vector<vector<double>> &X){
    //Declare X_diff 
    vector<vector<double>> X_diff; 
    size_t num_columns = X[0].size();

    //Compute Difference on each row of X 
    for (const auto& row : X)
    {
        vector<double> diff(num_columns - 1, 0.0);
        for (size_t n = 1; n < num_columns; ++n)
        {
            diff[n - 1] = row[n] - row[n - 1];
        }
        X_diff.push_back(diff);
    }

    return X_diff;
}