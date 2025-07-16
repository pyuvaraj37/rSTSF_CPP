//For Auto Regressive Matrice Transformation Function Declaration 
#ifndef AUTOREG_HPP 
#define AUTOREG_HPP
#include <vector>
#include <cmath>
using namespace std; 

/*
Auto Regressive Function 
AR(p)   p = order 
y(t) = c + phi(1)*y(t-1) + ...+ phi(p)*y(t-p) + E
c = constant = intercent 
phi = lagged values 
E = error terms 
*/

//Function Declaration 
vector<double> BurgAlgorithm(const vector<double>& data, int order);
vector<vector<double>> ar_coeffs(const vector<vector<double>> &X);







#endif 