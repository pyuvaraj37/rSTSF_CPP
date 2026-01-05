//For Auto Regressive Matrice Transformation Function Declaration 
#ifndef AUTOREG_HPP 
#define AUTOREG_HPP
#include <vector>
#include <cmath>
using namespace std; 

//Function Declaration 
vector<double> BurgAlgorithm(const vector<double>& data, int order);
vector<vector<double>> ar_coeffs(const vector<vector<double>> &X);

#endif 