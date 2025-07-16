//Put all header files for intervalBased Transfromation here 
#ifndef INTBASEDT_HPP 
#define INTBASEDT_HPP
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
using namespace std; 




//Get Interval Based Transform 
vector<vector<double>> getIntervalBasedTransform (vector<vector<double>> X, 
                                                vector<vector<double>> X_ar, 
                                                vector<vector<double>> X_per, 
                                                vector<vector<double>> X_diff, 
                                                vector<vector<double>> allCaf, //maybe change this later 
                                                vector<int> relevantCaf); 

#endif