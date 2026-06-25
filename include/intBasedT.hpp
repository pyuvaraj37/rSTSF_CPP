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
vector<vector<float>> getIntervalBasedTransform (vector<vector<float>> X, 
                                                vector<vector<float>> X_ar, 
                                                vector<vector<float>> X_per, 
                                                vector<vector<float>> X_diff, 
                                                vector<vector<vector<float>>> all_caf
                                                ); 

#endif