
#ifndef TREE_BASED_PREDICT_HPP
#define TREE_BASED_PREDICT_HPP

#include<iostream>
#include<vector>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <thread> //For Parallel 
#include <mutex> //For better organization, working in Parallel
using namespace std; 

//Function Prototypes 
vector<vector<Node>> readTreesFromFile(const string& filePath);
pair<int, vector<int>> partitionEstimators(int n_estimators, int numJobsRequested);
void accumulatePrediction(const vector<vector<double>>& X,
                          const vector<vector<Node>>& forestPart,
                          vector<vector<double>>& all_proba,
                          mutex& lock);
vector<vector<double>> predictProba(const vector<vector<double>>& X,
                                    const vector<vector<Node>>& forest);
vector<vector<double>> treeBasedPredict(const vector<vector<double>>& X);


#endif // TREE_BASED_PREDICT_HPP
