
#ifndef TREE_BASED_PREDICT_HPP
#define TREE_BASED_PREDICT_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <thread> //For Parallel 
#include <mutex>  //For memory locks
#include <stdexcept>    //For throwing errors 
using namespace std; 


struct Node {
    int feature;            //Index of the input data (the row of X to compare to)
    float threshold;       //Value of comparison for split (left or right)
    int left;               //Index of the left child          
    int right;              //Index of the right child 
    vector<float> values;   //class counts or regression values for each node
};



//Function Prototypes 
// vector<vector<Node>> readTreesFromFile(const string& filePath);
vector<vector<vector<float>>> getTreeProba(const vector<vector<Node>>& forest, 
                                            const vector<vector<float>>& X, 
                                            vector<vector<vector<float>>>& all_proba, 
                                            const int& numSamples, 
                                            const int& numClasses, 
                                            const int& numOut); 
vector<vector<vector<float>>> predictProba(const vector<vector<float>>& X,
                                    const vector<vector<Node>>& forest,
                                    const vector<int>& classes, 
                                    const int& numOut);

vector<int> treeBasedPredict(const vector<vector<float>>& X, vector<vector<Node>>& forest);


#endif // TREE_BASED_PREDICT_HPP
