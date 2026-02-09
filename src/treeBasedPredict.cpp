/**C++ Equivalent of self.extra_trees.predict() in rSTSF
 * extra_trees: declared in fit() of rstsf() class
 *              set to ExtraTreesClassifier(self.extra_trees = ExtraTreesClassifier(n_estimators=self.r,criterion='entropy',class_weight='balanced',max_features='sqrt')
 *              data was trained with fit(all_X_train_T, y_train)
 * Will be using self.extra_trees directly from python code (take from "self.extra_trees.txt")
 * Inference step (self.extra_trees.fit()) is not replicated..
 * */

 /** Explanation: self.extra_trees representation 
  * self.extra_trees is a list of individal decision TREE OBJECTS 
  * Tree Objects    :   each tree holds a certain number of NODES 
  * Nodes           :   each node has five features, all saved as type double 
  * 
  * struct Node is declared to handle these features 
  * vector<vector<Node>> Creates the matrix representation 
  * #rows = number of 
  * **/

#include "treeBasedPredict.hpp"
#include "extraFunctions.h"



vector<int> treeBasedPredict(const vector<vector<double>>& X){

    //PRESTEP: GET FORESTS AND OTHER VARS 
    vector<vector<Node>> forest = readTreesFromFile("/home/ccuev029/DATA/extraTrees.txt"); 
    vector<int> classes = {0,1};    //Caution 
    int numOut = 1;                 //Caution 
    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Printing the Forest 
    printForest(forest); 
   

    //3. GET THE PROBABILITIES OF TRANSFORMED MATRIX (Using the forest)
    vector<vector<vector<double>>> proba = predictProba(X, forest, classes, numOut); 
    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: 
    //Printing Proba  
    if (proba.empty()) {
    cerr << "ERROR: proba is empty!" << endl;
    }
    cout << "proba shape: " << proba.size() << " " << proba[0].size() << endl; 
    cout << "proba: " << endl; 
    for(size_t i=0; i<proba.size(); i++){
    cout << "["; 
        for(size_t j=0; j<proba[i].size(); j++){
        cout << "["; 
            for(size_t k=0; k<proba[i][j].size(); k++){
                cout << proba[i][j][k]; 
                if (k < classes.size() - 1)
                std::cout << ", ";

            }
        cout << "]" << endl; 
        }
    cout << "]" << endl; 
    }
    //Compare to original
    

    vector<int> predictions;
    //3. GET THE LARGEST PROBABILITIES //caution
    for(size_t i=0; i<proba.size(); i++){//For first (and only output)
            for(size_t j=0; j<proba[i].size(); j++){//Probability set of sample j 
                const vector<double>& classProbs = proba[i][j]; //Get the row of probabilities for that sample
                double maxValue = -1;  
                int classification = 0; 
                for(size_t k=0; k<proba[i][j].size(); k++){//Probability k 
                    if (classProbs[k]>maxValue){
                        maxValue = classProbs[k];
                        classification = classes[k]; 
                    }

                }
            predictions.push_back(classification);
            }
    }
    return predictions; 
}