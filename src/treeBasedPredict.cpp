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

//---------HELPER FUNCTIONS----START--------

/**readTreesFromFile
 * @param       filePath   ,file to be read from 
 * @return      allTrees   ,matrix of Nodes 
 * 
 * **/
vector<vector<Node>> readTreesFromFile(const string& filePath) {
    ifstream infile(filePath);
    string line;
    vector<vector<Node>> forest;

    if (!infile) {
        cerr << "Error opening file\n";
        return forest;
    }

    int numTrees = 0;
    while (getline(infile, line)) {
        if (line.find("Number of Trees:") != string::npos) {
            numTrees = stoi(line.substr(line.find(":") + 1));
            forest.reserve(numTrees);
        } else if (line.find("Tree:") != string::npos) {
            vector<Node> tree;

            // Read "NumNodes: N" line
            getline(infile, line);
            int numNodes = stoi(line.substr(line.find(":") + 1));
            tree.reserve(numNodes);

            for (int i = 0; i < numNodes; ++i) {
                getline(infile, line);

                Node node;
                size_t pos;

                // Parse feature=
                pos = line.find("feature=");
                size_t comma = line.find(',', pos);
                node.feature = stoi(line.substr(pos + 8, comma - pos - 8));

                // Parse threshold=
                pos = line.find("threshold=", comma);
                comma = line.find(',', pos);
                node.threshold = stod(line.substr(pos + 10, comma - pos - 10));

                // Parse left=
                pos = line.find("left=", comma);
                comma = line.find(',', pos);
                node.left = stoi(line.substr(pos + 5, comma - pos - 5));

                // Parse right=
                pos = line.find("right=", comma);
                comma = line.find(',', pos);
                node.right = stoi(line.substr(pos + 6, comma - pos - 6));

                // Parse values=[
                pos = line.find("values=[", comma);
                size_t endBracket = line.find(']', pos);
                string valuesStr = line.substr(pos + 8, endBracket - pos - 8);

                stringstream ss(valuesStr);
                string val;
                while (getline(ss, val, ',')) {
                    node.values.push_back(stod(val));
                }

                tree.push_back(node);
            }

            forest.push_back(tree);
        }
    }

    return forest;
}

/**check_is_fitted
 * @param   forest 
 * @return  none
 * checks if the forest is valid, throws an exception if not...
 * **/
void checkIsFitted(const vector<vector<Node>>& forest) {
    if (forest.empty()) {
        throw runtime_error("Forest is empty. Make sure trees are loaded.");
    }
}


/**validateXPredict (complete ai)
 * @param   X, original matrix 
 * @return  vector<vector<double>> X, validated matrix 
 * **/
vector<vector<double>> validateXPredict(const vector<vector<double>>& X, bool allow_nan = false) {
    // 1. Check if model is fitted
    // (Assuming you have a boolean flag like `is_fitted` in your class)
    // if (!is_fitted) throw runtime_error("Model not fitted yet.");

    // 2. Check if X contains any NaN or infinite values (if allow_nan == false)
    for (const auto& row : X) {
        for (double val : row) {
            if (!allow_nan && !isfinite(val)) {
                throw runtime_error("Input contains NaN or infinite values.");
            }
        }
    }

    // 3. Optionally check sparse matrix indices here (skip if not using sparse matrices)
    // ...

    // 4. Return validated X (you could also copy or preprocess here if needed)
    return X;
}


/**getTreeProba
 * @param       forest, X 
 * @return      matrice of probabilities
 * Traverses each tree of each forest with corresponding row of X 
 * **/
vector<vector<vector<double>>> getTreeProba(const vector<vector<Node>>& forest, 
                                            const vector<vector<double>>& X, 
                                            vector<vector<vector<double>>>& all_proba, 
                                            const int& numSamples, 
                                            const int& numClasses, 
                                            const int& numOut){
        for (size_t i=0; i<numOut; i++){ //Output i 
            for (size_t j = 0; j < 5; j++) {  // Tree j
                const auto& tree = forest[j];
                for (size_t k=0; k<numSamples; k++){   //Sample k 
                    //Get row of X 
                    vector<double> XRow = X[k];

                    //TREE TRAVERSAL STEP 
                    int nodeIdx = 0;           //Start at root node 
                    while (tree[nodeIdx].left != -1 && tree[nodeIdx].right != -1){
                        const Node& node = tree[nodeIdx]; 
                        double xValue = XRow[node.feature];
                        cout << "Tree " << j 
                            << " Sample " << k 
                            << " Node " << nodeIdx
                            << " Feature " << node.feature 
                            << " Threshold " << node.threshold 
                            << " Value " << xValue << endl;
                        if(xValue <= node.threshold){
                            nodeIdx = node.left; 
                        }else{
                            nodeIdx = node.right; 
                        }
                    }

                    //Once a leaf node is hit, get the values 
                    const vector<double>& probs = tree[nodeIdx].values;
                    cout << "Hit leaf at: " << nodeIdx << endl; 
                    cout << "Values to be added: [ "; 
                    for (int p=0; p<numClasses; p++){
                        cout << probs[p] << " " << endl; 
                    }
                    cout << "]" << endl; 

                    //Summing the probabilities and adding to all_proba 
                    for (size_t c = 0; c < probs.size(); ++c) {
                        cout << "Before Adding...Class " << c << " = " << all_proba[i][k][c];
                        all_proba[i][k][c] += probs[c];
                        cout << "After Adding " << probs[c] << "...Class " << c << " = " << all_proba[i][k][c] << endl; 
                    }
                }  
            }
        }
       
    return all_proba;
}


/**predictProba (for 1D Array)
 * @param   X, a matrix of doubles (xIntTrans)
 * @return  matrice of doubles 
 * Here is where we are getting the probabilities of the forest 
 */
vector<vector<vector<double>>> predictProba(const vector<vector<double>>& X_input, 
                                            const vector<vector<Node>>& forest, 
                                            const vector<int>& classes, 
                                            const int& numOut)
{

    //Prestep: Get required variables 
    int numSamples = X_input.size(); 
    int numClasses = classes.size(); 
   
    // //1. VALIDATE TRANSFORMED DATA 
    checkIsFitted(forest); 
    vector<vector<double>> X = validateXPredict(X_input);

    //3. ALLOCATE RESULT LIST OF MATRICES w/zeroes  
    int XRows = X.size(); 
    int XColumns = X[0].size(); 
    int selfNClasses = 2; 
    vector<vector<vector<double>>> allProba(numOut, vector<vector<double>>(numSamples, vector<double>(numClasses, 0.0)));

    //4. TRAVERSING THE TREES TO GET PROBABILITIES  
    allProba = getTreeProba(forest, X, allProba, numSamples, numClasses, numOut); 

    //5. AVERAGING THE PROBABILITIES
    int numTrees = forest.size();
    for (int k = 0; k < numSamples; ++k) {
        for (int i = 0; i < numOut; ++i) {
            for (int c = 0; c < numClasses; ++c) {
            allProba[i][k][c] /= static_cast<double>(numTrees);  
            }
        }
    }

    return allProba;
}
//---------HELPER FUNCTIONS  ----- END -----






//MAIN WRAPPER FUNCTION FOR 
/**treeBasedPredict
 *
 * The C++ equivalent of 
 * 
 * @param       X, a matrix of doubles (XIntTrans in intBasedTrans)
 * @return      
 * 
 * */

vector<int> treeBasedPredict(const vector<vector<double>>& X){

    //PRESTEP: GET FORESTS AND OTHER VARS 
    vector<vector<Node>> forest = readTreesFromFile("/home/ccuev029/DATA/extraTrees.txt"); 
    vector<int> classes = {0,1};    //Caution 
    int numOut = 1;                 //Caution 
    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Printing the Forest 
    cout << "Number of Trees: " << forest.size() << endl; 
    for (size_t t = 0; t < forest.size(); ++t) {
        cout << "Tree " << t << ":\n";
        for (size_t n = 0; n < forest[t].size(); ++n) {
            const Node& node = forest[t][n];
            cout << "  Node " << n << ": ";
            cout << "feature=" << node.feature << ", ";
            cout << "threshold=" << node.threshold << ", ";
            cout << "left=" << node.left << ", ";
            cout << "right=" << node.right << ", ";
            cout << "values=[";
            for (size_t v = 0; v < node.values.size(); ++v) {
                cout << node.values[v];
                if (v + 1 < node.values.size()) cout << ", ";
            }
            cout << "]\n";
        }
        cout << endl;
    }

    //3. GET THE PROBABILITIES OF TRANSFORMED MATRIX (Using the forest)
    vector<vector<vector<double>>> proba = predictProba(X, forest, classes, numOut); 
    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Printing Proba  
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

