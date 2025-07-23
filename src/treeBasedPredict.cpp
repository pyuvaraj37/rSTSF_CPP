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

#include "treeBasedPredict.hpp"; 

struct Node {
    int feature;            //Index of the input data (the row of X to compare to)
    double threshold;       //Value of comparison for split (left or right)
    int left;               //Index of the left child          
    int right;              //Index of the right child 
    vector<double> values;   //class counts or regression values for each node
};



//---------Other functions----start--------

/**readTreesFromFile
 * @param       filePath   ,file to be read from 
 * @return      allTrees   ,matrix of Nodes 
 * 
 * **/
vector<vector<Node>> readTreesFromFile(const string& filePath) {
    ifstream infile(filePath);
    string line;
    vector<vector<Node>> allTrees;
    vector<Node> currentTree;

    while (getline(infile, line)) {
        // Skip header/info lines
        if (line.empty() ||
            line.find("===") != string::npos ||
            line.find("NumNodes:") != string::npos ||
            line.find("Number of Trees:") != string::npos ||
            line.find("For Each Node:") != string::npos) {
            
            if (!currentTree.empty()) {
                allTrees.push_back(currentTree);
                currentTree.clear();
            }
            continue;
        }

        istringstream iss(line);
        Node node;

        // Read feature, threshold, left, right
        iss >> node.feature >> node.threshold >> node.left >> node.right;

        // Read remaining part of the line (class distribution)
        string valuesStr;
        getline(iss, valuesStr);  // everything after right

        // Clean whitespace and commas
        valuesStr.erase(0, valuesStr.find_first_not_of(" ,"));
        valuesStr.erase(valuesStr.find_last_not_of(" ,") + 1);

        stringstream vss(valuesStr);
        string val;
        while (getline(vss, val, ',')) {
            stringstream valStream(val);
            double num;
            if (valStream >> num) {
                node.values.push_back(num);
            }
        }

        currentTree.push_back(node);
    }

    // Append final tree if any
    if (!currentTree.empty()) {
        allTrees.push_back(currentTree);
    }

    return allTrees;
}


/**partitionEstimators 
 * Is used to distribute the number of trees amongst multiple "jobs"/CPU Cores
 * @param      number of estimators, number of jobs requested 
 * @return     a pair (number of jobs, how many trees per job)
 * **/
pair<int, vector<int>> partitionEstimators(int n_estimators, int numJobsRequested) {
    //Get number of CPU Cores/Threads for Job distribution 
    int maxThreads = thread::hardware_concurrency();
    if (maxThreads == 0) {maxThreads = 1;}  // fallback if undetectable


    int numJobs = numJobsRequested;
    if (numJobs <= 0 || numJobs > maxThreads) {
        numJobs = maxThreads;
    }

    if (numJobs > n_estimators) {
        numJobs = n_estimators;
    }

    vector<int> estimators_per_job(numJobs, n_estimators / numJobs);
    int remainder = n_estimators % numJobs;
    for (int i = 0; i < remainder; ++i) {
        estimators_per_job[i]++;
    }

    return {numJobs, estimators_per_job};
}


/**accumulatePrediction
 * @param       X, the input data 
 *              forestPart, a subset of the decision trees 
 *              all_proba, to store the generated predictions 
 *              lock, a mutex lock to make sure updates are safe 
 * **/
void accumulatePrediction(const vector<vector<double>>& X,
                            const vector<vector<Node>>& forestPart, 
                            const vector<vector<double>>& all_proba, 
                            mutex& lock){
    //For each tree in the given part of the forest...
    for(const auto&tree : forestPart){
        //Predict probability of the tree 
        vector<vector<double>> treeProba = predictProba(X, {tree}); 

        //Lock access to shared memory 
        lock_guard<mutex> guard(lock); 

        //Add the tree's prediction into all+proba
        for (size_t i = 0; i < all_proba.size(); ++i) {
            for (size_t j = 0; j < all_proba[i].size(); ++j) {
                all_proba[i][j] += tree_proba[i][j];  // safely update
            }
        }
    }
}


/**predictProba (for 1D Array)
 * @param   X, a matrix of doubles (xIntTrans)
 * @return  matrice of doubles 
 */
vector<double> predictProba(const vector<vector<double>>& X, 
                                    const vector<vector<Node>>& forest)
{
    //1. CHECK TO SEE IF X IS VALID *Caution 
    // check_is_fitted  ////Throws an exception if data not valid
    //X = validateXPredict(X) ////Returns a better version of X, if not good 

    //2. ASSIGN JOBS (HOW MANY TREES EACH CPU CORE WILL HANDLE) *Caution 
    int selfNJobs = 1;
    pair<int, vector<int>> partEst = partitionEstimators(forest.size(), selfNJobs); 
    int nJobs = partEst.first; 
    //checking
    cout << "nJobs: " << nJobs << endl; 


    //3. ALLOCATE RESULT MATRIX w/zeroes  
    int XRows = X.size(); 
    int XColumns = X[0].size(); 
    int selfNClasses = 2; 
    //For a matrix of doubles 
    // vector<vector<double>> all_proba(XRows, vector<double>(selfNClasses, 0.0));
    //List of Matrices 
    vector<vector<vector<double>>> all_proba;
    all_proba.push_back(vector<vector<double>>(XRows, vector<double>(selfNClasses, 0.0)));


    //4. OPTIONAL PARALLEL STEP 
    if (nJobs==1){
        //Declare Lock 
        mutex lock; 
        //Run in parallel 
    //4. IF ONLY ONE JOB, just do the for loop 
    }else{
        //just loop through self.estimators and accumulate predictions 
    }

    //5. LOOP THROUGH EACH MATRICE and divide by amount of self.estimators 
    for(size_t i=0; i<all_proba.size(); i++){
        //proba /= self.estimators.size(); 
    }


    //6. RETURN 
    if (all_proba.size() == 1){return all_proba[0];}
    else{return all_proba;} 
}







//---------Other functions -----end -----







/**treeBasedPredict
 *
 * The C++ equivalent of 
 * 
 * @param       X, a matrix of doubles (XIntTrans in intBasedTrans)
 * @return      
 * 
 * */

vector<double> treeBasedPredict(const vector<vector<double>>& X){

    //1. SAVE SELF.EXTRA_TREES FROM PY CODE 
    vector<vector<Node>> forest = readTreesFromFile("/home/ccuev029/rSTSF_CPP/DATA/self.extra_trees.txt"); 
        //Checking 
        cout << "Number of trees: " << forest.size() << endl;
    
    //2. GET THE PROBABILITIES OF X 
    // vector<vector<double>> proba = predictProba(X, forest); 

    //Checking: probabilities 

    //3. CHECKING STEP  

    //4. RETURN PREDICTIONS  
    // return yPred; 

}