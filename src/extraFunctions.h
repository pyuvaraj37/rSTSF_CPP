//For all my helper functions
#include "../include/main.hpp"
#include "autoreg.hpp"
#include "intBasedT.hpp"
#include "treeBasedPredict.hpp"


//✦•··MAIN.CPP HELPER FUNCTIONS - START ····················•✦•······················•✦•······················•✦

/*FOR READING/WRITING TO FILES... */

/**writeMatrixToFile
 * @param   matrix to write
 *          name of file 
 * Matriix -> File 
**/
inline void writeMatrixToFile(const vector<vector<double>>& matrix, const string& filename) {
    cout << "calling writeMatrixToFile..." << endl; 
    ofstream outFile(filename);
    //Make sure file is open 
    if (!outFile.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }
    //Writing to file 
    for (const auto& row : matrix) {
        for (const auto& val : row) {
            outFile << val <<  " ";
        }
        outFile << "\n"; // Newline after each row
    }
    //Close File 
    outFile.close();
}

/**readMatrix (dtype double)
 * @param           string name of txt file
 * @return          matrix from txt file 
 * File -> Matrix(double)
**/ 
inline vector<vector<double>> readMatrix(const string& filename){
    cout << "calling readMatrix.." << endl; 
    ifstream infile(filename);
    vector<vector<double>> matrix;
    string line;
    //Checking to see if file is Open 
    if (!infile.is_open()) 
    {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    //Reading each Line, 
    while (getline(infile, line)) {
        istringstream iss(line);
        vector<double> row;
        double val;

        while (iss >> val) {
            row.push_back(val);
        }

        // Optional: check column count consistency
        if (!matrix.empty() && row.size() != matrix[0].size()) {
            cerr << "Error: Inconsistent number of columns!" << endl;
            exit(1);
        }

        matrix.push_back(row);
    }
    return matrix;
}

/**readVector (int vector)
 * @param           string name of txt file
 * @return          vector from txt file 
**/ 
inline vector<int> readVector(const string& filename){
    cout << "calling readVector..." << endl; 
    vector<int> vector;
    ifstream infile(filename);
    //Checking if the file is open
    if (!infile.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }
    //Reading the file..
    int number;
    while (infile >> number) {
        vector.push_back(number);
    }
    infile.close();
    return vector;
}


/*FOR DEBUGGING*/

/** matrixMismatches 
 * catches the amount of mismatches between two matrices + writes mismatches to txt file 
 */
inline int writeMatrixMismatches(const vector<vector<double>>& matrixOne, 
                        const vector<vector<double>>& matrixTwo, 
                        const string& filename, const int& precision){

    //Open a file to write to 
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }
    //Checking for same size 
    cout << matrixOne.size() << " vs. " << matrixTwo.size() << endl;
    if (matrixOne.size() != matrixTwo.size() || matrixOne[0].size() != matrixTwo[0].size()) {
        throw runtime_error("Matrix size mismatch: Cannot compare.");
    }

    //Actual Comparison (with epsilon threshold)
    double epsilon = 1e-8;
    outFile << "With Epsilon Tolerance: " << epsilon << endl; 
    int errors = 0; 
    for(size_t i = 0; i<matrixOne.size(); i++){
        for(size_t j=0; j<matrixOne[0].size(); j++){
            double a = matrixOne[i][j];
            double b = matrixTwo[i][j];
            if (fabs(a - b) > epsilon) {
                errors++;
                outFile << "Mismatch at (" << i << ", " << j << ") : "
                        << fixed << setprecision(precision)
                        << a << " vs. " << b << endl;
            }

        }
    }
    outFile << "Errors: " << errors << " / " << matrixOne.size()*matrixOne[0].size() << endl; 
    return errors; 

    
}






/** compareVectors <int>
 * compares the contents of two inputted vectors, couts the 
 * number of errors, doesn't write to file**/
inline int compareVectorsAndErrors(const vector<int>& A, const vector<int>& B){
    cout << "calling compare vectors" << endl; 
    int errors = 0; 
    //Size comparison 
    int aSize = A.size(); 
    int bSize = B.size(); 
    if(aSize != bSize){
        cout << "Vectors are not the same size! Size of A: " << aSize
            << " Size of B: " << bSize << endl; 
    }
    //Check Individual Elements 
    for(size_t i=0; i<A.size(); i++) {
        if(A[i]!=B[i]){
            errors++; 
        }
    }
    //Return errors 
    cout << "Errors: " << errors << endl; 
    return errors; 
}

//✦•··MAIN.CPP HELPER FUNCTIONS - END ····················•✦•······················•✦•······················•✦





//---------HELPER FUNCTIONS----START--------from treebasedPredict

/**readTreesFromFile
 * @param       filePath   ,file to be read from 
 * @return      allTrees   ,matrix of Nodes 
 * 
 * **/
inline vector<vector<Node>> readTreesFromFile(const string& filePath) {
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
inline void checkIsFitted(const vector<vector<Node>>& forest) {
    if (forest.empty()) {
        throw runtime_error("Forest is empty. Make sure trees are loaded.");
    }
}


/**validateXPredict (complete ai)
 * @param   X, original matrix 
 * @return  vector<vector<double>> X, validated matrix 
 * **/
inline vector<vector<double>> validateXPredict(const vector<vector<double>>& X, bool allow_nan = false) {
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
inline vector<vector<vector<double>>> getTreeProba(const vector<vector<Node>>& forest, 
                                            const vector<vector<double>>& X, 
                                            vector<vector<vector<double>>>& all_proba, 
                                            const int& numSamples, 
                                            const int& numClasses, 
                                            const int& numOut){
        for (size_t i=0; i<numOut; i++){ //Output i 
            for (size_t j = 0; j < forest.size(); j++) {  // Tree j
                const auto& tree = forest[j];
                for (size_t k=0; k<numSamples; k++){   //Sample k 
                    //Get row of X 
                    vector<double> XRow = X[k];

                    //TREE TRAVERSAL STEP 
                    int nodeIdx = 0;           //Start at root node 
                    while (tree[nodeIdx].left != -1 && tree[nodeIdx].right != -1){
                        const Node& node = tree[nodeIdx]; 
                        double xValue = XRow[node.feature];
                        if(xValue <= node.threshold){
                            nodeIdx = node.left; 
                        }else{
                            nodeIdx = node.right; 
                        }
                    }

                    //Once a leaf node is hit, get the values 
                    const vector<double>& probs = tree[nodeIdx].values;

                    //Summing the probabilities and adding to all_proba 
                    for (size_t c = 0; c < probs.size(); ++c) {
                        all_proba[i][k][c] += probs[c];
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
inline vector<vector<vector<double>>> predictProba(const vector<vector<double>>& X_input, 
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

inline void printForest(vector<vector<Node>> forest){
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
}




//---------HELPER FUNCTIONS  ----- END -----