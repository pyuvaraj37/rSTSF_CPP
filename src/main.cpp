#include "../include/main.hpp"
#include "autoreg.hpp"
#include "intBasedT.hpp"
#include "treeBasedPredict.hpp"
#include "loader.hpp"


//✦•··MAIN.CPP HELPER FUNCTIONS - START ····················•✦•······················•✦•······················•✦

void writeMatrixToFile(const vector<vector<double>>& matrix, const string& filename) {
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

vector<vector<double>> readMatrix(const string& filename){
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

vector<int> readVector(const string& filename){
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


//catches amount of mismatches between two matrices, writes to file, returns number of errors
int writeMatrixMismatches(const vector<vector<double>>& matrixOne, 
                        const vector<vector<double>>& matrixTwo, 
                        const string& filename, const int& precision){
    cout << "calling writeMatrixMismatches..." << endl; 
    //Open a file to write to 
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }
    //Checking for same size 
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

// compares the contents of two inputted vectors, couts the number of errors, doesn't write to file
int compareVectorsAndErrors(const vector<int>& A, const vector<int>& B){
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



/*MAIN FUNCTION PREDICT EQUIVALENT in rSTSF*/ 
 


int main(int argc, char* argv[]) {
    // Example usage of Loader

    // 4 representations of the time series
    vector<vector<double>> X_Test;
    vector<vector<double>> X_Diff;
    vector<vector<double>> X_Per;
    vector<vector<double>> X_Ar;
    //Transform features
    vector<vector<vector<double>>> all_caf;
    vector<int> relevant_caf_idx; //Not used anymore
    //Tree 
    vector<vector<Node>> trees;
    // Expected outputs
    vector<int> y_test; // For compatibility with expected_outputs
    vector<int> y_pred; 

    //Load the test data from json 
    std::string test_filename = argv[1];
    std::cout << "Test filename: " << test_filename << std::endl;

    Loader loader;
    if (loader.load_test_data(test_filename, X_Test, X_Diff, X_Per, X_Ar, all_caf, trees, y_test, y_pred)) {
        cout << "Test inputs and expected outputs loaded successfully." << endl;
        cout << "Number of test samples: " << X_Test.size() << endl;
        cout << "X_Diff size: " << X_Diff.size() << " x " << X_Diff[0].size() << endl;
        cout << "X_Per size: " << X_Per.size() << " x " << X_Per[0].size() << endl;
        cout << "X_Ar size: " << X_Ar.size() << " x " << X_Ar[0].size() << endl;
        cout << "all_caf size: " << all_caf.size() << " x " << all_caf[0].size() << endl;
        cout << "Number of trees: " << trees.size() << endl;
        cout << "y_test size: " << y_test.size() << endl;


        //2. Interval Based Transform
        vector<vector<double>> XIntTrans = getIntervalBasedTransform(X_Test, X_Ar, X_Per, X_Diff, all_caf);
        cout << "getIntervalBasedTransform executed successfully." << endl;


        //3. Tree Based Predict
        vector<int> yPred = treeBasedPredict(XIntTrans, trees); //Error
        cout << "treeBasedPredict executed successfully." << endl;


        //Debug Print
        cout << "Debug Information:" << endl; 
        cout << "X_Test shape: " << X_Test.size() << " x " << X_Test[0].size() << endl;
        cout << "X_Diff shape" << X_Diff.size() << " x " << X_Diff[0].size() << endl;
        cout << "X_Per shape: " << X_Per.size() << " x " << X_Per[0].size() << endl;
        cout << "X_Ar shape: " << X_Ar.size() << " x " << X_Ar[0].size() << endl;
        cout << "all_caf shape: " << all_caf.size() << " x " << all_caf[0].size() << endl;
        cout << "XIntTrans shape: " << XIntTrans.size() << " x " << XIntTrans[0].size() << endl;
        cout << "YTest size: " << y_test.size() << endl;
        cout << "YPred size: " << yPred.size() << endl;
        cout << "Number of trees: " << trees.size() << endl;


        //4. Compare yPred with y_test
        int total = y_test.size();
        double errors = compareVectorsAndErrors(yPred, y_test);
        double accuracy = (double)(total - errors) / total;
        cout << "Accuracy: " << accuracy << endl;


    
    } else {
        cout << "Failed to load test data." << endl;
    }
    
    return 0;
}