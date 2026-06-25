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

// int main() {

//     //I. GETTING THE COPIED DATA NEEDED from r-STSF
//     //For actual use
//     vector<vector<double>> X_test = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/XtestData.txt");       //X_test: the original time series
//     vector<vector<double>> X_per = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/XperData.txt");        //Other Time Representations
//     vector<vector<double>> X_diff = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/XdiffData.txt");
//     vector<int> relevantCaf = readVector("/home/ccuev029/rSTSF_CPP/DATA/relevant_caf_idx.txt");      //relevantCaf
//     vector<vector<double>> allCaf = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/all_candidate_agg_feats.txt");
    
//     //For comparison/debugging 
//     vector<vector<double>> ar_X_test = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/XarData.txt");      //ar_X_test: for comparison with X_ar
//     vector<vector<double>> X_Test_T = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/X_test_T.txt");     //Transformed Matrix for comparison
//     vector<int> yTest = readVector("/home/ccuev029/rSTSF_CPP/DATA/y_test.txt"); 

//     //II. COMPUTING AR REPRESENTATION.... (different from ar_X_test, computed in rSTSF_CPP)
//     vector<vector<double>> X_ar = ar_coeffs(X_test); 
//     for (auto& row : X_ar) {                      //Flip the signs
//         for (auto& val : row) {
//             val = -val;
//         }
//     }

//     // //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Write X_ar to file for comparison 
//     // writeMatrixToFile(X_ar, "/home/ccuev029/rSTSF_CPP/DATA/X_ar_cpp.txt");

//     //Comparison:Original AR vs Computed AR 
//     cout << "\nSize of X_ar: " << X_ar.size() <<  " " << X_ar[0].size() << endl; 
//     cout << "Size of ar_X_Test: " << ar_X_test.size() << " " << ar_X_test[0].size() << endl; 
//     cout << "ar_X_test vs X_ar Errors " << writeMatrixMismatches(ar_X_test, X_ar, "/home/ccuev029/rSTSF_CPP/DEBUGGING/AR_Mismatch.txt", 16) << " / " << X_ar.size()*X_ar[0].size() << endl; 

//     //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Checking Sizes of Copied Data to make sure they match  
//     cout << "\nSize of X_test: " << X_test.size() <<  " " << X_test[0].size() << endl; 
//     cout << "Size of X_per: " << X_per.size() <<  " " << X_per[0].size() << endl; 
//     cout << "Size of X_diff: " << X_diff.size() <<  " " << X_diff[0].size() << endl;

//     //III. GET INTERVAL BASED TRANSFORMATION  
//     cout << "\nWith X_ar: " << endl; 
//     vector<vector<double>> XIntTrans = getIntervalBasedTransform(X_test, X_ar, X_per, X_diff, allCaf, relevantCaf);


//     //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Trying method with original X_ar Data to see if this works...
//     // cout << "\nWith ar_X_test: " << endl; 
//     // vector<vector<double>> XIntTrans = getIntervalBasedTransform(X_test, ar_X_test, X_per, X_diff, allCaf, relevantCaf);

//     //Write transformed matrice to file 
//     writeMatrixToFile(XIntTrans, "/home/ccuev029/rSTSF_CPP/DATA/XIntTransform.txt"); 
    
//     // //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Checking sizes...
//     // cout << "\nSize of xIntTrans: " << XIntTrans.size() << " " << XIntTrans[0].size() << endl; 
//     // cout << "Size of X_Test_T: " << X_Test_T.size() << " " << X_Test_T[0].size() << endl;

//     // //Comparison! X_int_T in Python vs XIntTrans in C++ - writing to file
//     // cout << "X_Test_T vs xIntTrans errors: " << matrixMismatches(X_Test_T, XIntTrans, "/home/ccuev029/rSTSF_CPP/DEBUGGING/xIntTrans_Mismatch.txt", 8) << " / " << XIntTrans.size() * XIntTrans[0].size() << endl; 


//     //IV. TREE BASED PREDICT (Y_PRED)
//     vector<int> yPred = treeBasedPredict(XIntTrans); 

//     // // //Debugging: Using Original Transformed Data
//     // vector<int> yPred = treeBasedPredict(X_Test_T); 

//     //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Printing yPred 
//     cout << "yPred Size: " << yPred.size() << endl; 
//     for (size_t i = 0; i < yPred.size(); ++i) {
//         std::cout << yPred[i] << " ";
        
//         if ((i + 1) % 30 == 0) {
//             std::cout << std::endl; // Start a new line after every 30 elements
//         }
//     }
//     // Optional: Final newline if total isn't a multiple of 30
//     if (yPred.size() % 30 != 0) {
//         cout << endl;
//     }

//     //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: yPred vs yTest (for accuracy)
//     int total = yTest.size(); 
//     double errors = compareVectorsAndErrors(yPred, yTest); 
//     double accuracy = (double)(total - errors) / total;
//     cout << "accuracy: " << accuracy << endl; 
// }
 


int main(int argc, char* argv[]) {
    // Example usage of Loader

    // 4 representations of the time series
    vector<vector<float>> X_Test;
    vector<vector<float>> X_Diff;
    vector<vector<float>> X_Per;
    vector<vector<float>> X_Ar;

    //Transform features
    vector<vector<vector<float>>> all_caf;
    vector<int> relevant_caf_idx;

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
        cout << "Number of expected outputs: " << y_test.size() << endl;
        cout << "Number of relevant CAF indices: " << relevant_caf_idx.size() << endl;


    //1. Getting X_Ar Representation 

    //2. Interval Based Transform 
    vector<vector<float>> XIntTrans = getIntervalBasedTransform(X_Test, X_Ar, X_Per, X_Diff, all_caf);//testing with no relCaf
    cout << "getIntervalBasedTransform executed successfully." << endl;

    //3. Tree Based Predict
    vector<int> yPred = treeBasedPredict(XIntTrans, trees); //Error
    cout << "treeBasedPredict executed successfully." << endl;

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