#include "../include/main.hpp"
#include "autoreg.hpp"
#include "intBasedT.hpp"

//--------------------------------------------//

/**writeMatrixToFile
 * @param   matrix to write
 *          name of file 
**/
void writeMatrixToFile(const vector<vector<double>>& matrix, const string& filename) {
    ofstream outFile(filename);

    if (!outFile.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }

    for (const auto& row : matrix) {
        for (const auto& val : row) {
            outFile << val <<  " ";
        }
        outFile << "\n"; // Newline after each row
    }

    outFile.close();
}

/**readMatrix (double matrices)
 * @param           string name of txt file
 * @return          matrix from txt file 
**/ 
vector<vector<double>> readMatrix(const string& filename){
    ifstream infile(filename);
    vector<vector<double>> matrix;
    string line;

    //Checking to see if file opened
    if (!infile.is_open()) 
    {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }


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

    //Check to see if matrix is empty 

    return matrix;
}

/**readVector (int vector)
 * @param           string name of txt file
 * @return          vector from txt file 
**/ 
vector<int> readVector(const string& filename){
    vector<int> vector;
    ifstream infile(filename);
    
    //Checking if the file is open
    if (!infile.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }

    int number;
    while (infile >> number) {
        vector.push_back(number);
    }

    infile.close();
    return vector;

}

/** matrixMismatches 
 * @param matrixOne, matrixTwo (the matrices to be compared)
 * @return the amount of errors  
 * catches the amount of mismatches between two matrices + writes mismatches to txt file 
 */
int matrixMismatches(const vector<vector<double>>& matrixOne, 
                        const vector<vector<double>>& matrixTwo, 
                        const string& filename, const int& precision){
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

//-----------Some Helper Functions -----------//





int main() {
    //I. GETTING THE COPIED DATA NEEDED from r-STSF
    vector<vector<double>> X_test = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/XtestData.txt");       //X_test: the original time series
    vector<vector<double>> ar_X_test = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/XarData.txt");      //ar_X_test: for comparison with X_ar
    vector<vector<double>> X_per = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/XperData.txt");        //Other Time Representations
    vector<vector<double>> X_diff = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/XdiffData.txt");
    vector<vector<double>> allCaf = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/all_candidate_agg_feats.txt");
    vector<int> relevantCaf = readVector("/home/ccuev029/rSTSF_CPP/DATA/relevant_caf_idx.txt");      //relevantCaf
    vector<vector<double>> X_Test_T = readMatrix("/home/ccuev029/rSTSF_CPP/DATA/X_test_T.txt");     //Transformed Matrix for comparison

    //II. COMPUTING AR REPRESENTATION.... (different from ar_X_test, computed in rSTSF_CPP)
    vector<vector<double>> X_ar = ar_coeffs(X_test); 
    for (auto& row : X_ar) {                      //Flip the signs
        for (auto& val : row) {
            val = -val;
        }
    }

    // //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Write X_ar to file for comparison 
    // writeMatrixToFile(X_ar, "/home/ccuev029/rSTSF_CPP/DATA/X_ar_cpp.txt");

    //Comparison:Original AR vs Computed AR 
    cout << "\nSize of X_ar: " << X_ar.size() <<  " " << X_ar[0].size() << endl; 
    cout << "Size of ar_X_Test: " << ar_X_test.size() << " " << ar_X_test[0].size() << endl; 
    cout << "ar_X_test vs X_ar Errors " << matrixMismatches(ar_X_test, X_ar, "/home/ccuev029/rSTSF_CPP/DEBUGGING/AR_Mismatch.txt", 16) << " / " << X_ar.size()*X_ar[0].size() << endl; 

    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Checking Sizes of Copied Data to make sure they match  
    cout << "\nSize of X_test: " << X_test.size() <<  " " << X_test[0].size() << endl; 
    cout << "Size of X_per: " << X_per.size() <<  " " << X_per[0].size() << endl; 
    cout << "Size of X_diff: " << X_diff.size() <<  " " << X_diff[0].size() << endl;

    //III. GET INTERVAL BASED TRANSFORMATION  
    cout << "\nWith X_ar: " << endl; 
    vector<vector<double>> XIntTrans = getIntervalBasedTransform(X_test, X_ar, X_per, X_diff, allCaf, relevantCaf);


    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Trying method with original X_ar Data to see if this works...
    // cout << "\nWith ar_X_test: " << endl; 
    // vector<vector<double>> XIntTrans = getIntervalBasedTransform(X_test, ar_X_test, X_per, X_diff, allCaf, relevantCaf);

    //Write transformed matrice to file 
    writeMatrixToFile(XIntTrans, "/home/ccuev029/rSTSF_CPP/DATA/XIntTransform.txt"); 
    
    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Checking sizes...
    cout << "\nSize of xIntTrans: " << XIntTrans.size() << " " << XIntTrans[0].size() << endl; 
    cout << "Size of X_Test_T: " << X_Test_T.size() << " " << X_Test_T[0].size() << endl;

    //Comparison! X_int_T in Python vs XIntTrans in C++ - writing to file
    cout << "X_Test_T vs xIntTrans errors: " << matrixMismatches(X_Test_T, XIntTrans, "/home/ccuev029/rSTSF_CPP/DEBUGGING/xIntTrans_Mismatch.txt", 8) << " / " << XIntTrans.size() * XIntTrans[0].size() << endl; 


    //IV. TREE BASED PREDICT (Y_PRED)
    yPred = treeBasedPredict(XIntTrans); 

    //Comparison: y_pred(py) vs yPred (c++)
}
 
