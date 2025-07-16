#include "../include/main.hpp"
#include "autoreg.hpp"
#include "intBasedT.hpp"

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
        throw std::runtime_error("Could not open file: " + filename);
    }

    int number;
    while (infile >> number) {
        vector.push_back(number);
    }

    infile.close();
    return vector;

}



int main() {
    //GETTING THE COPIED DATA NEEDED from r-STSF
    vector<vector<double>> X_test = readMatrix("/home/ccuev029/rSTSF_CPP/COPIED DATA/XtestData.txt");    //X_test: the original time series
    vector<vector<double>> ar_X_test = readMatrix("/home/ccuev029/rSTSF_CPP/COPIED DATA/XarData.txt");   //ar_X_test: for comparison with X_ar
    vector<vector<double>> X_per = readMatrix("/home/ccuev029/rSTSF_CPP/COPIED DATA/XperData.txt");      //Other Time Representations
    vector<vector<double>> X_diff = readMatrix("/home/ccuev029/rSTSF_CPP/COPIED DATA/XdiffData.txt");
    vector<vector<double>> allCaf = readMatrix("/home/ccuev029/rSTSF_CPP/COPIED DATA/all_candidate_agg_feats.txt");
    vector<int> relevantCaf = readVector("/home/ccuev029/rSTSF_CPP/COPIED DATA/relevant_caf_idx.txt");   //relevantCaf
    vector<vector<double>> X_Test_T = readMatrix("/home/ccuev029/rSTSF_CPP/COPIED DATA/X_test_T.txt");  //Transformed Matrix for comparison


    //GETTING THE AR REPRESENTATION.... (different from ar_X_test, computed in rSTSF_CPP)
    vector<vector<double>> X_ar = ar_coeffs(X_test); 
    for (auto& row : X_ar) {                      //Flip the signs
        for (auto& val : row) {
            val = -val;
        }
    }

    //Debugging! Checking Sizes of Copied Data 
    cout << "Size of X_test: " << X_test.size() <<  " " << X_test[0].size() << endl; 
    cout << "Size of X_ar: " << X_ar.size() <<  " " << X_ar[0].size() << endl; 
    cout << "Size of X_per: " << X_per.size() <<  " " << X_per[0].size() << endl; 
    cout << "Size of X_diff: " << X_diff.size() <<  " " << X_diff[0].size() << endl;

    //GET INTERVAL BASED TRANSFORMATION  
    vector<vector<double>> XIntTrans = getIntervalBasedTransform(X_test, X_ar, X_per, X_diff, allCaf, relevantCaf);

    //Debugging: 
    cout << "Size of xIntTrans: " << XIntTrans.size() << " " << XIntTrans[0].size() << endl; 
    cout << "Size of X_Test_T: " << X_Test_T.size() << " " << X_Test_T[0].size() << endl;

    //Check! Printing XIntTrans
    for (int i = 0; i < XIntTrans.size(); i++) {
        for (int j = 0; j < XIntTrans[0].size(); j++) {
            cout << XIntTrans[i][j] << " ";
        }
        cout << endl; // Move to the next line after printing a row
    }


    //Check! X_int_T in Python vs XIntTrans in C++
    int errors = 0; 
    for (size_t i = 0; i < X_Test_T.size(); ++i) {
        for (size_t j = 0; j < X_Test_T[i].size(); ++j) {
            if (X_Test_T[i][j] != XIntTrans[i][j]){
                errors++; 
            }
        }
    }
    cout << "Errors: " << errors << endl; 

}
 
