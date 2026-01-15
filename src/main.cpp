#include "../include/main.hpp"
#include "autoreg.hpp"
#include "intBasedT.hpp"
#include "treeBasedPredict.hpp"

//✦•··MAIN.CPP of rSTSF_CPP HELPER FUNCTIONS - START ···············•✦•·················•✦•·················•✦

/*FOR READING/WRITING TO FILES... */

/**writeMatrixToFile
 * @param   matrix to write
 *          name of file 
 * Matriix -> File 
**/
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

/**readMatrix (dtype double)
 * @param           string name of txt file
 * @return          matrix from txt file 
 * File -> Matrix(double)
**/ 
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

/**readVector (int vector)
 * @param           string name of txt file
 * @return          vector from txt file 
**/ 
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


/*FOR DEBUGGING*/

/** matrixMismatches 
 * @param matrixOne, matrixTwo, filename, precision 
 * @return errors  
 * catches the amount of mismatches between two matrices + writes mismatches to txt file 
 */
int compareMatrices(const vector<vector<double>>& matrixOne, 
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

//compares the contents of two inputted vectors, couts the 
//number of errors, doesn't write to file
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

//✦•··MAIN.CPP HELPER FUNCTIONS - END ·················•✦•···················•✦•···················•✦


/*MAIN FUNCTION PREDICT EQUIVALENT in rSTSF*/ 

int main() {

    //1. GETTING THE COPIED DATA NEEDED from r-STSF
    //For actual use
    vector<vector<double>> X_test = readMatrix("/home/ccuev029/DATA/XTest.txt");       
    vector<vector<double>> X_per = readMatrix("/home/ccuev029/DATA/XPer.txt");        
    vector<vector<double>> X_diff = readMatrix("/home/ccuev029/DATA/XDiff.txt");
    vector<int> relevantCaf = readVector("/home/ccuev029/DATA/relevant_CAF_idx.txt");  
    vector<vector<double>> allCaf = readMatrix("/home/ccuev029/DATA/allCAF.txt");
    //For comparison/debugging 
    vector<vector<double>> ar_X_test = readMatrix("/home/ccuev029/DATA/XAr.txt");      
    vector<vector<double>> X_Test_T = readMatrix("/home/ccuev029/DATA/XTest.txt");     
    vector<int> yTest = readVector("/home/ccuev029/DATA/YTest.txt"); 

    //1.1 𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Checking Sizes of Copied Data to make sure they match  
    cout << "\nChecking sizes: " << endl;
    cout << "Size of X_test: " << X_test.size() <<  "(rows) x " << X_test[0].size() << "(cols)" << endl; 
    cout << "Size of X_per: " << X_per.size() <<  "(rows) x " << X_per[0].size() << "(cols)" << endl; 
    cout << "Size of X_diff: " << X_diff.size() <<  "(rows) x " << X_diff[0].size() << "(cols)" << endl;
    cout << "Size of relevantCaf: " << relevantCaf.size() <<  endl;
    cout << "Size of allCaf: " << allCaf.size() << "(rows) x " << allCaf[0].size() << endl; 




    //2. COMPUTING XAr (CPP) (different from ar_X_test, computed in rSTSF_CPP)
    vector<vector<double>> X_ar = ar_coeffs(X_test); 
    //2.1Flip the Signs (for now)
    for (auto& row : X_ar) {                      
        for (auto& val : row) {
            val = -val;
        }
    }
    //2.2 Save generated XAr
    writeMatrixToFile(X_ar, "/home/ccuev029/rSTSF_CPP/DATA_CPP/XAr_cpp.txt");
    //2.3 𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: XAr (CPP) cs XAr (Py)
    cout << "This is XAr (Py)| Rows: " << ar_X_test.size() << " Cols: " << ar_X_test[0].size() << endl;
    cout << "This is XAr (CPP)| Rows: " << X_ar.size() << " Cols: " << X_ar[0].size() << endl;
    cout << "XAr (Py) vs XAr (CPP): " << compareMatrices(ar_X_test, X_ar, "/home/ccuev029/rSTSF_CPP/DEBUGGING/AR_Mismatch.txt", 16) << " / " << X_ar.size()*X_ar[0].size() << endl; 


    //3. GET INTERVAL BASED TRANSFORMATION  
    vector<vector<double>> XIntTrans = getIntervalBasedTransform(X_test, ar_X_test/*with Xar_PY*/, X_per, X_diff, allCaf, relevantCaf);
    //3.2 Save generated XIntTrans 
    writeMatrixToFile(XIntTrans, "/home/ccuev029/rSTSF_CPP/DATA_CPP/XIntTrans_cpp.txt");
    //3.3 𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: XIntTrans (PY) vs XIntTrans (CPP)
    cout << "This is XIntTrans (Py)| Rows: " << X_Test_T.size() << " Cols: " << X_Test_T[0].size() << endl;
    cout << "This is XIntTrans (CPP)| Rows: " << XIntTrans.size() << " Cols: " << XIntTrans[0].size() << endl;
    // cout << "XIntTrans (Py) vs XIntTrans (CPP): " << compareMatrices(X_Test_T, XIntTrans, "/home/ccuev029/rSTSF_CPP/DEBUGGING/xIntTrans_Mismatch.txt", 0) << endl;

    //IV. TREE BASED PREDICT (Y_PRED)
    vector<int> yPred = treeBasedPredict(XIntTrans); 
    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Printing yPred 
    cout << "yPred Size: " << yPred.size() << endl; 
    for (size_t i = 0; i < yPred.size(); ++i) {
        cout << yPred[i] << " ";
        if ((i + 1) % 30 == 0) {
            cout << endl; // Start a new line after every 30 elements
        }
    }
    // Optional: Final newline if total isn't a multiple of 30
    if (yPred.size() % 30 != 0) {
        cout << endl;
    }

    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: yPred vs yTest (for accuracy)
    int total = yTest.size(); 
    double errors = compareVectorsAndErrors(yPred, yTest); 
    double accuracy = (double)(total - errors) / total;
    cout << "accuracy: " << accuracy << endl;  //FINAL PRODUCT
}
 
