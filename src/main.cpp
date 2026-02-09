#include "extraFunctions.h"

/*MAIN FUNCTION PREDICT EQUIVALENT in rSTSF*/ 

int main() {

    //I. GETTING THE COPIED DATA NEEDED from r-STSF
    //For actual use
    vector<vector<double>> X_test = readMatrix("/home/ccuev029/DATA/XTest.txt");            //X_test: the original time series
    vector<vector<double>> X_per = readMatrix("/home/ccuev029/DATA/XPer.txt");              //Other Time Representations
    vector<vector<double>> X_diff = readMatrix("/home/ccuev029/DATA/XDiff.txt");
    vector<int> relevantCaf = readVector("/home/ccuev029/DATA/relevant_CAF_idx.txt");       //relevantCaf
    vector<vector<double>> allCaf = readMatrix("/home/ccuev029/DATA/allCAF.txt");
    //For comparison/debugging 
    vector<vector<double>> ar_X_test = readMatrix("/home/ccuev029/DATA/XAr.txt");           //ar_X_test: for comparison with X_ar
    vector<vector<double>> X_Test_T = readMatrix("/home/ccuev029/DATA/XIntTrans.txt");     //Transformed Matrix for comparison
    vector<int> yPred_og = readVector("/home/ccuev029/DATA/YPred.txt");                        //YPred
    vector<int> yTest = readVector("/home/ccuev029/DATA/YTest.txt");                        //YTest

    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: Checking Sizes of Copied Data to make sure they match  
    cout << "\nSize of X_test: " << X_test.size() <<  " " << X_test[0].size() << endl; 
    cout << "Size of X_per: " << X_per.size() <<  " " << X_per[0].size() << endl; 
    cout << "Size of X_diff: " << X_diff.size() <<  " " << X_diff[0].size() << endl;



    //II. COMPUTING AR REPRESENTATION.... (different from ar_X_test, computed in rSTSF_CPP)
    vector<vector<double>> X_ar = ar_coeffs(X_test); 
    for (auto& row : X_ar) {                      //Flip the signs
        for (auto& val : row) {
            val = -val;
        }
    }

    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: 
    //Write X_ar to file for visual comparison
    writeMatrixToFile(X_ar, "/home/ccuev029/rSTSF_CPP/DATA/X_ar_cpp.txt"); 
    //Comparison to original
    cout << "ar_X_test vs X_ar errors: " << writeMatrixMismatches(ar_X_test, X_ar, "/home/ccuev029/rSTSF_CPP/DEBUGGING/AR_Mismatch.txt", 16) << " / " << X_ar.size()*X_ar[0].size() << endl; 


    //III. GET INTERVAL BASED TRANSFORMATION  
    cout << "\nWith X_ar: " << endl; 
    vector<vector<double>> XIntTrans = getIntervalBasedTransform(X_test, ar_X_test, X_per, X_diff, allCaf, relevantCaf);

    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: 
    //Write to file for visual comparison 
    writeMatrixToFile(XIntTrans, "/home/ccuev029/rSTSF_CPP/DATA_CPP/XIntTrans_cpp.txt"); 
    //Comparison to original
    cout << "X_Test_T vs xIntTrans errors: " << writeMatrixMismatches(X_Test_T, XIntTrans, "/home/ccuev029/rSTSF_CPP/DEBUGGING/xIntTrans_Mismatch.txt", 8) << " / " << XIntTrans.size() * XIntTrans[0].size() << endl; 




    
    //IV. TREE BASED PREDICT (Y_PRED)
    vector<int> yPred = treeBasedPredict(XIntTrans); 

    //𝓓𝓮𝓫𝓾𝓰𝓰𝓲𝓷𝓰 𓆣⊹ ࣪ 𖢥: 
    //Printing yPred 
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
    //Comparison to original
    compareVectorsAndErrors(yPred, yPred_og);

    //FINAL ACCURACY 
    int total = yTest.size(); 
    double errors = compareVectorsAndErrors(yPred, yTest); 
    double accuracy = (double)(total - errors) / total;
    cout << "accuracy: " << accuracy << endl; 
}
 