#include "intBasedT.hpp"


//Inner IQR Equivalent 
double inner_iqr(vector<double> a) {
    size_t n = a.size();
    if (n == 0) return 0.0;

    sort(a.begin(), a.end());

    if (n % 2 != 0) {
        size_t median_idx = n / 2;
        size_t q1_idx = median_idx / 2;
        size_t q3_idx = ((n - median_idx) / 2) + median_idx;
        return a[q3_idx] - a[q1_idx];
    } else {
        size_t median_idx_lower = (n - 1) / 2;
        size_t median_idx_upper = n / 2;
        size_t q1_idx = median_idx_lower / 2;
        size_t q3_idx = (median_idx_upper / 2) + median_idx_upper;
        return a[q3_idx] - a[q1_idx];
    }
}

//Inner Slope Equivalent (ai)
double innerSlope(const vector<double>& Y, const vector<int>& x) {
    size_t n = x.size();
    if (n != Y.size() || n == 0) return 0.0;

    double mean_x = 0.0, mean_y = 0.0;
    for (size_t i = 0; i < n; ++i) {
        mean_x += x[i];
        mean_y += Y[i];
    }
    mean_x /= n;
    mean_y /= n;

    double numerator = 0.0;
    double denominator = 0.0;
    for (size_t i = 0; i < n; ++i) {
        numerator   += (x[i] - mean_x) * (Y[i] - mean_y);
        denominator += (x[i] - mean_x) * (x[i] - mean_x);
    }

    if (denominator == 0.0) return 0.0;

    double slope = numerator / denominator;
    return slope;
}

//Fast_Slope equivalent **CHECK THIS LATER** 
vector<double> fastSlope(const vector<vector<double>>& Y){

    //#Take # of rows and columns in the subinterval 
    int rows = Y.size(); 
    int columns = Y[0].size(); 

    //np.arange equivalent (create array of indices)
    vector<int> x(columns);
    for (int i = 0; i < columns; ++i) {
        x[i] = i;
    }

    //np.zeroes eqivalent: create array of 0's (this is return value); 
    vector<double> _X(rows, 0.0);
    //Calculate inner slope for every row against indices  
    for(size_t i=0; i<Y.size(); i++){
        _X[i] = innerSlope(Y[i],x); 
    }

    //return appended matrix 
    return _X; 

}

//Inner Mean Equivalent 
double innerMean(const vector<double>& Y) {
    size_t columns = Y.size();
    double accum = 0.0;
    for (size_t i = 0; i < columns; ++i) {
        accum += Y[i];
    }
    double sum = accum/columns; 
    return sum;
}

//Fast_Mean equivalent **Check this later** + Change Variable later
vector<double> fastMean(const vector<vector<double>>& Y){
    //Get number of rows and columns of subInt
    int rows = Y.size(); 
    int columns = Y[0].size(); 
    //Declare result array of 0's (Will return this)
    vector<double> _X(rows, 0.0); 
    //For each row of subInterval, calculate inner mean and append! 
    for(size_t i=0; i<rows; i++){
        _X[i] = innerMean(Y[i]); 
    }
    //return fastSlope vector 
    return _X;
}

//Inner Std Equivalent 
double innerStd(const vector<double>& X) {
    //Get number of columns 
    int columns = X.size(); 

    //Set accum for sums 
    double accum = 0; 

    //Get the mean of X (the row of the subinterval)
    double X_mean = innerMean(X); 

    // Compute sum of squared deviations
    for(int i = 0; i < columns; i++) {
        accum += (X[i] - X_mean) * (X[i] - X_mean);
    }

    // Compute standard deviation
    double std = sqrt(accum / columns);
    return std; 
}

//Fast_std equivalent 
vector<double> fastStd(const vector<vector<double>>& Y){
    //Get number of rows and columns within the sub interval 
    int rows = Y.size(); 
    int columns = Y[0].size(); 

    //Allocate an array of zeroes with rows same # as sub interval (return value)
    vector<double> _X(rows, 0.0); 

    //For each row in the sub interval...
    for(size_t i=0; i<rows; i++){
        //Calculate the innerStd and append to result matrix
        _X[i] = innerStd(Y[i]); 
    }
    
    //return the vector of std! 
    return _X; 
}

//Fast_Iqr equivalent 
vector<double> fast_iqr(const vector<vector<double>>& Y){
    //Declare the vector
    int rows = Y.size(); 
    vector<double> fast_iqr(rows, 0.0);
    //Loop through sub interval, calculate inner_iqr for each row, append to fast_iqr 
    for(size_t i=0; i<fast_iqr.size(); i++){
        fast_iqr[i] = inner_iqr(Y[i]); 
    }
    return fast_iqr; 
}

//Count_mean_crossing equivalent 
vector<double> count_mean_crossing(const vector<vector<double>>& X) {
    size_t nrows = X.size();
    vector<double> X_(nrows, 0.0);

    for (size_t i = 0; i < nrows; ++i) {
        const vector<double>& row = X[i];
        size_t n = row.size();

        // Compute mean
        double sum = accumulate(row.begin(), row.end(), 0.0);
        double mean = sum / n;

        // Create boolean sequence: row[j] > mean ? 1 : 0
        vector<int> above_mean(n);
        for (size_t j = 0; j < n; ++j) {
            above_mean[j] = row[j] > mean ? 1 : 0;
        }

        // Count number of sign changes (crossings)
        int crossings = 0;
        for (size_t j = 1; j < n; ++j) {
            if (above_mean[j] != above_mean[j - 1]) {
                crossings++;
            }
        }

        X_[i] = static_cast<double>(crossings);
    }

    return X_;
}

//Count Values above mean equivalent 
vector<double> count_values_above_mean(const vector<vector<double>>& X) {
    size_t nrows = X.size();
    vector<double> X_(nrows, 0.0);

    for (size_t i = 0; i < nrows; ++i) {
        const vector<double>& row = X[i];
        size_t n = row.size();

        if (n == 0) {
            X_[i] = 0;
            continue;
        }

        // Compute mean
        double sum = accumulate(row.begin(), row.end(), 0.0);
        double mean = sum / n;

        // Count values above mean
        int count = 0;
        for (size_t j = 0; j < n; ++j) {
            if (row[j] > mean) {
                count++;
            }
        }

        X_[i] = static_cast<double>(count);
    }

    return X_;
}

//innerMedian, computing the median of the row  
double innerMedian(vector<double> Y) {
    //Get the number of elements in row of subinterval; 
    int n = Y.size(); 

    //Sort the vector 
    sort(Y.begin(), Y.end()); 

    //If odd number of elements in row, return middle element
    if(n%2 == 1){
        double median = Y[n / 2];
        return median; 
    //If even number of elements in row, average of two middle elements 
    } else {
        double median = ((Y[n/2-1] + Y[n/2])/2.0); 
        return median;
    }
    
}

//Median, getting a vector of medians
vector<double> median(const vector<vector<double>>& X){
    //Getting #rows and columns of subinterval 
    int rows = X.size(); 
    int columns = X[0].size(); 

    //Declare the result matrix 
    vector<double> _X(rows, 0.0); 

    //For each row within the sub interval, calculate median and append to result matrix 
    for(size_t i=0; i<rows; i++){
        _X[i] = innerMedian(X[i]); 
    }

    //return result matrix 
    return _X; 
}

//innerMin : Getting the minimum value in a given row 
double innerMin(const vector<double>& Y){
    //Compute the min value of the given row directly 
    double min = Y[0]; 
    for(size_t i=0; i<Y.size(); i++){
        if(Y[i]<min){
            min = Y[i];
        }
    }
    //return min value 
    return min; 
}

//Min : Getting a vector of minimum values  
vector<double> min(const vector<vector<double>>& X) {
    //Get #rows and columns of subinterval 
    int rows = X.size(); 
    int columns = X[0].size(); 

    //Declare result array 
    vector<double> _X(rows, 0.0);

    //For each row of subinterval, compute min and append to result array 
    for(size_t i=0; i<rows; i++){
        _X[i] = innerMin(X[i]); 
    }

    return _X; 
}

//innerMax : Max value in a given row 
double innerMax(const vector<double>& Y){
    //Compute the max value of the given row directly 
    double max = Y[0]; 
    for(size_t i=0; i<Y.size(); i++){
        if(Y[i]>max){
            max = Y[i];
        }
    }
    //return min value 
    return max; 
}

//Max : Getting a vector of max values 
vector<double> max(const vector<vector<double>>& X) {
    //Get rows and columns of subinterval matrix 
    int rows = X.size(); 
    int columns = X[0].size(); 

    //Declare result array (the return)
    vector<double> _X(rows, 0.0); 

    //For each row of subInt, calculate max value, append to array 
    for(size_t i=0; i<rows; i++){
        _X[i] = innerMax(X[i]); 
    }

    return _X; 
}

//**CHECK ALL OF THE ABOVE LATER */



/**getIntervalFeature
 * @param       sub_interval, agg_fn 
 * @return      vector of doubles with corresponding calculation
 * **/
vector<double> getIntervalFeature(vector<vector<double>> sub_interval, int agg_fn){
    if(agg_fn == 0){ 
        //0. polyfit -> fastSlope
        return fastSlope(sub_interval);

    }else if(agg_fn == 1){ 
        //1. mean -> fastMean
        return fastMean(sub_interval);

    }else if(agg_fn == 2){ 
        //2. std -> fastStd
        return fastStd(sub_interval);

    }else if(agg_fn == 3){ 
        //3. median -> calculate directly
        return median(sub_interval); 

    }else if(agg_fn == 4){ 
        //4. min -> calculate directly 
        return min(sub_interval);

    }else if(agg_fn == 5){ 
        //5. max -> calculate directly 
        return max(sub_interval);

    }else if(agg_fn == 6){ 
        //6. iqr -> fastIqr
        return fast_iqr(sub_interval); 

    }else if(agg_fn == 7){ 
        //7. percentile -> countMeanCrossing
        return count_mean_crossing(sub_interval);

    }else if(agg_fn == 8){ 
        //8. quantile -> countValuesAbove Mean 
        return count_values_above_mean(sub_interval); 

    }else{
        throw invalid_argument("Invalid agg_fn in getIntervalFeature");
    }
}   


/** getIntervalBasedTransform
 * @param    X, X_ar, X_per, X_diff
 * @param    allCaf, relevantCaf  //??
 *             #an empty list
 * @return   X_test_T: _Array[tuple[int, int], floating[_32Bit]] 
 *          in C++, this is a 2D dynamic array? 
 * 
 * 
 * **/
vector<vector<double>> getIntervalBasedTransform(vector<vector<double>> X,
                                                 vector<vector<double>> X_ar,
                                                 vector<vector<double>> X_per,
                                                 vector<vector<double>> X_diff,
                                                 vector<vector<double>> allCaf,
                                                 vector<int> relevantCaf)
{
    //Declare matrix with rows same number as original time series 
    //and columns same number as allCaf, this is the transformed matrix 
    size_t numRows = X.size();
    size_t numColumns = allCaf.size();
    vector<vector<double>> XIntTrans(numRows, vector<double>(numColumns, 0.0));

    //For each row in relevantCaf....
    for (size_t i = 0; i < relevantCaf.size(); ++i) {

        //Save each element of each row into corresponding variables 
        int j = relevantCaf[i];
        double w        = allCaf[j][0];
        double score    = allCaf[j][1];
        int li = static_cast<int>(allCaf[j][2]);     //*Caution
        int ls = static_cast<int>(allCaf[j][3]);     //*Caution
        int agg_fn      = static_cast<int>(allCaf[j][4]);
        int repr_type   = static_cast<int>(allCaf[j][5]);

        //Debugging 
        cout << "li: " << li << " ls: " << ls << endl; 


        vector<vector<double>> X_temp;
        if (repr_type == 1)      { X_temp = X; }
        else if (repr_type == 2) { X_temp = X_per; }
        else if (repr_type == 3) { X_temp = X_ar; }
        else if (repr_type == 4) { X_temp = X_diff; }
        else { throw invalid_argument("Invalid repr_type in getIntervalBasedTransform"); } 


        // //Debugging: for subinterval indices
        // for (size_t row = 0; row < numRows; ++row) {
        //     if (li >= X_temp[row].size() || ls > X_temp[row].size() || li >= ls) {
        //         throw out_of_range("Invalid sub-interval indices");
        //     }
        // }

        //Getting the subintervals 
        vector<vector<double>> sub_interval(numRows);
        //For each subinterval 
        for (size_t row = 0; row < numRows; ++row) {
            sub_interval[row] = vector<double>(X_temp[row].begin() + li, X_temp[row].begin() + ls);
        }

        //Debugging: Size of sub interval 
        cout << "Size of sub_interval: " << sub_interval.size() << " " << sub_interval[0].size() << endl; 


        vector<double> to_add = getIntervalFeature(sub_interval, agg_fn);

        for (size_t row = 0; row < numRows; ++row) {
            XIntTrans[row][j] = to_add[row];
        }


    }

    return XIntTrans;
}
