#include "intBasedT.hpp"


//np.arange C++ equialent *Caution
vector<int> arange(int start, int stop, int step = 1) {
    vector<int> result;
    for (int i = start; i < stop; i += step) {
        result.push_back(i);
    }
    return result;
}

//Inner Mean Equivalent **Corrected**
float innerMean(const vector<float>& X) {
    //number of columns = length of row of subinterval 
    int ncols = X.size(); 

    //accumulator to add elements  *Caution: float or int? 
    float accum = 0; 

    //add each element in the row together 
    for(size_t i=0; i<ncols; i++){
        accum+=X[i]; 
    }

    //mean = sum of row elements/num of elements 
    return accum/ncols; 
    
}

//Fast_Mean equivalent **Corrected**
vector<float> fastMean(const vector<vector<float>>& X){
    //Get number of rows and columns of the subinterval 
    int nrows = X.size(); 
    int _ = X[0].size(); 

    //Allocate the result vector with 0's 
    vector<float> _X(nrows, 0.0); 

    //For each row of the subinterval, calculate innermean of that row 
    for(size_t i=0; i<nrows; i++){
        _X[i] = innerMean(X[i]); 
    }

    //Return the result 
    return _X; 
}

//Inner Std Equivalent **Corrected**
float innerStd(const vector<float>& X) {
    //Get number of columns 
    int ncols = X.size(); 

    //Initialize accumulator 
    float accum = 0; 

    //Calculate the mean of that row 
    float X_mean = innerMean(X); 

    //std calculation  *Caution 
    for(size_t i; i<ncols; i++){
        accum += (X[i] - X_mean) * (X[i] - X_mean);
    }

    //return result 
    return sqrt(accum / ncols);

}

//Fast_std equivalent *Corrected*
vector<float> fastStd(const vector<vector<float>>& X){
    //Get number of rows and columns of subinterval 
    int nrows = X.size(); 
    int _ = X[0].size(); 

    //Allocate result vector, filled with zeroes
    vector<float> _X(nrows, 0.0); 

    //For each row in subinterval, calculate std of that row 
    for(size_t i=0; i<nrows; i++){
        _X[i] = innerStd(X[i]); 
    }

    //return result 
    return _X; 
}

//Inner Slope Equivalent *CORRECTED**
float innerSlope(const vector<float>& X, const vector<int>& indices) {
    //Get number of columns/elements in subinterval row 
    int ncols = X.size(); 

    //Initialize SUMS *Caution: Use float or integers? 
    float SUMx = 0; 
    float SUMy = 0; 
    float SUMxy = 0; 
    float SUMxx = 0; 

    //For each element in the subinterval row.... 
    for(size_t i=0; i<ncols; i++){
        SUMx = SUMx + indices[i];
        SUMy = SUMy + X[i];
        SUMxy = SUMxy + indices[i]*X[i];
        SUMxx = SUMxx + indices[i]*indices[i];
    }

    //Return the value 
    return ( SUMx*SUMy - ncols*SUMxy ) / ( SUMx*SUMx - ncols*SUMxx );
}

//Fast_Slope equivalent **CORRECTED** 
vector<float> fastSlope(const vector<vector<float>>& Y){

    //#Take # of rows and columns in the subinterval 
    int r = Y.size(); 
    int c = Y[0].size(); 

    //Create 1D array of integers(starts at 0, stops before c
    vector<int> x = arange(0, c); 

    //Create 1D array of 0's, of length r 
    vector<float> _X(r, 0.0); 

    //For each row in subinterval [], calculate slope of row elements 
    for (size_t i = 0; i < r; ++i) {
        _X[i] = innerSlope(Y[i], x); // Pass row i of Y and x
    }

    //Return Value 
    return _X;
}

//Inner IQR Equivalent 
float inner_iqr(vector<float> a) {
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

//Fast_Iqr equivalent 
vector<float> fast_iqr(const vector<vector<float>>& Y){
    //Declare the vector
    int rows = Y.size(); 
    vector<float> fast_iqr(rows, 0.0);
    //Loop through sub interval, calculate inner_iqr for each row, append to fast_iqr 
    for(size_t i=0; i<fast_iqr.size(); i++){
        fast_iqr[i] = inner_iqr(Y[i]); 
    }
    return fast_iqr; 
}

//Count_mean_crossing equivalent 
vector<float> count_mean_crossing(const vector<vector<float>>& X) {
    size_t nrows = X.size();
    vector<float> X_(nrows, 0.0);

    for (size_t i = 0; i < nrows; ++i) {
        const vector<float>& row = X[i];
        size_t n = row.size();

        // Compute mean
        float sum = accumulate(row.begin(), row.end(), 0.0);
        float mean = sum / n;

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

        X_[i] = static_cast<float>(crossings);
    }

    return X_;
}

//Count Values above mean equivalent 
vector<float> count_values_above_mean(const vector<vector<float>>& X) {
    size_t nrows = X.size();
    vector<float> X_(nrows, 0.0);

    for (size_t i = 0; i < nrows; ++i) {
        const vector<float>& row = X[i];
        size_t n = row.size();

        if (n == 0) {
            X_[i] = 0;
            continue;
        }

        // Compute mean
        float sum = accumulate(row.begin(), row.end(), 0.0);
        float mean = sum / n;

        // Count values above mean
        int count = 0;
        for (size_t j = 0; j < n; ++j) {
            if (row[j] > mean) {
                count++;
            }
        }

        X_[i] = static_cast<float>(count);
    }

    return X_;
}

//ELSE COMPUTATIONS


//innerMedian, computing the median of the row  *Caution
float innerMedian(vector<float> Y) {
    //Get the number of elements in row of subinterval; 
    int n = Y.size(); 

    //Sort the vector 
    sort(Y.begin(), Y.end()); 

    //If odd number of elements in row, return middle element
    if(n%2 == 1){
        float median = Y[n / 2];
        return median; 
    //If even number of elements in row, average of two middle elements 
    } else {
        float median = ((Y[n/2-1] + Y[n/2])/2.0); 
        return median;
    }
    
}

//Median, getting a vector of medians *Caution
vector<float> median(const vector<vector<float>>& X){
    //Getting #rows and columns of subinterval 
    int rows = X.size(); 
    int columns = X[0].size(); 

    //Declare the result matrix 
    vector<float> _X(rows, 0.0); 

    //For each row within the sub interval, calculate median and append to result matrix 
    for(size_t i=0; i<rows; i++){
        _X[i] = innerMedian(X[i]); 
    }

    //return result matrix 
    return _X; 
}

//innerMin : Getting the minimum value in a given row 
float innerMin(const vector<float>& Y){
    //Compute the min value of the given row directly 
    float min = Y[0]; 
    for(size_t i=0; i<Y.size(); i++){
        if(Y[i]<min){
            min = Y[i];
        }
    }
    //return min value 
    return min; 
}

//Min : Getting a vector of minimum values  
vector<float> min(const vector<vector<float>>& X) {
    //Get #rows and columns of subinterval 
    int rows = X.size(); 
    int columns = X[0].size(); 

    //Declare result array 
    vector<float> _X(rows, 0.0);

    //For each row of subinterval, compute min and append to result array 
    for(size_t i=0; i<rows; i++){
        _X[i] = innerMin(X[i]); 
    }

    return _X; 
}

//innerMax : Max value in a given row 
float innerMax(const vector<float>& Y){
    //Compute the max value of the given row directly 
    float max = Y[0]; 
    for(size_t i=0; i<Y.size(); i++){
        if(Y[i]>max){
            max = Y[i];
        }
    }
    //return min value 
    return max; 
}

//Max : Getting a vector of max values 
vector<float> max(const vector<vector<float>>& X) {
    //Get rows and columns of subinterval matrix 
    int rows = X.size(); 
    int columns = X[0].size(); 

    //Declare result array (the return)
    vector<float> _X(rows, 0.0); 

    //For each row of subInt, calculate max value, append to array 
    for(size_t i=0; i<rows; i++){
        _X[i] = innerMax(X[i]); 
    }

    return _X; 
}

//**CHECK ALL OF THE ABOVE LATER */



/**getIntervalFeature
 * @param       sub_interval, agg_fn 
 * @return      vector of floats with corresponding calculation
 * **/
vector<float> getIntervalFeature(vector<vector<float>> sub_interval, int agg_fn){
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
 * @param    allCaf  //??
 *             #an empty list
 * @return   X_test_T: _Array[tuple[int, int], floating[_32Bit]] 
 *          in C++, this is a 2D dynamic array? 
 * 
 * 
 * **/
vector<vector<float>> getIntervalBasedTransform(vector<vector<float>> X,
                                                 vector<vector<float>> X_ar,
                                                 vector<vector<float>> X_per,
                                                 vector<vector<float>> X_diff,
                                                 vector<vector<vector<float>>> allCaf
                                                 )
{


    //Allocate the Result Matrix with zeroes (rows of X, columns of allCaf)
    size_t numRows = X.size();
    size_t numColumns = allCaf[0].size();
    vector<vector<float>> XIntTrans(numRows, vector<float>(numColumns, 0.0));

    //For each jth item in relavantCaf (actual value) SF (Save features)
    for(/*int j : relevantCaf*/ int j=0; j<allCaf.size(); j++){//*Caution with j 

        //Save each element of each row into corresponding variables 
        // float w        = allCaf[0][j][0];
        // float score    = allCaf[0][j][1];
        int li = static_cast<int>(allCaf[0][j][0]);     //*Caution
        int ls = static_cast<int>(allCaf[0][j][1]);     //*Caution
        int agg_fn      = static_cast<int>(allCaf[0][j][2]);
        int repr_type   = static_cast<int>(allCaf[0][j][3]);

        vector<vector<float>> X_temp;
        if (repr_type == 0)      { X_temp = X; }
        else if (repr_type == 1) { X_temp = X_per; }
        else if (repr_type == 2) { X_temp = X_ar; }
        else if (repr_type == 3) { X_temp = X_diff; }
        else { throw invalid_argument("Invalid repr_type in getIntervalBasedTransform at " + to_string(j) + " value: " + to_string(repr_type)); }

        //Getting the subintervals 
        vector<vector<float>> sub_interval(numRows);
        //For each subinterval 
        for (size_t row = 0; row < numRows; ++row) {
            sub_interval[row] = vector<float>(X_temp[row].begin() + li, X_temp[row].begin() + ls);
        }


        //getIntervalFeature on subinterval
        vector<float> to_add = getIntervalFeature(sub_interval, agg_fn);
        //put into return matrix
        for (size_t row = 0; row < numRows; ++row) {
            XIntTrans[row][j] = to_add[row];
        }
        
        // //Debugging: Keeping Track of subInterval Selections...
        // cout << "Processing column (j): " << j << std::endl;
        // cout << "Size of sub_interval: " << sub_interval.size() << " " << sub_interval[0].size() << endl; 
        // cout << "li: " << li << " ls: " << ls <<endl;
        // cout << "Chosen representations: " << repr_type << endl; 
        // cout <<"Chosen agg_fn: " << agg_fn << "\n" << endl; 

        // //Setting precision 
        // for (auto& row : XIntTrans) {
        //     for (auto& value : row) {
        //         // value = static_cast<int64_t>(value * 1e8) / 1e8; //Truncate without rounding
        //         value = round(value * 1e8) / 1e8; //*Caution 
        //     }
        // }
    }
    
    //Return the result matrix (the transformed X_test)
    return XIntTrans;
}
