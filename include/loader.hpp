#ifndef MINIROCKET_HLS_TESTBENCH_LOADER_H
#define MINIROCKET_HLS_TESTBENCH_LOADER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "treeBasedPredict.hpp"


// Testbench-only loader (NOT synthesizable - for simulation only)
class Loader {
private:
    // Simple JSON parser for arrays
    std::vector<double> parse_double_array(const std::string& content, const std::string& key);
    std::vector<int> parse_int_array(const std::string& content, const std::string& key);
    std::vector<std::vector<int>> parse_2d_int_array(const std::string& content, const std::string& key);
    std::vector<std::vector<double>> parse_2d_double_array(const std::string& content, const std::string& key);
    std::vector<std::vector<std::vector<double>>> parse_3d_double_array(const std::string& content, const std::string& key);
    std::vector<std::vector<Node>> parse_trees(const std::string& content);
    int parse_int_value(const std::string& content, const std::string& key);
    
    std::string read_file(const std::string& filename);
    void trim_whitespace(std::string& str);
    
public:
    // Load test data for verification (testbench only)
    bool load_test_data(
        const std::string& test_filename, 
        std::vector<std::vector<double>>& X_test,
        std::vector<std::vector<double>>& X_Diff,
        std::vector<std::vector<double>>& X_Per,
        std::vector<std::vector<double>>& X_Ar,
        std::vector<std::vector<std::vector<double>>>& all_caf,
        std::vector<std::vector<Node>>& trees,
        std::vector<int>& y_test,
        std::vector<int>& y_pred
    );
};

#endif // MINIROCKET_HLS_TESTBENCH_LOADER_H