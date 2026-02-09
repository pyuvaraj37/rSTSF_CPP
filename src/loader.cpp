#include "loader.hpp"

std::string Loader::read_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return "";
    }
    
    std::string content, line;
    while (std::getline(file, line)) {
        content += line;
    }
    return content;
}

void Loader::trim_whitespace(std::string& str) {
    str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
    str.erase(str.find_last_not_of(" \t\n\r\f\v") + 1);
}

int Loader::parse_int_value(const std::string& content, const std::string& key) {
    size_t key_pos = content.find("\"" + key + "\"");
    if (key_pos == std::string::npos) return -1;
    
    size_t colon_pos = content.find(":", key_pos);
    size_t value_start = colon_pos + 1;
    size_t value_end = content.find_first_of(",}", value_start);
    
    std::string value_str = content.substr(value_start, value_end - value_start);
    trim_whitespace(value_str);
    
    return std::stoi(value_str);
}

std::vector<int> Loader::parse_int_array(const std::string& content, const std::string& key) {
    std::vector<int> result;
    
    size_t key_pos = content.find("\"" + key + "\"");
    if (key_pos == std::string::npos) return result;
    
    size_t array_start = content.find("[", key_pos);
    size_t array_end = content.find("]", array_start);
    
    if (array_start == std::string::npos || array_end == std::string::npos) return result;
    
    std::string array_content = content.substr(array_start + 1, array_end - array_start - 1);
    
    std::stringstream ss(array_content);
    std::string item;
    
    while (std::getline(ss, item, ',')) {
        trim_whitespace(item);
        if (!item.empty()) {
            result.push_back(std::stoi(item));
        }
    }
    
    return result;
}

std::vector<float> Loader::parse_float_array(const std::string& content, const std::string& key) {
    std::vector<float> result;
    
    size_t key_pos = content.find("\"" + key + "\"");
    if (key_pos == std::string::npos) return result;
    
    size_t array_start = content.find("[", key_pos);
    size_t array_end = content.find("]", array_start);
    
    if (array_start == std::string::npos || array_end == std::string::npos) return result;
    
    std::string array_content = content.substr(array_start + 1, array_end - array_start - 1);
    
    std::stringstream ss(array_content);
    std::string item;
    
    while (std::getline(ss, item, ',')) {
        trim_whitespace(item);
        if (!item.empty()) {
            result.push_back(std::stof(item));
        }
    }
    
    return result;
}

std::vector<std::vector<int>> Loader::parse_2d_int_array(const std::string& content, const std::string& key) {
    std::vector<std::vector<int>> result;
    
    size_t key_pos = content.find("\"" + key + "\"");
    if (key_pos == std::string::npos) return result;
    
    size_t array_start = content.find("[", key_pos);
    size_t array_end = array_start;
    int bracket_count = 0;
    
    // Find matching closing bracket
    for (size_t i = array_start; i < content.length(); i++) {
        if (content[i] == '[') bracket_count++;
        if (content[i] == ']') bracket_count--;
        if (bracket_count == 0) {
            array_end = i;
            break;
        }
    }
    
    std::string array_content = content.substr(array_start + 1, array_end - array_start - 1);
    
    // Parse each sub-array
    size_t pos = 0;
    while (pos < array_content.length()) {
        size_t sub_start = array_content.find("[", pos);
        if (sub_start == std::string::npos) break;
        
        size_t sub_end = array_content.find("]", sub_start);
        if (sub_end == std::string::npos) break;
        
        std::string sub_array = array_content.substr(sub_start + 1, sub_end - sub_start - 1);
        
        std::vector<int> sub_result;
        std::stringstream ss(sub_array);
        std::string item;
        
        while (std::getline(ss, item, ',')) {
            trim_whitespace(item);
            if (!item.empty()) {
                sub_result.push_back(std::stoi(item));
            }
        }
        
        result.push_back(sub_result);
        pos = sub_end + 1;
    }
    
    return result;
}

std::vector<std::vector<float>> Loader::parse_2d_float_array(const std::string& content, const std::string& key) {
    std::vector<std::vector<float>> result;
    
    size_t key_pos = content.find("\"" + key + "\"");
    if (key_pos == std::string::npos) return result;
    
    size_t array_start = content.find("[", key_pos);
    size_t array_end = array_start;
    int bracket_count = 0;
    
    // Find matching closing bracket
    for (size_t i = array_start; i < content.length(); i++) {
        if (content[i] == '[') bracket_count++;
        if (content[i] == ']') bracket_count--;
        if (bracket_count == 0) {
            array_end = i;
            break;
        }
    }
    
    std::string array_content = content.substr(array_start + 1, array_end - array_start - 1);
    
    // Parse each sub-array
    size_t pos = 0;
    while (pos < array_content.length()) {
        size_t sub_start = array_content.find("[", pos);
        if (sub_start == std::string::npos) break;
        
        size_t sub_end = array_content.find("]", sub_start);
        if (sub_end == std::string::npos) break;
        
        std::string sub_array = array_content.substr(sub_start + 1, sub_end - sub_start - 1);
        
        std::vector<float> sub_result;
        std::stringstream ss(sub_array);
        std::string item;
        
        while (std::getline(ss, item, ',')) {
            trim_whitespace(item);
            if (!item.empty()) {
                sub_result.push_back(std::stof(item));
            }
        }
        
        result.push_back(sub_result);
        pos = sub_end + 1;
    }
    
    return result;
}


std::vector<std::vector<std::vector<float>>>
Loader::parse_3d_float_array(const std::string& content,
                             const std::string& key)
{
    std::vector<std::vector<std::vector<float>>> result;

    size_t key_pos = content.find("\"" + key + "\"");
    if (key_pos == std::string::npos) return result;

    size_t array_start = content.find("[", key_pos);
    if (array_start == std::string::npos) return result;

    size_t array_end = array_start;
    int bracket_count = 0;

    // Find matching closing bracket for the whole 3D array
    for (size_t i = array_start; i < content.size(); i++) {
        if (content[i] == '[') bracket_count++;
        else if (content[i] == ']') bracket_count--;

        if (bracket_count == 0) {
            array_end = i;
            break;
        }
    }

    std::string array_content =
        content.substr(array_start + 1, array_end - array_start - 1);

    // ---- parse 2D slices ----
    size_t pos = 0;
    while (pos < array_content.size()) {

        // find next [ ... ] (this is a 2D slice)
        size_t slice_start = array_content.find("[", pos);
        if (slice_start == std::string::npos) break;

        size_t slice_end = slice_start;
        bracket_count = 0;

        for (size_t i = slice_start; i < array_content.size(); i++) {
            if (array_content[i] == '[') bracket_count++;
            else if (array_content[i] == ']') bracket_count--;

            if (bracket_count == 0) {
                slice_end = i;
                break;
            }
        }

        std::string slice =
            array_content.substr(slice_start + 1, slice_end - slice_start - 1);

        // ---- parse rows inside slice ----
        std::vector<std::vector<float>> slice_result;

        size_t row_pos = 0;
        while (row_pos < slice.size()) {

            size_t row_start = slice.find("[", row_pos);
            if (row_start == std::string::npos) break;

            size_t row_end = slice.find("]", row_start);
            if (row_end == std::string::npos) break;

            std::string row =
                slice.substr(row_start + 1, row_end - row_start - 1);

            std::vector<float> row_vals;
            std::stringstream ss(row);
            std::string item;

            while (std::getline(ss, item, ',')) {
                trim_whitespace(item);
                if (!item.empty()) {
                    row_vals.push_back(std::stof(item));
                }
            }

            slice_result.push_back(row_vals);
            row_pos = row_end + 1;
        }

        result.push_back(slice_result);
        pos = slice_end + 1;
    }

    return result;
}

std::vector<std::vector<Node>>
Loader::parse_trees(const std::string& content)
{
    std::vector<std::vector<Node>> result;

    size_t key_pos = content.find("\"trees\"");
    if (key_pos == std::string::npos) return result;

    size_t array_start = content.find("[", key_pos);
    if (array_start == std::string::npos) return result;

    // ---- find matching closing bracket for the entire trees array ----
    int bracket_count = 0;
    size_t array_end = array_start;

    for (size_t i = array_start; i < content.size(); i++) {
        if (content[i] == '[') bracket_count++;
        else if (content[i] == ']') bracket_count--;

        if (bracket_count == 0) {
            array_end = i;
            break;
        }
    }

    std::string trees_content =
        content.substr(array_start + 1, array_end - array_start - 1);

    // ---- parse each tree ----
    size_t pos = 0;
    while (pos < trees_content.size()) {

        size_t tree_start = trees_content.find("[", pos);
        if (tree_start == std::string::npos) break;

        // find end of this tree
        bracket_count = 0;
        size_t tree_end = tree_start;

        for (size_t i = tree_start; i < trees_content.size(); i++) {
            if (trees_content[i] == '[') bracket_count++;
            else if (trees_content[i] == ']') bracket_count--;

            if (bracket_count == 0) {
                tree_end = i;
                break;
            }
        }

        std::string tree_content =
            trees_content.substr(tree_start + 1, tree_end - tree_start - 1);

        std::vector<Node> tree_nodes;

        // ---- parse nodes inside tree ----
        size_t npos = 0;
        while (npos < tree_content.size()) {

            size_t node_start = tree_content.find("[", npos);
            if (node_start == std::string::npos) break;

            size_t node_end = tree_content.find("]", node_start);
            if (node_end == std::string::npos) break;

            std::string node_str =
                tree_content.substr(node_start + 1, node_end - node_start - 1);

            // split first 4 scalars + last array
            size_t values_start = node_str.find("[");
            std::string header = node_str.substr(0, values_start);
            std::string values =
                node_str.substr(values_start + 1,
                                node_str.find("]", values_start) - values_start - 1);

            // ---- parse header ----
            std::stringstream hs(header);
            std::string token;
            Node node;

            std::getline(hs, token, ',');
            node.feature = std::stoi(token);

            std::getline(hs, token, ',');
            node.threshold = std::stod(token);

            std::getline(hs, token, ',');
            node.left = std::stoi(token);

            std::getline(hs, token, ',');
            node.right = std::stoi(token);

            // ---- parse values ----
            std::stringstream vs(values);
            node.values.clear();

            while (std::getline(vs, token, ',')) {
                trim_whitespace(token);
                if (!token.empty())
                    node.values.push_back(std::stod(token));
            }

            tree_nodes.push_back(node);
            npos = node_end + 1;
        }

        result.push_back(tree_nodes);
        pos = tree_end + 1;
    }

    return result;
}



bool Loader::load_test_data(
    const std::string& test_filename, 
    std::vector<std::vector<float>>& X_test,
    std::vector<std::vector<float>>& X_Diff,
    std::vector<std::vector<float>>& X_Per,
    std::vector<std::vector<float>>& X_Ar,
    std::vector<std::vector<std::vector<float>>>& all_caf,
    std::vector<std::vector<Node>>& trees,
    std::vector<int>& y_test,
    std::vector<int>& y_pred
) {
    std::string content = read_file(test_filename);
    if (content.empty()) return false;
    
    X_test = parse_2d_float_array(content, "X_test");
    X_Diff = parse_2d_float_array(content, "X_Diff");
    X_Per = parse_2d_float_array(content, "X_Per");
    X_Ar = parse_2d_float_array(content, "X_Ar");
    all_caf = parse_3d_float_array(content, "all_candidate_agg_feats");
    trees = parse_trees(content);
    y_test = parse_int_array(content, "y_test");
    y_pred = parse_int_array(content, "y_pred");

    std::cout << "Test data loaded: " << X_test.size() << " samples" << std::endl;
    return true;
}