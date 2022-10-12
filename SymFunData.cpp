#include "SymFun.hpp"
#include <cstring>
#include <iostream>
#include <fstream>

#define MY_PI 3.1415926535897932

#define DIM 3

SymmetryFunctionParams::SymmetryFunctionParams() : has_three_body_(false) {
}

SymmetryFunctionParams::~SymmetryFunctionParams() {}

inline void SymmetryFunctionParams::set_species(std::vector<std::string> &species) {
    species_.resize(species.size());
    std::copy(species.begin(), species.end(), species_.begin());
};

inline void SymmetryFunctionParams::get_species(std::vector<std::string> &species) {
    species.resize(species_.size());
    std::copy(species_.begin(), species_.end(), species.begin());
};

inline int SymmetryFunctionParams::get_num_species() { return species_.size(); }

void SymmetryFunctionParams::set_cutoff(char const *name,
                                        std::size_t const Nspecies,
                                        double const *rcut_2D) {
    (void) name;   // to avoid unused warning
    rcut_2D_.resize(Nspecies, Nspecies, rcut_2D);
}

inline double SymmetryFunctionParams::get_cutoff(int const iCode, int const jCode) {
    return rcut_2D_(iCode, jCode);
};

void SymmetryFunctionParams::add_descriptor(char const *name,
                                            double const *values,
                                            int const row,
                                            int const col) {
    // Enzyme string comparison workaround
    if (strcmp(name, "g1") == 0) { name_.push_back(1); };
    if (strcmp(name, "g2") == 0) { name_.push_back(2); };
    if (strcmp(name, "g3") == 0) { name_.push_back(3); };
    if (strcmp(name, "g4") == 0) { name_.push_back(4); };
    if (strcmp(name, "g5") == 0) { name_.push_back(5); };

    Array2D<double> params(row, col, values);
    params_.push_back(std::move(params));

    auto sum = std::accumulate(num_param_sets_.begin(), num_param_sets_.end(), 0);
    starting_index_.push_back(sum);

    num_param_sets_.push_back(row);
    num_params_.push_back(col);

    if (strcmp(name, "g4") == 0 || strcmp(name, "g5") == 0) {
        has_three_body_ = true;
    }
}

int SymmetryFunctionParams::get_num_descriptors() {
    return std::accumulate(num_param_sets_.begin(), num_param_sets_.end(), 0);
}


void SymmetryFunctionParams::init(std::string file_name) {
    std::fstream file_ptr(file_name);
    std::string placeholder_string;
    int n_species;

    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    n_species = std::stoi(placeholder_string);

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');

    double * cutoff_matrix = new double[n_species * n_species];
    for (int i=0;i<n_species;i++){
        for (int j=0;j<n_species;j++){
            auto pos = placeholder_string.find(' ');
            *( cutoff_matrix + n_species * i + j) = std::stod(placeholder_string.substr(0, pos));
            if (pos != std::string::npos) placeholder_string.erase(0, pos + 1);
        }
        std::getline(file_ptr, placeholder_string);
    }

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    std::string cutoff_function = placeholder_string;

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    width = std::stoi(placeholder_string);
    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');
    int n_func = std::stoi(placeholder_string);

    std::vector<std::string> sym_func_list;
    for (int i = 0; i < n_func; i++) {
        std::getline(file_ptr, placeholder_string);
        sym_func_list.push_back(placeholder_string);
    }

    std::vector<int> sym_func_lengths;
    for (int i = 0; i < n_func; i++) {
        std::getline(file_ptr, placeholder_string);
        sym_func_lengths.push_back(std::stoi(placeholder_string));
    }

    std::vector<std::vector<double>> sym_func_elements;
    for (int i = 0; i < n_func; i++) {
        //blank line
        std::getline(file_ptr, placeholder_string);
        std::vector<double> tmp_desc_list;
        for (int j = 0; j < sym_func_lengths[i]; j++){
            std::getline(file_ptr, placeholder_string);
            tmp_desc_list.push_back(std::stod(placeholder_string));
        }
        sym_func_elements.push_back(std::move(tmp_desc_list));
    }

    for (int i = 0; i < n_func ; i++){
        if (sym_func_list[i] == "g2"){
            for (int j = 0; j < sym_func_elements[i].size(); j = j+2) sym_func_elements[i][j] /= (bhor2ang * bhor2ang);
        } else if (sym_func_list[i] == "g4") {
            for (int j = 2; j < sym_func_elements[i].size(); j = j+3) sym_func_elements[i][j] /= (bhor2ang * bhor2ang);
        }
    }

    // blank line
    std::getline(file_ptr, placeholder_string);
    // Ignore comments
    do {
        std::getline(file_ptr, placeholder_string);
    } while (placeholder_string[0] == '#');

    std::vector<int> dims;
    for (int i=0; i < n_func * 2; i++){
        dims.push_back(std::stoi(placeholder_string));
        std::getline(file_ptr, placeholder_string);
    }

    set_cutoff(cutoff_function.c_str(), n_species, cutoff_matrix);

    for (int i =0 ; i< n_func; i++){
        add_descriptor(sym_func_list[i].c_str(), sym_func_elements[i].data(), dims[2 *i], dims[2*i+1] );
    }
    delete[] cutoff_matrix;
}

#undef DIM
#undef MY_PI
