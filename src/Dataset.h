#ifndef __DATASET_H__
#define __DATASET_H__

#include "Tuple.h"

// dominant graph
class DG {
public:
    int n;
    vector<vector<int>> layers;
    unordered_map<int, unordered_set<int>> inN;
    unordered_map<int, unordered_set<int>> outN;

    DG(const int& n);
    virtual ~DG();
    
};

class Dataset {
public:
    int dim;
    int n;
    double center;
    string data_path;
    vector<Tuple> tuples;

    Dataset(const int& dim, const int& n, const double& center);
    virtual ~Dataset();

    void gen_ind_data();
    void gen_anti_data();
    void gen_corr_data();
    void load_data(const char* data_path);
    void write_data(const char* data_path);
    void print_data();
    DG construct_graph();
};

#endif