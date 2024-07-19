#ifndef __UCSKY_H__
#define __UCSKY_H__

#include "Dataset.h"
#include "unordered_map"
#include "Heap.h"
#include "fptree.h"

class UCSKY_Solver {
private:
    // dataset
    int dim;
    int n;
    double center;
    string data_path;
    vector<Tuple> tuples;
    // uc-sky
    int l;
    vector<int> ucsky;
    BigFloat ucsky_prob;
    // dominant graph
    // int nG;
    // vector<int> L0;
    // unordered_map<int, unordered_set<int>> inN;
    // unordered_map<int, unordered_set<int>> outN;

    bool dominate(const vector<int>& set, const int& t);
    bool conflict(const vector<int>& set, const int& t);
    void enum_l_subset(const vector<pair<int, double>>& T, const vector<int>& sky, int s, int j, vector<int>& subset, int i);
    void enum_l_subset(int j, vector<int>& subset, int i, BigFloat& prob);
    void enum_l_subset(const vector<int>& T, int j, vector<int>& subset, int i, map<vector<int>, int>& freq);
public:
    UCSKY_Solver(const int& dim, const int& n, const double& center, const int& l);
    virtual ~UCSKY_Solver();
    // data operation
    void gen_ind_data();
    void gen_anti_data();
    void gen_corr_data();
    void load_data(const char* data_path);
    void write_data(const char* data_path);
    void print_data();
    
    void construct_dg();

    // algorithm
    void loop_bsl();
    void basic_dp();
    void group_dp();

    BigFloat greedy();
    BigFloat alpha_greedy(const double& step);
    void branch_bound(const double& step);
    void divide(vector<pair<int, double>>& T);
    void bb(const vector<int>& T, const vector<int>& S, const BigFloat& prob, const int& k);
    void sample_fim(const int& theta);

    void print_ucsky();
    void check_prob();
    BigFloat compute_prob(const set<int>& s);
    BigFloat compute_prob(const vector<int>& s);
    void skyline(vector<pair<int, double>>& T, vector<int>& sky); // sky is removed from T
};

#endif