#ifndef __UCSKY_H__
#define __UCSKY_H__

#include "Tuple.h"
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
    bool sorted;
    int nR;
    // uc-sky
    int l;
    vector<int> ucsky;
    BigFloat ucsky_prob;
    // dominant graph
    int nG;
    vector<vector<int>> inN;
    vector<vector<int>> outN;
    vector<BigFloat> beta;


public:
    UCSKY_Solver(const int& dim, const int& n, const double& center, const int& l);
    virtual ~UCSKY_Solver();
    // data operation
    void gen_ind_data();
    void gen_ind_data_2d(const int& m);
    void gen_anti_data();
    // void gen_anti_data_2d(const int& m);
    void gen_corr_data();
    // void gen_corr_data_2d(const int& m);
    void load_data(const char* data_path);
    void write_data(const char* data_path);
    void print_data();

    // algorithm
    void loop_bsl();
    void basic_dp();
    void basic_dp2(); // without data reduction
    void basic_dp3(); // vertically dp (gamma space from n^2 to n, further from n to 1)
    void group_dp();
    void group_dp2(); // without data reduction, for testing partition
    void group_dp3(); // group_dp2 plus pruning, vertically dp
    bool greedy();
    bool greedy_plus(const double& step); // memory too high
    int greedy_plus2(const double& step); // tighter lb after sky_2n
    bool sky_2n();
    void branch_bound();
    void sample_FIM(const int& theta);
    void sample_FIM2(const int& theta); // direct sample + skyline
    int reduce_data(const double& step);
    int reduce_data_bsl();

    // helper
    void print_ucsky();
    void check_ucsky();
    void enum_l_subset(int j, vector<int>& subset, int i, BigFloat& prob);
    bool conflict(const vector<int>& S, const int& t);
    BigFloat compute_prob(const vector<int>& S);
    BigFloat compute_prob(const set<int>& S);
    void divide(vector<pair<int, double>>& T);
    void enum_subset(const vector<pair<int, double>>& T, const vector<pair<int, double>>& sky, int s, int j, vector<int>& subset, int i, BigFloat& prob);
    void skyline(vector<pair<int, double>>& T, vector<pair<int, double>>& sky);
    void search(const vector<int>& T, const vector<int>& S, const BigFloat& prob, const int& k);
    bool contain(const vector<int>& T, const vector<int>& subset);
};

#endif