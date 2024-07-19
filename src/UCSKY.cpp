#include "UCSKY.h"

UCSKY_Solver::UCSKY_Solver(const int& dim, const int& n, const double& center, const int& l) {
    srand((unsigned)time(nullptr));
    srand48((unsigned)time(nullptr));
    this->dim = dim;
    this->n = n;
    this->center = center;
    this->l = l;
    this->tuples.resize(n, Tuple(-1, dim, nullptr, 0));
    ucsky_prob = BigFloat(0);
    ucsky.resize(l);
    stringstream stream;
    stream << fixed << setprecision(2) << center;
    data_path = to_string(dim) + "_" + stream.str();
}

UCSKY_Solver::~UCSKY_Solver() {

}

void UCSKY_Solver::gen_ind_data() {
    double delta = min(center, 1 - center);
    double lower = center - delta, upper = center + delta;
    if (dim == 2) { // non-repeating

    } else {
        for (int i = 0; i < n; ++ i) {
            tuples[i].id = i;
            for (int j = 0; j < dim; ++ j) {
                tuples[i].coord[j] = rand_uniform(0, 1);
            }
            tuples[i].prob = rand_uniform(lower, upper);
        }
    }
}

void UCSKY_Solver::gen_anti_data() {
    double delta = min(center, 1 - center);
    double lower = center - delta, upper = center + delta;
    double x[dim];
    if (dim == 2) {

    } else {
        for (int i = 0; i < n; ++ i) {
            tuples[i].id = i;
            double range = 0.5*dim + rand_normal(0, 0.05);
            tuples[i].coord[0] = rand_uniform(0, 1)*min(1.0, range);
            for (int j = 1; j < dim; ++ j) {
                range -= tuples[i].coord[j - 1];
                tuples[i].coord[j] = rand_uniform(0, 1)*min(1.0, range);
                if (j == dim - 1) tuples[i].coord[j] = min(1.0, range);
            }
            tuples[i].prob = rand_uniform(lower, upper);
        }
    }
}

void UCSKY_Solver::gen_corr_data() {
    double delta = min(center, 1 - center);
    double lower = center - delta, upper = center + delta;
    double x[dim];
    if (dim == 2) {

    } else {
        for (int i = 0; i < n; ++ i) {
            tuples[i].id = i;
            tuples[i].coord[0] = rand_uniform(0, 1);
            for (int j = 1; j < dim; ++ j) {
                do {
                    tuples[i].coord[j] = rand_normal(tuples[i].coord[0], 0.05);
                } while (tuples[i].coord[j] < 0 || tuples[i].coord[j] > 1);
            }
            tuples[i].prob = rand_uniform(lower, upper);
        }
    }
}

void UCSKY_Solver::load_data(const char* data_path) {
    if (strcmp(data_path, "data/iip.dat") == 0) {
        ifstream file("data/iip.dat", ios::in);
        if (!file.is_open()) {
            printf("Fail in opening data file.\n");
            exit(1);
        }
        for (int i = 0; i < n; ++ i) {
            tuples[i].id = i;
            file >> tuples[i].coord[0];
            file >> tuples[i].coord[1];
            file >> tuples[i].prob;
        }
        file.close();
    } else if (strcmp(data_path, "data/car.dat") == 0) {

    } else if (strcmp(data_path, "data/nba.dat") == 0) {

    } else {
        ifstream file((string(data_path) + this->data_path + ".dat").c_str(), ios::in);
        if (!file.is_open()) {
            printf("Fail in opening data file.\n");
            exit(1);
        }
        int nMax = 0;
        file >> nMax;
        if (n > nMax) {
            printf("n exceeds nMax.\n");
            exit(1);
        }
        for (int i = 0; i < n; ++ i) {
            tuples[i].id = i;
            for (int j = 0; j < dim; ++ j) {
                file >> tuples[i].coord[j];
            }
            file >> tuples[i].prob;
        }
        file.close();
    }
}

void UCSKY_Solver::write_data(const char* data_path) {
    ofstream file((string(data_path) + this->data_path + ".dat").c_str(), ios::out);
    file << n << " ";
    for (int i = 0; i < n; ++ i) {
        for (int j = 0; j < dim; ++ j) {
            file << tuples[i].coord[j] << " ";
        }
        file << tuples[i].prob << " ";
    }
    file.close();
}

void UCSKY_Solver::print_data() {
    for (int i = 0; i < tuples.size(); ++ i) {
        cout << i << "," << tuples[i] << endl;
    }
}

void UCSKY_Solver::construct_dg() {
    int nG;
    vector<int> L0;
    unordered_map<int, unordered_set<int>> inN;
    unordered_map<int, unordered_set<int>> outN;

    nG = n;
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        return t.sum() < s.sum();
    });
    for (int i = 0; i < nG; ++ i) {
        if (!dominate(L0, i)) {
            L0.push_back(i);
        }
        inN[i] = unordered_set<int>{};
        outN[i] = unordered_set<int>{};
        for (int j = 0; j < i; ++ j) {
            if (tuples[j].dominate(tuples[i])) {
                outN[j].insert(i);
                inN[i].insert(j);
            }
        }
    }
    // for (int i = 0; i < nG; ++ i) {
    //     cout << i << " degree = " << nG - inN[i].size() - outN[i].size() - 1 << endl;
    // }
}

bool UCSKY_Solver::dominate(const vector<int>& set, const int& t) {
    for (auto s : set) {
        if (tuples[s].dominate(tuples[t])) return true;
    }
    return false;
}

bool UCSKY_Solver::conflict(const vector<int>& set, const int& t) {
    for (auto s : set) {
        if (s != -1 && 
        (tuples[s].dominate(tuples[t]) || tuples[t].dominate(tuples[s]))) {
            return true;
        }
    }
    return false;
}

void UCSKY_Solver::enum_l_subset(int j, vector<int>& subset, int i, BigFloat& prob) {
    if (i == l) {
        BigFloat prob_real = prob;
        for (int p = 0; p < n; ++ p) {
            for (int q : subset) {
                if (tuples[p].dominate(tuples[q])) {
                    // cout << tuples[p] << endl;
                    prob_real = prob_real * (1 - tuples[p].prob);
                    break;
                }
            }
            if (prob_real < ucsky_prob) return;
        }
        if (prob_real > ucsky_prob) {
            // cout << "prob = " << prob_real << " ucsky_prob = " << ucsky_prob << endl;
            ucsky_prob = prob_real;
            for (int k = 0; k < l; ++ k) {
                ucsky[k] = subset[k];
            }
        }
        return;
    }
    if (j >= n) {
        return;
    }
    if (prob < ucsky_prob) {
        return;
    }
    if (conflict(subset, j)) {
        enum_l_subset(j + 1, subset, i, prob);
    } else {
        subset[i] = j;
        prob = prob * tuples[j].prob;
        // include tuples[j]
        enum_l_subset(j + 1, subset, i + 1, prob);
        // exclude tuples[j]
        subset[i] = -1;
        prob = prob / tuples[j].prob;
        enum_l_subset(j + 1, subset, i, prob);
    }
}

void UCSKY_Solver::loop_bsl() {
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        return t.sum() < s.sum();
    });
    vector<int> subset(l);
    for (int i = 0; i < l; ++ i) subset[i] = -1;
    BigFloat prob(1);
    enum_l_subset(0, subset, 0, prob);
    // cout << "end enum\n";
    // output id
    /* for (int i = 0; i < l; ++ i) {
        int index = ucsky[i];
        ucsky[i] = tuples[index].id;
    } */
}

BigFloat UCSKY_Solver::greedy() {
    vector<int> res(l, 0);
    BigFloat prob = 1;
    // construct_dg();
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        return t.sum() < s.sum();
    });

    Heap<BigFloat> h = Heap<BigFloat>();
    unordered_map<int, unordered_set<int>> inN;
    unordered_map<int, unordered_set<int>> outN;
    for (int i = 0; i < n; ++ i) {
        inN[i] = unordered_set<int>{};
        outN[i] = unordered_set<int>{};
        BigFloat beta = tuples[i].prob;
        for (int j = 0; j < i; ++ j) {
            if (tuples[j].dominate(tuples[i])) {
                outN[j].insert(i);
                inN[i].insert(j);
                beta = beta * (1 - tuples[j].prob);
            }
        }
        h.push(make_pair(beta, i));
        // cout << beta << endl;
    }
    // cout << "pre done\n";
    vector<bool> inH(n, true);
    for (int i = 0; i < l && !h.empty(); ++ i) {
        // cout << h.size() << endl;
        auto t = h.top();
        h.pop();
        inH[t.second] = false;
        res[i] = t.second;
        prob = prob * t.first;
        // cout << i << '\t' << t.second << '\t' << t.first << endl;
        for (auto s : outN[t.second]) {
            if (inH[s]) {
                inH[s] = false;
                h.remove(s);
            }
        }
        for (auto s : inN[t.second]) {
            if (inH[s]) {
                inH[s] = false;
                h.remove(s);
                for (auto w : outN[s]) {
                    if (inH[w]) {
                        h.divide(w, tuples[s].prob);
                    }
                }
            }
        }
    }
    // cout << prob << endl;
    return prob;
}

BigFloat UCSKY_Solver::alpha_greedy(const double& step) {
    // construct_dg();
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        return t.sum() < s.sum();
    });

    unordered_map<int, unordered_set<int>> inN;
    unordered_map<int, unordered_set<int>> outN;
    vector<BigFloat> beta(n, 0);
    vector<int> degree(n, 0);
    for (int i = 0; i < n; ++ i) {
        inN[i] = unordered_set<int>{};
        outN[i] = unordered_set<int>{};
        beta[i] = tuples[i].prob;
        for (int j = 0; j < i; ++ j) {
            if (tuples[j].dominate(tuples[i])) {
                outN[j].insert(i);
                ++ degree[j];
                inN[i].insert(j);
                ++ degree[i];
                beta[i] = beta[i] * (1 - tuples[j].prob);
            }
        }
    }

    for (double alpha = 0; alpha <= 1; alpha += step) {
        vector<int> res(l, 0);
        BigFloat prob = 1;
        // cout << alpha << endl;
        Heap<double> h = Heap<double>();
        for (int i = 0; i < n; ++ i) {
            double key1 = (1 - alpha)*beta[i].log();
            double key2 = alpha*degree[i];
            // cout << key1 << '\t' << key2 << '\t' << inN[i].size() + outN[i].size() << endl;
            h.push(make_pair(key1 - key2, i));
        }
        vector<bool> inH(n, true);
        int i = 0;
        for ( ; i < l && !h.empty(); ++ i) {
            // cout << h.size() << endl;
            auto t = h.top();
            h.pop();
            inH[t.second] = false;
            res[i] = t.second;
            prob = prob * tuples[t.second].prob;
            // cout << i << '\t' << t.second << '\t' << t.first << endl;
            for (auto s : outN[t.second]) {
                if (inH[s]) {
                    inH[s] = false; h.remove(s);
                    for(auto w : outN[s]) if (inH[w]) h.add(w, alpha);
                    for(auto w : inN[s]) if (inH[w]) h.add(w, alpha);
                }
            }
            for (auto s : inN[t.second]) {
                if (inH[s]) {
                    inH[s] = false; h.remove(s);
                    for (auto w : outN[s]) {
                        if (inH[w]) {
                            h.add(w, (alpha - 1.0)*log10(1 - tuples[s].prob));
                            h.add(w, alpha);
                        }
                    }
                    for (auto w : inN[s]) if (inH[w]) h.add(w, alpha);
                }
            }
        }

        if (i == l) {
            for (int i = 0; i < n; ++ i) {
                for (int j : res) {
                    if (tuples[i].dominate(tuples[j])) {
                        prob = prob * (1 - tuples[i].prob);
                        break;
                    }
                }
            }
            // cout << prob << endl;

            // for (int i = 0; i < l; ++ i) {
            //     for (int j = 0; j < i; ++ j) {
            //         if (tuples[res[i]].dominate(tuples[res[j]]) || tuples[res[j]].dominate(tuples[res[i]]))
            //             cout << "conflict\n";
            //     }
            // }

            ucsky_prob = prob;
            ucsky = res;

            return prob;
        }
    }

    return 0;
}

void UCSKY_Solver::basic_dp() {
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        if (t.coord[0] < s.coord[0]) return true;
        else if (t.coord[0] == s.coord[0]) {
            if (t.coord[1] < s.coord[1]) return true;
        }
        return false;
    }); 
    vector<unordered_map<int, BigFloat>> gamma(n);
    for (int i = 1; i < n; ++ i) {
        for (int j = 0; j < i; ++ j) {
            if (!tuples[j].dominate(tuples[i])) {
                BigFloat g(1);
                for (int k = 0; k < i; ++ k) {
                    if (tuples[k].dominate(tuples[i]) && !tuples[k].dominate(tuples[j])) {
                        g = g*(1 - tuples[k].prob);
                    }
                }
                // gamma.insert({make_pair(tuples[i].id, tuples[j].id), g});
                gamma[i].insert({j, g});
            }
        }
    }

    /* for (int i = 0; i < n; ++ i) {
        for (auto it = gamma[i].begin(); it != gamma[i].end(); ++ it) {
            cout << tuples[i].id << ", " << tuples[it->first] << ", " << it->second << endl;
        }
    } */

    // cout << "gamma done\n";

    int** Y = new int*[n];
    BigFloat** cost = new BigFloat*[n];
    for (int i = 0; i < n; ++ i) {
        Y[i] = new int[l + 1];
        cost[i] = new BigFloat[l + 1];
    }
    // dp
    for (int k = 1; k <= l; ++ k) {
        // cout << "done " << k << endl;
        if (k == 1) {
            for (int i = 0; i < n; ++ i) {
                Y[i][k] = i;
                cost[i][k] = BigFloat(tuples[i].prob);
                for (int j = 0; j < i; ++ j) {
                    if (tuples[j].dominate(tuples[i])) {
                        cost[i][k] = cost[i][k] * (1 - tuples[j].prob);
                    }
                }
            }
        } else {
            for (int i = 0; i < n; ++ i) {
                Y[i][k] = -1;
                cost[i][k] = BigFloat(0);
                for (int j = 0; j < i; ++ j) {
                    if (!tuples[j].dominate(tuples[i]) && Y[j][k - 1] != -1) {
                        BigFloat temp = tuples[i].prob;
                        // temp = temp * gamma[make_pair(tuples[i].id, tuples[j].id)];
                        temp = temp * gamma[i][j];
                        temp = temp * cost[j][k - 1];
                        if (temp > cost[i][k]) {
                            cost[i][k] = temp;
                            Y[i][k] = j;
                        }
                    }
                }
            }
        }
    }

    /* for (int i = 0; i < n; ++ i) {
        for (int k = 1; k <= l; ++ k) {
            cout << "[" << i << ", " << k << ", " << cost[i][k] << ", " << Y[i][k] << "]";
        }
        cout << endl;
    } */
    // collect result
    int index = -1;
    for (int i = 0; i < n; ++ i) {
        if (cost[i][l] > ucsky_prob) {
            ucsky_prob = cost[i][l];
            index = i;
        }
    }
    for (int i = l; i >= 1; -- i) {
        // output id
        // ucsky[i - 1] = tuples[index].id;
        ucsky[i - 1] = index;
        index = Y[index][i];
    }
    for (int i = 0; i < n; ++ i) {
        delete[] Y[i];
        delete[] cost[i];
    }
    delete[] Y;
    delete[] cost;
}

void UCSKY_Solver::group_dp() {
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        if (t.coord[0] < s.coord[0]) return true;
        else if (t.coord[0] == s.coord[0]) {
            if (t.coord[1] < s.coord[1]) return true;
        }
        return false;
    });
    // group tuples
    vector<vector<int>> groups;
    groups.push_back(vector<int>{0});
    vector<double> group2x;
    group2x.push_back(tuples[0].coord[0]);
    for (int i = 1; i < n; ++ i) {
        if (tuples[i].coord[0] != tuples[i - 1].coord[0]) {
            groups.push_back(vector<int>{i});
            group2x.push_back(tuples[i].coord[0]);
        } else {
            groups.back().push_back(i);
        }
    }
    struct aux{
        double y;
        int index1;
        int index2;
    };
    vector<vector<aux>> A(groups.size());
    for (int q = 0; q < groups.size(); ++ q) {
        reverse(groups[q].begin(), groups[q].end());
        A[q].resize(groups[q].size());
        for (int i = 0; i < groups[q].size(); ++ i) {
            A[q][i] = {tuples[groups[q][i]].coord[1], -1, -1};
        }
    }
    
    vector<unordered_map<int, BigFloat>> gamma(n);
    for (int q = 0; q < groups.size(); ++ q) {
        for (int i : groups[q]) {
            for (int p = 0; p < q; ++ p) {
                BigFloat g(1);
                Tuple virtual_t(-1, 2, nullptr, 0);
                virtual_t.coord[0] = group2x[p];
                virtual_t.coord[1] = tuples[i].coord[1];
                for (int k = 0; k < i; ++ k) {
                    if (tuples[k].dominate(tuples[i]) && !tuples[k].dominate(virtual_t)) {
                        g = g*(1 - tuples[k].prob);
                    }
                }
                gamma[i].insert({p, g});
            }
        }
    }
    // cout << "gamma rdy.\n";
    /* for (int q = 0; q < groups.size(); ++ q) {
        for (int i : groups[q]) {
            for (int p = 0; p < q; ++ p) {
                cout << tuples[i].id << ", " << group2x[p] << ", " << gamma[i][p] << endl;
            }
        }
    } */

    int** Y = new int*[n];
    BigFloat** cost = new BigFloat*[n];
    for (int i = 0; i < n; ++ i) {
        Y[i] = new int[l + 1];
        cost[i] = new BigFloat[l + 1];
    }

    for (int k = 1; k <= l; ++ k) {
        for (int q = 0; q < groups.size(); ++ q) {
            for (int i = 0; i < groups[q].size(); ++ i) {
                int t = groups[q][i]; // index of current in tuples
                if (k == 1) {
                    Y[t][k] = t;
                    cost[t][k] = BigFloat(tuples[t].prob);
                    for (int j = 0; j < t; ++ j) {
                        if (tuples[j].dominate(tuples[t])) {
                            cost[t][k] = cost[t][k] * (1 - tuples[j].prob);
                        }
                    }
                } else {
                    Y[t][k] = -1;
                    cost[t][k] = BigFloat(0);
                    for (int p = 0; p < q; ++ p) {
                        int l = 0, r = A[p].size();
                        int index = -1;
                        while (l < r) {
                            index = (l + r)/2;
                            if (A[p][index].y <= tuples[t].coord[1]) {
                                r = index;
                            } else {
                                l = index + 1;
                            }
                        }
                        -- l;
                        if (l >= 0 && A[p][l].y > tuples[t].coord[1]) {
                            int j = k%2 == 0 ? A[p][l].index2 : A[p][l].index1;
                            if (Y[j][k - 1] != -1) {
                                BigFloat temp = tuples[t].prob;
                                temp = temp * gamma[t][p];
                                temp = temp * cost[j][k - 1];    
                                if (temp > cost[t][k]) {
                                    cost[t][k] = temp;
                                    Y[t][k] = j;
                                }
                            }
                        }
                    }
                }
                // update A
                if (k%2 == 0) {
                    if (i == 0) { A[q][i].index1 = t; }
                    else {
                        if (cost[t][k] > cost[A[q][i - 1].index1][k]) {
                            A[q][i].index1 = t;
                        } else {
                            A[q][i].index1 = A[q][i - 1].index1;
                        }
                    }
                } else {
                    if (i == 0) { A[q][i].index2 = t; }
                    else {
                        if (cost[t][k] > cost[A[q][i - 1].index2][k]) {
                            A[q][i].index2 = t;
                        } else {
                            A[q][i].index2 = A[q][i - 1].index2;
                        }
                    }
                }
            }
        }
    }

    /* for (int i = 0; i < n; ++ i) {
        for (int k = 1; k <= l; ++ k) {
            cout << "[" << i << ", " << k << ", " << cost[i][k] << ", " << Y[i][k] << "]";
        }
        cout << endl;
    } */
    

    int index = -1;
    for (int i = 0; i < n; ++ i) {
        if (cost[i][l] > ucsky_prob) {
            ucsky_prob = cost[i][l];
            index = i;
        }
    }
    for (int i = l; i >= 1; -- i) {
        // output id
        // ucsky[i - 1] = tuples[index].id;
        ucsky[i - 1] = index;
        index = Y[index][i];
    }
    for (int i = 0; i < n; ++ i) {
        delete[] Y[i];
        delete[] cost[i];
    }
    delete[] Y;
    delete[] cost;
}

void UCSKY_Solver::sample_fim(const int& theta) {
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        return t.sum() < s.sum();
    });

    vector<int> L0;
    unordered_map<int, unordered_set<int>> inN;
    unordered_map<int, unordered_set<int>> outN;
    for (int i = 0; i < n; ++ i) {
        if (!dominate(L0, i)) {
            L0.push_back(i);
        }
        inN[i] = unordered_set<int>{};
        outN[i] = unordered_set<int>{};
        for (int j = 0; j < i; ++ j) {
            if (tuples[j].dominate(tuples[i])) {
                outN[j].insert(i);
                inN[i].insert(j);
            }
        }
    }

    vector<vector<int>> S(theta, vector<int>{});
    unordered_map<int, int> freq;
    bool flag = false;
    for (int i = 0; i < theta; ++ i) {
        vector<int> degree(n, 0);
        for (int j = 0; j < n; ++ j) degree[j] = inN[j].size();
        vector<int> C = L0;
        while (!C.empty()) {
            int t = C.back();
            C.pop_back();
            double p = rand_uniform(0, 1);
            if (p <= tuples[t].prob) {
                S[i].push_back(t);
            } else {
                for (auto s : outN[t]) {
                    -- degree[s];
                    if (degree[s] == 0) {
                        C.push_back(s);
                    }
                }
            }
        }
        for (auto t : S[i]) {
            freq[t] ++;
        }
        if (S[i].size() >= l) flag = true;
        /* cout << S[i].size() << endl;
        for (auto t : S[i]) {
            for (auto s : S[i]) {
                if (t != s) {
                    if (tuples[t].dominate(tuples[s]) || tuples[s].dominate(tuples[t])) {
                        cout << "conflict\n";
                    }
                }
            }
        } */
    }
    if (!flag) {
        cout << "no pattern with length " << l << endl;
        return;
    }
    // mining frequent itemset of size l
    int min_sup = theta;
    for (auto iter = freq.begin(); iter != freq.end(); ++ iter) {
        min_sup = min(min_sup, iter->second);
    }
    bool found = false;
    while (!found) {
        FPTree fptree(S, min_sup);
        set<Pattern> patterns = fptree_growth(fptree);
        int max_freq = 0;
        unordered_map<int, vector<set<int>>> l_set;
        // cout << patterns.size() << endl;
        for (auto p : patterns) {
            if (p.first.size() == l) {
                found = true;
                max_freq = max(max_freq, (int)p.second);
                l_set[p.second].push_back(p.first);
            }
        }
        if (found) {
            // cout << l_set[max_freq].size() << endl;
            for (auto s : l_set[max_freq]) {
                BigFloat prob = compute_prob(s);
                if (prob > ucsky_prob) {
                    ucsky_prob = prob;
                    int j = 0;
                    for (auto k : s) {
                        ucsky[j ++] = k;
                    }
                }
            }
            return;
        }
        min_sup /= 2;
    }
}

void UCSKY_Solver::branch_bound(const double& step) {
    // double step = 0.002;
    BigFloat lb = alpha_greedy(step);

    vector<pair<double, int>> prob_id;
    prob_id.resize(n);
    for (int i = 0; i < n; ++ i) {
        prob_id[i] = make_pair(tuples[i].prob, i);
    }
    sort(prob_id.begin(), prob_id.end(), [](const pair<double, int>& t, const pair<double, int>& s){
        return t.first > s.first;
    });

    vector<int> reduced_D;
    // vector<pair<int, double>> reduced_D;
    if (lb > 0) { // preprocess
        for (int i = 0; i < n; ++ i) {
            BigFloat beta = tuples[i].prob;
            for (int j = 0; j < i; ++ j) {
                if (tuples[j].dominate(tuples[i])) {
                    beta = beta * (1 - tuples[j].prob);
                }
            }
            if (beta > lb) {
                int j = l;
                int k = 0;
                while (j > 0 && k < n) {
                    int id = prob_id[k].second;
                    if (!tuples[id].dominate(tuples[i]) && !tuples[i].dominate(tuples[id])) {
                        beta = beta * tuples[id].prob;
                        -- j;
                    }
                    ++ k;
                }
                if (beta > lb) {
                    reduced_D.push_back(i);
                    // reduced_D.push_back(make_pair(i, tuples[i].sum()));
                }
            }
        }
    }
    // lb maintained by ucsky_prob
    // cout << reduced_D.size() << endl;
    vector<int> S;
    BigFloat prob(0);
    bb(reduced_D, S, prob, l);
    /* vector<pair<int, double>> sky;
    cout << reduced_D.size() << endl;
    skyline(reduced_D, sky);
    for (auto t : sky) {
        cout << tuples[t.first] << endl;
    }
    cout << reduced_D.size() << '\t' << sky.size() << endl; */
    // divide(reduced_D);
}

void UCSKY_Solver::divide(vector<pair<int, double>>& T) {
    vector<int> sky;
    skyline(T, sky);
    for (int i = 0; i < min(l, (int)sky.size()); ++ i) {
        if (i == 0 && T.size() >= l) divide(T);
        else {
            vector<int> subset;
            enum_l_subset(T, sky, i, 0, subset, 0);
        }
    }
}

void UCSKY_Solver::skyline(vector<pair<int, double>>& T, vector<int>& sky) {
    sort(T.begin(), T.end(), [](const pair<int, double>& t, const pair<int, double>& s){
        return t.second < s.second;
    });
    vector<int> index;
    for (int i = 0; i < T.size(); ++ i) {
        bool flag = false;
        for (int j = 0; j < sky.size(); ++ j) {
            if (tuples[sky[j]].dominate(tuples[T[i].first])) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            sky.push_back(T[i].first);
            index.push_back(i);
        }
    }
    for (int i : index) {
        T[i] = T.back();
        T.pop_back();
    }
}

void UCSKY_Solver::enum_l_subset(const vector<pair<int, double>>& T, const vector<int>& sky, int s, int j, vector<int>& subset, int i) {
    if (i == s) {
        vector<int> new_T;
        for (auto t : T) {
            for (auto s : subset) {
                if (!tuples[t.first].dominate(tuples[s]) && !tuples[s].dominate(tuples[t.first])) {
                    new_T.push_back(t.first);
                }
            }
        }
        BigFloat prob = compute_prob(subset);
        bb(new_T, subset, prob, l - s);
    }
    if (j >= sky.size()) return;
    subset[i] = sky[j];
    enum_l_subset(T, sky, s, j + 1, subset, i + 1);
    enum_l_subset(T, sky, s, j + 1, subset, i);
}

void UCSKY_Solver::bb(const vector<int>& T, const vector<int>& S, const BigFloat& prob, const int& k) {
    if (k == 0) {
        if (prob > ucsky_prob) {
            ucsky_prob = prob;
            ucsky = S;
        }
        return;
    }
    if (T.size() < k) return;
    if (prob < ucsky_prob) return;
    
    vector<double> prob_T;
    prob_T.resize(T.size());
    for (int i = 0; i < T.size(); ++ i) {
        prob_T[i] = tuples[T[i]].prob;
    }
    sort(prob_T.begin(), prob_T.end(), [](const double& t, const double& s){
        return t > s;
    });
    BigFloat beta = prob;
    for (int i = 0; i < k; ++ i) {
        beta = beta * prob_T[i];
    }
    if (beta < ucsky_prob) return;

    for (int i = 0; i < T.size(); ++ i) {
        vector<int> new_T;
        for (int j = i + 1; j < T.size(); ++ j) {
            if (!tuples[T[i]].dominate(tuples[T[j]]) && !tuples[T[j]].dominate(tuples[T[i]])) {
                new_T.push_back(T[j]);
            }
        }
        vector<int> new_S = S;
        new_S.push_back(T[i]);
        BigFloat new_S_prob = compute_prob(new_S);
        bb(new_T, new_S, new_S_prob, k - 1);
    }
}

void UCSKY_Solver::print_ucsky() {
    cout << "UCSKY: ";
    for (int i : ucsky) {
        cout << tuples[i].id << ", ";
        // cout << tuples[i] << endl;
    }
    cout << "\t Prob: " << ucsky_prob << endl;
}

void UCSKY_Solver::check_prob() {
    BigFloat prob(1);
    for (int i = 0; i < n; ++ i) {
        for (int j : ucsky) {
            if (tuples[i].dominate(tuples[j])) {
                prob = prob * (1 - tuples[i].prob);
                break;
            }
        }
    }
    for (int j : ucsky) {
        prob = prob * tuples[j].prob;
    }
    cout << prob << '\t' << ucsky_prob << endl;
}


BigFloat UCSKY_Solver::compute_prob(const set<int>& s) {
    BigFloat prob(1);
    for (int i = 0; i < n; ++ i) {
        for (int j : s) {
            if (tuples[i].dominate(tuples[j])) {
                prob = prob * (1 - tuples[i].prob);
                break;
            }
        }
    }
    for (int j : s) {
        prob = prob * tuples[j].prob;
    }
    return prob;
}

BigFloat UCSKY_Solver::compute_prob(const vector<int>& s) {
    BigFloat prob(1);
    for (int i = 0; i < n; ++ i) {
        for (int j : s) {
            if (tuples[i].dominate(tuples[j])) {
                prob = prob * (1 - tuples[i].prob);
                break;
            }
        }
    }
    for (int j : s) {
        prob = prob * tuples[j].prob;
    }
    return prob;
}
