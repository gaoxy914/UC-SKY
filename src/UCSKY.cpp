#include "UCSKY.h"

UCSKY_Solver::UCSKY_Solver(const int& dim, const int& n, const double& center, const int& l) {
    srand((unsigned)time(nullptr));
    srand48((unsigned)time(nullptr));
    this->dim = dim;
    this->n = n;
    this->center = center;
    this->l = l;
    this->tuples.resize(n, Tuple(-1, dim, nullptr, 0, false));
    sorted = false;
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
    for (int i = 0; i < n; ++ i) {
        tuples[i].id = i;
        for (int j = 0; j < dim; ++ j) {
            tuples[i].coord[j] = rand_uniform(0, 1);
        }
        tuples[i].prob = rand_uniform(lower, upper);
    }
}

void UCSKY_Solver::gen_ind_data_2d(const int& m) {
    double delta = min(center, 1 - center);
    double lower = center - delta, upper = center + delta;
    vector<double> x(m, 0);
    for (int i = 0; i < m; ++ i) {
        x[i] = rand_uniform(0, 1);
    }
    for (int i = 0; i < n; ++ i) {
        tuples[i].id = i;
        tuples[i].coord[0] = x[i%m];
        tuples[i].coord[1] = rand_uniform(0, 1);
        tuples[i].prob = rand_uniform(lower, upper);
    }
}

void UCSKY_Solver::gen_anti_data() { // fix generation
    double delta = min(center, 1 - center);
    double lower = center - delta, upper = center + delta;
    double rho = 0.5;
    for (int i = 0; i < n; ++ i) {
        tuples[i].id = i;
        // double range = 0.5*dim + rand_normal(0, 0.05);
        // tuples[i].coord[0] = rand_uniform(0, 1)*min(1.0, range);
        // for (int j = 1; j < dim; ++ j) {
            // range -= tuples[i].coord[j - 1];
            // tuples[i].coord[j] = rand_uniform(0, 1)*min(1.0, range);
            // if (j == dim - 1) tuples[i].coord[j] = min(1.0, range);
        // }
        tuples[i].coord[0] = rand_uniform(0, 1);
        for (int j = 1; j < dim; ++ j) {
            tuples[i].coord[j] = roundoff((rho*(1 - tuples[i].coord[j - 1]) + sqrt(1 - rho*rho)*rand_uniform(0, 1))/(rho + sqrt(1 - rho*rho)), PREC);
        }
        tuples[i].prob = rand_uniform(lower, upper);
    }
}

void UCSKY_Solver::gen_corr_data() {// fix generation
    double delta = min(center, 1 - center);
    double lower = center - delta, upper = center + delta;
    double rho = 0.5;
    for (int i = 0; i < n; ++ i) {
        tuples[i].id = i;
        // tuples[i].coord[0] = rand_uniform(0, 1);
        // for (int j = 1; j < dim; ++ j) {
            // do {
                // tuples[i].coord[j] = rand_normal(tuples[i].coord[0], 0.05);
            // } while (tuples[i].coord[j] < 0 || tuples[i].coord[j] > 1);
        // }
        tuples[i].coord[0] = rand_uniform(0, 1);
        for (int j = 1; j < dim; ++ j) {
            tuples[i].coord[j] = roundoff((rho*tuples[i].coord[j - 1] + sqrt(1 - rho*rho)*rand_uniform(0, 1))/(rho + sqrt(1 - rho*rho)), PREC);
        }
        tuples[i].prob = rand_uniform(lower, upper);
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
        ifstream file("data/car.dat", ios::in);
        if (!file.is_open()) {
            printf("Fail in opening data file.\n");
            exit(1);
        }
        for (int i = 0; i < n; ++ i) {
            tuples[i].id = i;
            file >> tuples[i].coord[0];
            file >> tuples[i].coord[1];
            file >> tuples[i].coord[2];
            file >> tuples[i].coord[3];
            file >> tuples[i].prob;
        }
        file.close();
    } else if (strcmp(data_path, "data/nba.dat") == 0) {
        ifstream file("data/nba.dat", ios::in);
        if (!file.is_open()) {
            printf("Fail in opening data file.\n");
            exit(1);
        }
        for (int i = 0; i < n; ++ i) {
            tuples[i].id = i;
            int j = 0;
            double temp = 0;
            for ( ; j < dim; ++ j) {
                file >> tuples[i].coord[j];
            }
            for ( ; j < 8; ++ j) {
                file >> temp;
            }
            file >> tuples[i].prob;
        }
        file.close();
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

void UCSKY_Solver::loop_bsl() {
    if (!sorted) {
        sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
            return t.sum() < s.sum();
        });
        sorted = true;
    }
    vector<int> subset(l);
    for (int i = 0; i < l; ++ i) subset[i] = -1;
    BigFloat prob(1);
    enum_l_subset(0, subset, 0, prob);
}

void UCSKY_Solver::enum_l_subset(int j, vector<int>& subset, int i, BigFloat& prob) {
    if (i == l) {
        BigFloat prob_real = prob;
        for (int p = 0; p < n; ++ p) {
            for (int q : subset) {
                if (tuples[p].dominate(tuples[q])) {
                    prob_real = prob_real * (1 - tuples[p].prob);
                    break;
                }
            }
            if (prob_real < ucsky_prob) return;
        }
        if (prob_real > ucsky_prob) {
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
    if (conflict(subset, j) || tuples[j].reduced) {
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

bool UCSKY_Solver::conflict(const vector<int>& S, const int& t) {
    for (auto s : S) {
        if (s != -1 && 
        (tuples[s].dominate(tuples[t]) || tuples[t].dominate(tuples[s]))) {
            return true;
        }
    }
    return false;
}

bool UCSKY_Solver::greedy() {
    if (!sorted) {
        sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
            return t.sum() < s.sum();
        });
        sorted = true;
    }
    Heap<BigFloat> h = Heap<BigFloat>(n);
    inN.resize(n);
    outN.resize(n);
    for (int i = 0; i < n; ++ i) {
        inN[i] = vector<int>{};
        outN[i] = vector<int>{};
        BigFloat beta = tuples[i].prob;
        for (int j = 0; j < i; ++ j) {
            if (tuples[j].dominate(tuples[i])) {
                outN[j].push_back(i);
                inN[i].push_back(j);
                beta = beta * (1 - tuples[j].prob);
            }
        }
        h.push(make_pair(beta, i));
    }
    vector<bool> inH(n, true);
    int i = 0;
    for ( ; i < l && !h.empty(); ++ i) {
        // cout << h.size() << endl;
        auto t = h.top();
        h.pop();
        inH[t.second] = false;
        ucsky[i] = t.second;
        ucsky_prob = ucsky_prob * t.first;
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
    if (i < l) {
        return false;
    }
    return true;
}

// collect in-time
bool UCSKY_Solver::greedy_plus(const double& step) {
    if (!sorted) {
        sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
            return t.sum() < s.sum();
        });
        sorted = true;
    }
    // inN.resize(n);
    // outN.resize(n);
    beta.resize(n);
    vector<int> degree(n, 0);
    for (int i = 0; i < n; ++ i) {
        // inN[i] = vector<int>{};
        // outN[i] = vector<int>{};
        beta[i] = tuples[i].prob;
        for (int j = 0; j < i; ++ j) {
            if (tuples[j].dominate(tuples[i])) {
                // outN[j].push_back(i);
                // inN[i].push_back(j);
                ++ degree[i];
                ++ degree[j];
                beta[i] = beta[i] * (1 - tuples[j].prob);
            }
        }
    }
    for (double alpha = 0; alpha <= 1; alpha += step) {
        ucsky_prob = 1;
        Heap<double> h = Heap<double>(n);
        for (int i = 0; i < n; ++ i) {
            double key1 = (1 - alpha)*beta[i].log();
            // double key2 = alpha*(inN[i].size() + outN[i].size());
            double key2 = alpha*degree[i];
            h.push(make_pair(key1 - key2, i));
        }
        // vector<bool> inH(n, true);
        int i = 0;
        for ( ; i < l && !h.empty(); ++ i) {
            auto t = h.top();
            h.pop();
            // inH[t.second] = false;
            ucsky[i] = t.second;
            ucsky_prob = ucsky_prob * tuples[t.second].prob;
            
            /* for (auto s : outN[t.second]) {
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
            } */

            vector<int> dom;
            vector<int> bedom;
            for (auto s : h.heap) {
                if (tuples[s.second].dominate(tuples[t.second])) {
                    bedom.push_back(s.second);
                } else if (tuples[t.second].dominate(tuples[s.second])) {
                    dom.push_back(s.second);
                }
            }
            for (auto s : dom) h.remove(s);
            for (auto s : bedom) h.remove(s);
            for (auto s : dom) {
                for (auto w : h.heap) {
                    if (tuples[s].dominate(tuples[w.second]) || tuples[w.second].dominate(tuples[s])) {
                        h.add(w.second, alpha);
                    }
                }
            }
            for (auto s : bedom) {
                for (auto w : h.heap) {
                    if (tuples[s].dominate(tuples[w.second])) {
                        h.add(w.second, (alpha - 1.0)*log10(1 - tuples[s].prob));
                        h.add(w.second, alpha);
                    } else if (tuples[w.second].dominate(tuples[s])) {
                        h.add(w.second, alpha);
                    }
                }
            }

        } 
        
        if (i == l) {
            for (int j = 0; j < n; ++ j) {
                for (int k : ucsky) {
                    if (tuples[j].dominate(tuples[k])) {
                        ucsky_prob = ucsky_prob * (1 - tuples[j].prob);
                        break;
                    }
                }
            }
            return true;
        }
    }
    // cout << "fail to find valid solution.\n";
    return false;
}

int UCSKY_Solver::greedy_plus2(const double& step) {
    if (!sorted) {
        sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
            return t.sum() < s.sum();
        });
        sorted = true;
    }
    vector<int> reduced_tuples;
    for (int i = 0; i < n; ++ i) {
        if (!tuples[i].reduced) reduced_tuples.push_back(i);
    }
    assert(reduced_tuples.size() == nR);
    int m = nR < n ? nR : n;
    inN.resize(m);
    outN.resize(m);
    // degree.resize(m);
    for (int i = 0; i < m; ++ i) {
        inN[i] = vector<int>{};
        outN[i] = vector<int>{};
        for (int j = 0; j < i; ++ j) {
            if (tuples[reduced_tuples[j]].dominate(tuples[reduced_tuples[i]])) {
                outN[j].push_back(i);
                // ++ degree[j];
                inN[i].push_back(j);
                // ++ degree[i];
            }
        }
    }
    for (double alpha = 0; alpha <= 1; alpha += step) {
        BigFloat ucsky_prob_temp = 1;
        vector<int> ucsky_temp(l);
        Heap<double> h = Heap<double>(m);
        for (int i = 0; i < m; ++ i) {
            double key1 = (1 - alpha)*beta[reduced_tuples[i]].log();
            double key2 = alpha*(inN[i].size() + outN[i].size());
            h.push(make_pair(key1 - key2, i));
        }
        vector<bool> inH(m, true);
        int i = 0;
        for ( ; i < l && !h.empty(); ++ i) {
            auto t = h.top();
            h.pop();
            inH[t.second] = false;
            ucsky_temp[i] = reduced_tuples[t.second];
            ucsky_prob_temp = ucsky_prob_temp * tuples[reduced_tuples[t.second]].prob;
            for (auto s : outN[t.second]) {
                if (inH[s]) {
                    inH[s] = false; h.remove(s);
                    for (auto w : outN[s]) if (inH[w]) h.add(w, alpha);
                    for (auto w : inN[s]) if (inH[w]) h.add(w, alpha);
                }
            }
            for (auto s : inN[t.second]) {
                if (inH[s]) {
                    inH[s] = false; h.remove(s);
                    for (auto w : outN[s]) {
                        if (inH[w]) {
                            h.add(w, (alpha - 1.0)*log10(1 - tuples[reduced_tuples[s]].prob));
                            h.add(w, alpha);
                        }
                    }
                    for (auto w : inN[s]) if (inH[w]) h.add(w, alpha);
                }
            }
        }
        if (i == l) {
            for (int j = 0; j < n; ++ j) {
                for (int k : ucsky_temp) {
                    if (tuples[j].dominate(tuples[k])) {
                        ucsky_prob_temp = ucsky_prob_temp * (1 - tuples[j].prob);
                        break;
                    }
                }
            }
            if (ucsky_prob_temp > ucsky_prob) {
                ucsky_prob = ucsky_prob_temp;
                for (int j = 0; j < l; ++ j) ucsky[j] = ucsky_temp[j];
                return 2;
            }
            return 1;
        }
    }
    // cout << "fail to find valid solution.\n";
    return 0;
}

bool UCSKY_Solver::sky_2n() {
    if (!sorted) {
        sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
            return t.sum() < s.sum();
        });
    }
    vector<pair<int, double>> sky;
    bool flag = false;
    for (int i = 0; i < n; ++ i) {
        if (tuples[i].prob > 0.5) {
            flag = false;
            for (auto j : sky) {
                if (tuples[j.first].dominate(tuples[i])) {
                    flag = true;
                    break;
                }
            }
            if (!flag) {
                sky.push_back(make_pair(i, tuples[i].prob));
            }
        }
    }
    cout << sky.size() << endl;
    if (sky.size() < l) {
        // cout << "fail to find valid solution.\n";
        return false;
    }
    sort(sky.begin(), sky.end(), [](const pair<int, double>& t, const pair<int, double>& s){
        return t.second > s.second;
    });
    ucsky_prob = 1;
    for (int i = 0; i < l; ++ i) {
        ucsky[i] = sky[i].first;
        ucsky_prob = ucsky_prob * sky[i].second;
    }
    for (int i = 0; i < n; ++ i) {
        for (int j : ucsky) {
            if (tuples[i].dominate(tuples[j])) {
                ucsky_prob = ucsky_prob * (1 - tuples[i].prob);
                break;
            }
        }
    }
    return true;
}

int UCSKY_Solver::reduce_data(const double& step) {
    nR = n;
    beta.resize(n);
    if (sky_2n()) {
        for (int i = 0; i < n; ++ i) {
            tuples[i].reduced = false;
            beta[i] = tuples[i].prob;
            for (int j = 0; j < i; ++ j) {
                if (tuples[j].dominate(tuples[i])) {
                    beta[i] = beta[i] * (1 - tuples[j].prob);
                    if (beta[i] < ucsky_prob) {
                        tuples[i].reduced = true;
                        -- nR;
                        break;
                    }
                }
            }
        }
        // cout << ucsky_prob << '\t' << nR << endl;
        if (greedy_plus2(step) == 2) {
            vector<pair<double, int>> prob_id;
            prob_id.resize(nR);
            int j = 0;
            for (int i = 0; i < n; ++ i) {
                if (!tuples[i].reduced) {
                    prob_id[j ++] = make_pair(tuples[i].prob, i);
                }
            }
            sort(prob_id.begin(), prob_id.end(), [](const pair<double, int>& t, const pair<double, int>& s){
                return t.first > s.first;
            });
            for (int i = 0; i < n; ++ i) {
                if (!tuples[i].reduced) {
                    if (beta[i] < ucsky_prob) {
                        tuples[i].reduced = true;
                        -- nR;
                    } else {
                        int j = l - 1, k = 0;
                        BigFloat temp = beta[i];
                        while (j > 0 && k < nR) {
                            int id = prob_id[k].second;
                            if (!tuples[id].dominate(tuples[i]) && !tuples[i].dominate(tuples[id])) {
                                temp = temp * tuples[id].prob;
                                -- j;
                            }
                            ++ k;
                        }
                        if (temp < ucsky_prob) {
                            tuples[i].reduced = true;
                            -- nR;
                        }
                    }
                }
            }
        }   
    } else if (greedy_plus(step)) {
        vector<pair<double, int>> prob_id;
        prob_id.resize(n);
        for (int i = 0; i < n; ++ i) {
            prob_id[i] = make_pair(tuples[i].prob, i);
        }
        sort(prob_id.begin(), prob_id.end(), [](const pair<double, int>& t, const pair<double, int>& s){
            return t.first > s.first;
        });
        nR = 0;
        for (int i = 0; i < n; ++ i) {
            tuples[i].reduced = true;
            if (beta[i] > ucsky_prob) {
                int j = l, k = 0;
                BigFloat temp = beta[i];
                while (j > 0 && k < n) {
                    int id = prob_id[k].second;
                    if (!tuples[id].dominate(tuples[i]) && !tuples[i].dominate(tuples[id])) {
                        temp = temp * tuples[id].prob;
                        -- j;
                    }
                    ++ k;
                }
                if (temp > ucsky_prob) {
                    ++ nR;
                    // reduced_tuples.push_back(tuples[i]);
                    tuples[i].reduced = false;
                }
            }
        }
    }
    return nR;
}

int UCSKY_Solver::reduce_data_bsl() {
    if (!sorted) {
        sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
            return t.sum() < s.sum();
        });
        sorted = true;
    }
    beta.resize(n);
    vector<pair<BigFloat, int>> temp(n);
    for (int i = 0; i < n; ++ i) {
        tuples[i].reduced = false;
        beta[i] = tuples[i].prob;
        for (int j = 0; j < i; ++ j) {
            if (tuples[j].dominate(tuples[i])) {
                beta[i] = beta[i] * (1 - tuples[j].prob);
            }
        }
        temp[i] = make_pair(beta[i], i);
    }
    sort(temp.begin(), temp.end(), [](const pair<BigFloat, int>& p1, const pair<BigFloat, int>& p2){
        return p1.first > p2.first;
    });
    vector<int> no_conflic;
    int i = 0;
    vector<int> topl;
    while (i < n && topl.size() < l) {
        bool flag = false;
        for (auto j : topl) {
            if (tuples[j].dominate(tuples[temp[i].second]) || tuples[temp[i].second].dominate(tuples[j])) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            topl.push_back(temp[i].second);
        }
        ++ i;
    }
    // for (int i : topl) {
        // cout << tuples[i] << endl;
    // }
    if (topl.size() < l) {
        return n;
    }
    ucsky_prob = 1;
    for (int i = 0; i < l; ++ i) {
        ucsky[i] = topl[i];
        ucsky_prob = ucsky_prob * tuples[topl[i]].prob;
    }
    for (int i = 0; i < n; ++ i) {
        for (int j : ucsky) {
            if (tuples[i].dominate(tuples[j])) {
                ucsky_prob = ucsky_prob * (1 - tuples[i].prob);
                break;
            }
        }
    }
    int nR = n;
    for (int i = 0; i < n; ++ i) {
        if (beta[i] < ucsky_prob) {
            tuples[i].reduced = true;
            -- nR;
        }
    }
    return nR;
}

void UCSKY_Solver::basic_dp() {
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        if (t.coord[0] < s.coord[0]) return true;
        else if (t.coord[0] == s.coord[0]) {
            if (t.coord[1] < s.coord[1]) return true;
        }
        return false;
    });
    vector<int> reduced_tuples;
    for (int i = 0; i < n; ++ i) {
        if (!tuples[i].reduced) reduced_tuples.push_back(i);
    }
    assert(reduced_tuples.size() == nR);
    int m = nR < n ? nR : n;
    BigFloat** gamma = new BigFloat*[m];
    for (int i = 0; i < m; ++ i) gamma[i] = new BigFloat[m];
    for (int i = 0; i < m; ++ i) {
        for (int j = 0; j < i; ++ j) {
            int t = reduced_tuples[i], s = reduced_tuples[j];
            if (!tuples[s].dominate(tuples[t])) {
                gamma[i][j] = 1;
                for (int k = s + 1; k < t; ++ k) {
                    if (tuples[k].dominate(tuples[t])) {
                        gamma[i][j] = gamma[i][j] * (1 - tuples[k].prob);
                    }
                }
            }
        }
    }
    // cout << "precompute gamma\n";
    int** Y = new int*[m];
    BigFloat** cost = new BigFloat*[m];
    for (int i = 0; i < m; ++ i) {
        Y[i] = new int[l + 1];
        cost[i] = new BigFloat[l + 1];
    }
    for (int k = 1; k <= l; ++ k) {
        if (k == 1) {
            for (int i = 0; i < m; ++ i) {
                int t = reduced_tuples[i];
                Y[i][k] = i;
                cost[i][k] = BigFloat(tuples[t].prob);
                for (int j = 0; j < t; ++ j) {
                    if (tuples[j].dominate(tuples[t])) {
                        cost[i][k] = cost[i][k] * (1 - tuples[j].prob);
                    }
                }
            }
        } else {
            for (int i = 0; i < m; ++ i) {
                Y[i][k] = -1;
                cost[i][k] = BigFloat(0);
                int t = reduced_tuples[i];
                for (int j = 0; j < i; ++ j) {
                    int s = reduced_tuples[j];
                    if (!tuples[s].dominate(tuples[t]) && Y[j][k - 1] != -1) {
                        BigFloat temp = tuples[t].prob;
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
    // cout << "dp computation\n";
    int index = -1;
    for (int i = 0; i < m; ++ i) {
        if (cost[i][l] >= ucsky_prob) {
            ucsky_prob = cost[i][l];
            index = i;
        }
    }
    for (int i = l; i >= 1; -- i) {
        ucsky[i - 1] = reduced_tuples[index];
        index = Y[index][i];
    }
    for (int i = 0; i < m; ++ i) {
        delete[] Y[i];
        delete[] cost[i];
        delete[] gamma[i];
    }
    delete[] Y;
    delete[] cost;
    delete[] gamma;
}

void UCSKY_Solver::basic_dp2() {
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        if (t.coord[0] < s.coord[0]) return true;
        else if (t.coord[0] == s.coord[0]) {
            if (t.coord[1] < s.coord[1]) return true;
        }
        return false;
    });
    BigFloat** gamma = new BigFloat*[n];
    for (int i = 0; i < n; ++ i) gamma[i] = new BigFloat[i + 1];
    for (int i = 0; i < n; ++ i) {
        for (int j = 0; j < i; ++ j) {
            if (!tuples[j].dominate(tuples[i])) {
                gamma[i][j] = 1;
                for (int k = j + 1; k < i; ++ k) {
                    if (tuples[k].dominate(tuples[i])) {
                        gamma[i][j] = gamma[i][j] * (1 - tuples[k].prob);
                    }
                }
            }
        }
    }
    // cout << "precompute gamma\n";
    int** Y = new int*[n];
    BigFloat** cost = new BigFloat*[n];
    for (int i = 0; i < n; ++ i) {
        Y[i] = new int[l + 1];
        cost[i] = new BigFloat[l + 1];
    }
    for (int k = 1; k <= l; ++ k) {
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
    // cout << "dp computation\n";
    int index = -1;
    for (int i = 0; i < n; ++ i) {
        if (cost[i][l] >= ucsky_prob) {
            ucsky_prob = cost[i][l];
            index = i;
        }
    }
    for (int i = l; i >= 1; -- i) {
        ucsky[i - 1] = index;
        index = Y[index][i];
    }
    for (int i = 0; i < n; ++ i) {
        delete[] Y[i];
        delete[] cost[i];
        delete[] gamma[i];
    }
    delete[] Y;
    delete[] cost;
    delete[] gamma;
}

void UCSKY_Solver::basic_dp3() {
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        if (t.coord[0] < s.coord[0]) return true;
        else if (t.coord[0] == s.coord[0]) {
            if (t.coord[1] < s.coord[1]) return true;
        }
        return false;
    });
    // BigFloat* gamma = new BigFloat[n];
    int** Y = new int*[n];
    BigFloat** cost = new BigFloat*[n];
    for (int i = 0; i < n; ++ i) {
        Y[i] = new int[l + 1];
        cost[i] = new BigFloat[l + 1];
    }
    for (int i = 0; i < n; ++ i) {
        // for (int k = 1; k <= l; ++ k) {
        //     if (k == 1) {
        //         Y[i][k] = i;
        //         cost[i][k] = BigFloat(tuples[i].prob);
        //         for (int j = 0; j < i; ++ j) {
        //             if (tuples[j].dominate(tuples[i])) {
        //                 cost[i][k] = cost[i][k] * (1 - tuples[j].prob);
        //             }
        //         }
        //     } else {
        //         Y[i][k] = -1;
        //         cost[i][k] = BigFloat(0);
        //         for (int j = 0; j < i; ++ j) {
        //             if (!tuples[j].dominate(tuples[i])) {
        //                 if (k == 2) {
        //                     gamma[j] = 1;
        //                     for (int w = j + 1; w < i; ++ w) {
        //                         if (tuples[w].dominate(tuples[i])) {
        //                             gamma[j] = gamma[j] * (1 - tuples[w].prob);
        //                         }
        //                     }
        //                 }
        //                 if (Y[j][k - 1] != -1) {
        //                     BigFloat temp = tuples[i].prob;
        //                     temp = temp * gamma[j];
        //                     temp = temp * cost[j][k - 1];
        //                     if (temp > cost[i][k]) {
        //                         cost[i][k] = temp;
        //                         Y[i][k] = j;
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }
        Y[i][1] = i;
        cost[i][1] = BigFloat(tuples[i].prob);
        for (int j = 0; j < i; ++ j) {
            if (tuples[j].dominate(tuples[i])) {
                cost[i][1] = cost[i][1] * (1 - tuples[j].prob);
            }
        }
        for (int k = 2; k <= l; ++ k) {
            Y[i][k] = -1;
            cost[i][k] = 0;
        }
        for (int j = 0; j < i; ++ j) {
            if (!tuples[j].dominate(tuples[i])) {
                BigFloat gamma = 1;
                for (int k = 2; k <= l; ++ k) {
                    if (k == 2) {
                        for (int w = j + 1; w < i; ++ w) {
                            if (tuples[w].dominate(tuples[i])) {
                                gamma = gamma * (1 - tuples[w].prob);
                            }
                        }
                    }
                    if (Y[j][k - 1] != -1) {
                        BigFloat temp = tuples[i].prob;
                        temp = temp * gamma;
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
    // cout << "dp computation\n";
    int index = -1;
    for (int i = 0; i < n; ++ i) {
        if (cost[i][l] >= ucsky_prob) {
            ucsky_prob = cost[i][l];
            index = i;
        }
    }
    for (int i = l; i >= 1; -- i) {
        ucsky[i - 1] = index;
        index = Y[index][i];
    }
    for (int i = 0; i < n; ++ i) {
        delete[] Y[i];
        delete[] cost[i];
    }
    delete[] Y;
    delete[] cost;
    // delete[] gamma;
}

void UCSKY_Solver::group_dp() {
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        if (t.coord[0] < s.coord[0]) return true;
        else if (t.coord[0] == s.coord[0]) {
            if (t.coord[1] < s.coord[1]) return true;
        }
        return false;
    });
    vector<int> reduced_tuples;
    for (int i = 0; i < n; ++ i) {
        if (!tuples[i].reduced) reduced_tuples.push_back(i);
    }
    assert(reduced_tuples.size() == nR);
    int m = nR < n ? nR : n;
    // group tuples
    vector<vector<int>> groups;
    groups.push_back(vector<int>{0});
    vector<double> group2x;
    group2x.push_back(tuples[reduced_tuples[0]].coord[0]);
    for (int i = 1; i < m; ++ i) {
        if (tuples[reduced_tuples[i]].coord[0] != tuples[reduced_tuples[i - 1]].coord[0]) {
            groups.push_back(vector<int>{i});
            group2x.push_back(tuples[reduced_tuples[i]].coord[0]);
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
            A[q][i] = {tuples[reduced_tuples[groups[q][i]]].coord[1], -1, -1};
        }
    }

    BigFloat** gamma = new BigFloat*[m];
    for (int i = 0; i < m; ++ i) gamma[i] = new BigFloat[groups.size()];
    for (int q = 0; q < groups.size(); ++ q) {
        for (int i : groups[q]) {
            int t = reduced_tuples[i];
            for (int p = 0; p < q; ++ p) {
                gamma[i][p] = 1;
                Tuple virtual_t(-1, 2, nullptr, 0, false);
                virtual_t.coord[0] = group2x[p];
                virtual_t.coord[1] = tuples[t].coord[1];
                for (int k = 0; k < t; ++ k) {
                    if (tuples[k].dominate(tuples[t]) && !tuples[k].dominate(virtual_t)) {
                        gamma[i][p] = gamma[i][p] * (1 - tuples[k].prob);
                    }
                }
            }
        }
    }
    int** Y = new int*[m];
    BigFloat** cost = new BigFloat*[m];
    for (int i = 0; i < m; ++ i) {
        Y[i] = new int[l + 1];
        cost[i] = new BigFloat[l + 1];
    }
    for (int k = 1; k <= l; ++ k) {
        for (int q = 0; q < groups.size(); ++ q) {
            for (int u = 0; u < groups[q].size(); ++ u) {
                int i = groups[q][u];
                int t = reduced_tuples[i];
                if (k == 1) {
                    Y[i][k] = i;
                    cost[i][k] = BigFloat(tuples[t].prob);
                    for (int j = 0; j < t; ++ j) {
                        if (tuples[j].dominate(tuples[t])) {
                            cost[i][k] = cost[i][k] * (1 - tuples[j].prob);
                        }
                    }
                } else {
                    Y[i][k] = -1;
                    cost[i][k] = 0;
                    for (int p = 0; p < q; ++ p) {
                        int left = 0, right = A[p].size();
                        int mid = -1;
                        while (left < right) {
                            mid = (left + right)/2;
                            if (A[p][mid].y <= tuples[t].coord[1]) {
                                right = mid;
                            } else {
                                left = mid + 1;
                            }
                        }
                        -- left;
                        if (left >= 0 && A[p][left].y > tuples[t].coord[1]) {
                            int j = k%2 == 0 ? A[p][left].index2 : A[p][left].index1;
                            if (Y[j][k - 1] != -1) {
                                BigFloat temp = tuples[t].prob;
                                temp = temp * gamma[i][p];
                                temp = temp * cost[j][k - 1];
                                if (temp > cost[i][k]) {
                                    cost[i][k] = temp;
                                    Y[i][k] = j;
                                }
                            }
                        }
                    }
                }
                if (k%2 == 0) {
                    if (u == 0) { A[q][u].index1 = i; }
                    else {
                        if (cost[i][k] > cost[A[q][u - 1].index1][k]) {
                            A[q][u].index1 = i;
                        } else {
                            A[q][u].index1 = A[q][u - 1].index1;
                        }
                    }
                } else {
                    if (u == 0) { A[q][u].index2 = i; }
                    else {
                        if (cost[i][k] > cost[A[q][u - 1].index2][k]) {
                            A[q][u].index2 = i;
                        } else {
                            A[q][u].index2 = A[q][u - 1].index2;
                        }
                    }
                }
            }
        }
    }
    int index = -1;
    for (int i = 0; i < m; ++ i) {
        if (cost[i][l] >= ucsky_prob) {
            ucsky_prob = cost[i][l];
            index = i;
        }
    }
    for (int i = l; i >= 1; -- i) {
        ucsky[i - 1] = reduced_tuples[index];
        index = Y[index][i];
    }
    for (int i = 0; i < m; ++ i) {
        delete[] Y[i];
        delete[] cost[i];
        delete[] gamma[i];
    }
    delete[] Y;
    delete[] cost;
    delete[] gamma;
}

void UCSKY_Solver::group_dp2() {
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

    BigFloat** gamma = new BigFloat*[n];
    // for (int i = 0; i < n; ++ i) gamma[i] = new BigFloat[groups.size()];
    for (int q = 0; q < groups.size(); ++ q) {
        for (int i : groups[q]) gamma[i] = new BigFloat[q + 1];
    }
    for (int q = 0; q < groups.size(); ++ q) {
        for (int i : groups[q]) {
            for (int p = 0; p < q; ++ p) {
                gamma[i][p] = 1;
                Tuple virtual_t(-1, 2, nullptr, 0, false);
                virtual_t.coord[0] = group2x[p];
                virtual_t.coord[1] = tuples[i].coord[1];
                for (int k = 0; k < i; ++ k) {
                    if (tuples[k].dominate(tuples[i]) && !tuples[k].dominate(virtual_t)) {
                        gamma[i][p] = gamma[i][p] * (1 - tuples[k].prob);
                    }
                }
            }
        }
    }
    int** Y = new int*[n];
    BigFloat** cost = new BigFloat*[n];
    for (int i = 0; i < n; ++ i) {
        Y[i] = new int[l + 1];
        cost[i] = new BigFloat[l + 1];
    }
    for (int k = 1; k <= l; ++ k) {
        for (int q = 0; q < groups.size(); ++ q) {
            for (int u = 0; u < groups[q].size(); ++ u) {
                int i = groups[q][u];
                if (k == 1) {
                    Y[i][k] = i;
                    cost[i][k] = BigFloat(tuples[i].prob);
                    for (int j = 0; j < i; ++ j) {
                        if (tuples[j].dominate(tuples[i])) {
                            cost[i][k] = cost[i][k] * (1 - tuples[j].prob);
                        }
                    }
                } else {
                    Y[i][k] = -1;
                    cost[i][k] = 0;
                    for (int p = 0; p < q; ++ p) {
                        int left = 0, right = A[p].size();
                        int mid = -1;
                        while (left < right) {
                            mid = (left + right)/2;
                            if (A[p][mid].y <= tuples[i].coord[1]) {
                                right = mid;
                            } else {
                                left = mid + 1;
                            }
                        }
                        -- left;
                        if (left >= 0 && A[p][left].y > tuples[i].coord[1]) {
                            int j = k%2 == 0 ? A[p][left].index2 : A[p][left].index1;
                            if (Y[j][k - 1] != -1) {
                                BigFloat temp = tuples[i].prob;
                                temp = temp * gamma[i][p];
                                temp = temp * cost[j][k - 1];
                                if (temp > cost[i][k]) {
                                    cost[i][k] = temp;
                                    Y[i][k] = j;
                                }
                            }
                        }
                    }
                }
                if (k%2 == 0) {
                    if (u == 0) { A[q][u].index1 = i; }
                    else {
                        if (cost[i][k] > cost[A[q][u - 1].index1][k]) {
                            A[q][u].index1 = i;
                        } else {
                            A[q][u].index1 = A[q][u - 1].index1;
                        }
                    }
                } else {
                    if (u == 0) { A[q][u].index2 = i; }
                    else {
                        if (cost[i][k] > cost[A[q][u - 1].index2][k]) {
                            A[q][u].index2 = i;
                        } else {
                            A[q][u].index2 = A[q][u - 1].index2;
                        }
                    }
                }
            }
        }
    }

    int index = -1;
    for (int i = 0; i < n; ++ i) {
        if (cost[i][l] >= ucsky_prob) {
            ucsky_prob = cost[i][l];
            index = i;
        }
    }
    for (int i = l; i >= 1; -- i) {
        ucsky[i - 1] = index;
        index = Y[index][i];
    }
    for (int i = 0; i < n; ++ i) {
        delete[] Y[i];
        delete[] cost[i];
        delete[] gamma[i];
    }
    delete[] Y;
    delete[] cost;
    delete[] gamma;
}

void UCSKY_Solver::group_dp3() {
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
            reverse(groups.back().begin(), groups.back().end());
            groups.push_back(vector<int>{i});
            group2x.push_back(tuples[i].coord[0]);
        } else {
            groups.back().push_back(i);
        }
    }
    reverse(groups.back().begin(), groups.back().end());
    struct aux{
        double y;
        vector<int> index;
    };
    vector<vector<aux>> A;
    A.resize(groups.size());
    int** Y = new int*[n];
    BigFloat** cost = new BigFloat*[n];
    for (int i = 0; i < n; ++ i) {
        Y[i] = new int[l + 1];
        cost[i] = new BigFloat[l + 1];
    }
    int index = -1;
    for (int q = 0; q < groups.size(); ++ q) {
        for (int u = 0; u < groups[q].size(); ++ u) {
            int i = groups[q][u];
            Y[i][1] = i;
            cost[i][1] = BigFloat(tuples[i].prob);
            for (int j = 0; j < i; ++ j) {
                if (tuples[j].dominate(tuples[i])) {
                    cost[i][1] = cost[i][1] * (1 - tuples[j].prob);
                }
            }
            if (cost[i][1] > ucsky_prob) {
                for (int k = 2; k <= l; ++ k) {
                    Y[i][k] = -1;
                    cost[i][k] = 0;
                }
                for (int p = 0; p < q; ++ p) {
                    BigFloat gamma = 1;
                    Tuple virtual_t(-1, 2, nullptr, 0, false);
                    virtual_t.coord[0] = group2x[p];
                    virtual_t.coord[1] = tuples[i].coord[1];
                    for (int j = 0; j < i; ++ j) {
                        if (tuples[j].dominate(tuples[i]) && !tuples[j].dominate(virtual_t)) {
                            gamma = gamma * (1 - tuples[j].prob);
                        }
                    }
                    int left = 0, right = A[p].size();
                    int mid = -1;
                    while (left < right) {
                        mid = (left + right)/2;
                        if (A[p][mid].y <= tuples[i].coord[1]) {
                            right = mid;
                        } else {
                            left = mid + 1;
                        }
                    }
                    -- left;
                    if (left >= 0 && A[p][left].y > tuples[i].coord[1]) {
                        for (int k = 2; k <= l; ++ k) {
                            int j = A[p][left].index[k - 1];
                            if (Y[j][k - 1] != -1) {
                                BigFloat temp = tuples[i].prob;
                                temp = temp * gamma;
                                temp = temp * cost[j][k - 1];
                                if (temp > cost[i][k]) {
                                    cost[i][k] = temp;
                                    Y[i][k] = j;
                                }
                            }
                        }
                    }
                }
                A[q].push_back({tuples[i].coord[1], vector<int>{}});
                A[q].back().index.resize(l + 1);
                if (A[q].size() == 1) {
                    for (int k = 1; k <= l; ++ k) {
                        A[q].back().index[k] = i;
                    }
                } else {
                    for(int k = 1; k <= l; ++ k) {
                        int s = A[q].size() - 2;
                        if (cost[i][k] > cost[A[q][s].index[k]][k]) {
                            A[q].back().index[k] = i;
                        } else {
                            A[q].back().index[k] = A[q][s].index[k];
                        }
                    }
                }
                if (cost[i][l] >= ucsky_prob) {
                    ucsky_prob = cost[i][l];
                    index = i;
                }
            }
        }
        // for (auto a : A[q]) {
            // cout << a.y << '\t';
            // int i = a.index[2];
            // cout << i << '\t' << tuples[i].coord[1] << '\t' << cost[i][2] << endl;
        // }
    }
    for (int i = l; i >= 1; -- i) {
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

void UCSKY_Solver::branch_bound() {
    /* vector<int> reduced_tuples;
    for (int i = 0; i < n; ++ i) {
        if (!tuples[i].reduced) reduced_tuples.push_back(i);
    }
    assert(reduced_tuples.size() == nR);
    int m = nR < n ? nR : n;
    vector<int> S;
    BigFloat prob(0);
    search(reduced_tuples, S, prob, l); */
    vector<pair<int, double>> reduced_tuples;
    for (int i = 0; i < n; ++ i) {
        if (!tuples[i].reduced) {
            reduced_tuples.push_back(make_pair(i, tuples[i].sum()));
        }
    }
    assert(reduced_tuples.size() == nR);
    int m = nR < n ? nR : n;
    divide(reduced_tuples);
}

void UCSKY_Solver::divide(vector<pair<int, double>>& T) {
    vector<pair<int, double>> sky;
    skyline(T, sky);
    sort(sky.begin(), sky.end(), [](const pair<int, double>& t, const pair<int, double>& s){
        return t.second > s.second;
    });
    BigFloat prob = 1;
    for (int i = 1; i <= min(l, (int)sky.size()); ++ i) {
        vector<int> subset(i);
        prob = prob * sky[i - 1].second;
        if (prob > ucsky_prob) {
            BigFloat prob2 = 1;
            enum_subset(T, sky, i, 0, subset, 0, prob2);
        } else {
            break;
        }
    }
    if (T.size() >= l) divide(T);
}

void UCSKY_Solver::skyline(vector<pair<int, double>>& T, vector<pair<int, double>>& sky) {
    sort(T.begin(), T.end(), [](const pair<int, double>& t, const pair<int, double>& s){
        return t.second < s.second;
    });
    vector<pair<int, double>> new_T;
    for (int i = 0; i < T.size(); ++ i) {
        bool flag = false;
        for (int j = 0; j < sky.size(); ++ j) {
            if (tuples[sky[j].first].dominate(tuples[T[i].first])) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            sky.push_back(make_pair(T[i].first, tuples[T[i].first].prob));
        } else {
            new_T.push_back(T[i]);
        }
    }
    T = new_T;
}

void UCSKY_Solver::enum_subset(const vector<pair<int, double>>& T, const vector<pair<int, double>>& sky, int s, int j, vector<int>& subset, int i, BigFloat& prob) {
    if (i == s) {
        vector<int> new_T;
        for (auto t : T) {
            bool flag = false;
            for (auto s : subset) {
                if (tuples[s].dominate(tuples[t.first])) {
                    flag = true;
                    break;
                }
            }
            if (!flag) new_T.push_back(t.first);
        }
        BigFloat prob = compute_prob(subset);
        if (new_T.size() >= l - s && prob > ucsky_prob) {
            search(new_T, subset, prob, l - s);
        }
        return;
    }
    if (j >= sky.size()) { return; }
    subset[i] = sky[j].first;
    prob = prob * tuples[subset[i]].prob;
    if (prob > ucsky_prob) {
        enum_subset(T, sky, s, j + 1, subset, i + 1, prob);
    }
    prob = prob / tuples[subset[i]].prob;
    if (prob > ucsky_prob) {
        enum_subset(T, sky, s, j + 1, subset, i, prob);
    }
}

void UCSKY_Solver::search(const vector<int>& T, const vector<int>& S, const BigFloat& prob, const int& k) {
    if (k == 0) {
        if (prob > ucsky_prob) {
            ucsky_prob = prob;
            ucsky = S;
        }
        return;
    }
    if (T.size() < k) { return; }
    if (prob < ucsky_prob) { return;}

    vector<double> prob_T(T.size());
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
    if (beta < ucsky_prob) { return; }

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
        search(new_T, new_S, new_S_prob, k - 1);
    }
}

void UCSKY_Solver::sample_FIM(const int& theta) {
    vector<int> reduced_tuples;
    for (int i = 0; i < n; ++ i) {
        if (!tuples[i].reduced) reduced_tuples.push_back(i);
    }
    assert(reduced_tuples.size() == nR);
    int m = nR < n ? nR : n;
    bool flag = false;
    vector<vector<int>> S(theta, vector<int>{});
    int min_sup = 0;
    vector<int> subset(l);
    vector<int> degree(m);
    for (int i = 0; i < theta; ++ i) {
        vector<int> C;
        for (int j = 0; j < m; ++ j) {
            degree[j] = inN[j].size();
            if (degree[j] == 0) C.push_back(j);
        }
        while (!C.empty()) {
            int t = reduced_tuples[C.back()];
            C.pop_back();
            double p = rand_uniform(0, 1);
            if (p <= tuples[t].prob) {
                if (!tuples[t].reduced) S[i].push_back(t);
            } else {
                for (auto s : outN[t]) {
                    -- degree[s];
                    if (degree[s] == 0) {
                        C.push_back(s);
                    }
                }
            }
        }
        if (S[i].size() >= l) {
            flag = true;
            if (contain(S[i], ucsky)) {
                ++ min_sup;
            }
        }
    }
    if (!flag) {
        cout << "no pattern with length " << l << endl;
        return;
    }
    min_sup = max(min_sup, 1);
    FPTree fptree(S, min_sup);
    // set<Pattern> patterns = fptree_growth(fptree, l);
    set<Pattern> patterns = fptree_growth(fptree);
    int max_freq = 0;
    unordered_map<int, vector<set<int>>> l_set;
    // cout << patterns.size() << endl;
    for (auto p : patterns) {
        if (p.first.size() == l) {
            max_freq = max(max_freq, (int)p.second);
            l_set[p.second].push_back(p.first);
        }
    }
    ucsky_prob = 0;
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
}

void UCSKY_Solver::sample_FIM2(const int& theta) {
    sort(tuples.begin(), tuples.end(), [](const Tuple& t, const Tuple& s){
        // return t.sum() < s.sum();
        if (t.min() < s.min()) return true;
        else if (t.min() == s.min()) return t.sum() < s.sum();
        else return false;
    });
    vector<vector<int>> S(theta, vector<int>{});
    double min_sky = __DBL_MAX__;
    bool flag = false;
    for (int i = 0; i < theta; ++ i) {
        for (int j = 0; j < n; ++ j) {
            double p = rand_uniform(0, 1);
            if (p <= tuples[j].prob) {
                if (tuples[j].min() > min_sky) {
                    break;
                } else {
                    bool dominate = false;
                    for (int k : S[i]) {
                        if (tuples[k].dominate(tuples[j])) {
                            dominate = true;
                            break;
                        }
                    }
                    if (!dominate) {
                        S[i].push_back(j);
                        min_sky = min(min_sky, tuples[j].max());
                    }
                }
            }
        }
        if (S[i].size() >= l) {
            flag = true;
        }
        // cout << S[i].size() << endl;
    }
    if (!flag) {
        cout << "no pattern with length " << l << endl;
        // return;
        exit(0);
    }
    // min_sup = max(min_sup, 1);
    FPTree fptree(S, 1);
    // set<Pattern> patterns = fptree_growth(fptree, l);
    set<Pattern> patterns = fptree_growth(fptree);
    int max_freq = 0;
    unordered_map<int, vector<set<int>>> l_set;
    // cout << patterns.size() << endl;
    for (auto p : patterns) {
        if (p.first.size() == l) {
            max_freq = max(max_freq, (int)p.second);
            l_set[p.second].push_back(p.first);
        }
    }
    ucsky_prob = 0;
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
}


void UCSKY_Solver::check_ucsky() {
    for (int i = 0; i < l; ++ i) {
        for (int j = 0; j < i; ++ j) {
            if (tuples[ucsky[i]].dominate(tuples[ucsky[j]]) || tuples[ucsky[j]].dominate(tuples[ucsky[i]])) {
                cout << "conflict (" << tuples[ucsky[i]].id << ", " << tuples[ucsky[j]].id << ")\n";
            }
        }
    }
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
    cout << "real: " << prob << "\t prob: " << ucsky_prob << endl;
}

void UCSKY_Solver::print_ucsky() {
    cout << "UCSKY: ";
    for (int i : ucsky) {
        cout << tuples[i] << ", ";
        // cout << tuples[i] << endl;
    }
    cout << "\t Prob: " << ucsky_prob << endl;
}

BigFloat UCSKY_Solver::compute_prob(const vector<int>& S) {
    BigFloat prob(1);
    for (int i = 0; i < n; ++ i) {
        for (int j : S) {
            if (tuples[i].dominate(tuples[j])) {
                prob = prob * (1 - tuples[i].prob);
                break;
            }
        }
    }
    for (int j : S) {
        prob = prob * tuples[j].prob;
    }
    return prob;
}

BigFloat UCSKY_Solver::compute_prob(const set<int>& S) {
    BigFloat prob(1);
    for (int i = 0; i < n; ++ i) {
        for (int j : S) {
            if (tuples[i].dominate(tuples[j])) {
                prob = prob * (1 - tuples[i].prob);
                break;
            }
        }
    }
    for (int j : S) {
        prob = prob * tuples[j].prob;
    }
    return prob;
}

bool UCSKY_Solver::contain(const vector<int>& T, const vector<int>& subset) {
    for (auto s : subset) {
        bool found = false;
        for (auto t : T) {
            if (s == t) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}