#include "Dataset.h"

DG::DG(const int& n) {

}

DG::~DG() {
    
}

Dataset::Dataset(const int& dim, const int& n, const double& center) {
    srand((unsigned)time(nullptr));
    srand48((unsigned)time(nullptr));
    this->dim = dim;
    this->n = n;
    this->center = center;
    this->tuples.resize(n, Tuple(-1, dim, nullptr, 0));
    stringstream stream;
    stream << fixed << setprecision(2) << center;
    data_path = to_string(dim) + "_" + stream.str();
}

Dataset::~Dataset() {

}

void Dataset::gen_ind_data() {
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

void Dataset::gen_anti_data() {
    double delta = min(center, 1 - center);
    double lower = center - delta, upper = center + delta;
    double x[dim];
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

void Dataset::gen_corr_data() {
    double delta = min(center, 1 - center);
    double lower = center - delta, upper = center + delta;
    double x[dim];
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

void Dataset::load_data(const char* data_path) {
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

void Dataset::write_data(const char* data_path) {
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

void Dataset::print_data() {
    for (int i = 0; i < tuples.size(); ++ i) {
        cout << i << "," << tuples[i] << endl;
    }
}

DG Dataset::construct_graph() {

}