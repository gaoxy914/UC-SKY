#include "Tuple.h"

Tuple::Tuple() {
    id = -1;
    dim = 0;
    coord = nullptr;
    prob = 0;
    reduced = false;
}

Tuple::Tuple(const int& id, const int& dim, const double* coord, const double& prob, const bool& reduced) {
    this->id = id;
    this->dim = dim;
    this->prob = prob;
    this->reduced = reduced;
    this->coord = new double[this->dim];
    if (coord != nullptr) {
        memcpy(this->coord, coord, this->dim*sizeof(double));
    } else {
        memset(this->coord, 0, this->dim*sizeof(double));
    }
}

Tuple::Tuple(const Tuple& other) {
    id = other.id;
    dim = other.dim;
    prob = other.prob;
    reduced = other.reduced;
    coord = new double[dim];
    memcpy(coord, other.coord, dim*sizeof(double));
}

Tuple::~Tuple() {
    if (coord != nullptr) {
        delete[] coord;
        coord = nullptr;
    }
    dim = 0;
    prob = 0;
    reduced = false;
}

Tuple& Tuple::operator= (const Tuple& other) {
    if (&other != this) {
        id = other.id;
        prob = other.prob;
        reduced = other.reduced;
        if (dim != other.dim) {
            dim = other.dim;
            if (coord != nullptr) {
                delete[] coord;
            }
            coord = new double[dim];
        }
        memcpy(coord, other.coord, dim*sizeof(double));
    }
    return *this;
}

double Tuple::operator[] (const int& index) {
    return coord[index];
}

bool Tuple::dominate(const Tuple& other) {
    bool unique = false;
    for (int i = 0; i < dim; ++ i) {
        if (coord[i] > other.coord[i]) {
            return false;
        } else if (coord[i] < other.coord[i]) {
            unique = true;
        }
    }
    return unique;
}

bool Tuple::dominate(const Tuple& other) const {
    bool unique = false;
    for (int i = 0; i < dim; ++ i) {
        if (coord[i] > other.coord[i]) {
            return false;
        } else if (coord[i] < other.coord[i]) {
            unique = true;
        }
    }
    return unique;
}

double Tuple::sum() const {
    double res = coord[0];
    for (int i = 1; i < dim; ++ i) {
        res += coord[i];
    }
    return res;
}

double Tuple::min() const {
    double res = coord[0];
    for (int i = 1; i < dim; ++ i) {
        res = std::min(res, coord[i]);
    }
    return res;
}

double Tuple::max() const {
    double res = coord[0];
    for (int i = 1; i < dim; ++ i) {
        res = std::max(res, coord[i]);
    }
    return res;
}

ostream& operator<< (ostream& out, const Tuple& t) {
    out << "[" << t.id;
    out << " (";
    for (int i = 0; i < t.dim - 1; ++ i) {
        out << t.coord[i] << ", ";
    }
    out << t.coord[t.dim - 1] << "), " << t.prob << "]";
    return out;
}