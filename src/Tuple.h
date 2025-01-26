#ifndef __TUPLE_H__
#define __TUPLE_H__

#include "BigFloat.h"

class Tuple {
public:
    int id;
    int dim;
    double *coord;
    double prob;
    bool reduced;

    Tuple();
    Tuple(const int& id, const int& dim, const double* coord, const double& prob, const bool& reduced);
    Tuple(const Tuple& other);
    virtual ~Tuple();
    Tuple& operator= (const Tuple& other);
    double operator[] (const int& index);
    bool dominate(const Tuple& other);
    bool dominate(const Tuple& other) const;
    double sum() const;
    double min() const;
    double max() const;
    friend ostream& operator<< (ostream& out, const Tuple& t);
};

#endif