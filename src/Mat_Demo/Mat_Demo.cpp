//
// Created by Chiahwa Young on 2023/3/17.
//

#include "Mat_Demo.h"
#include <cmath>
#include <iostream>
#include <algorithm>

HwaUtil::Mat_Demo::Mat_Demo() {
    nrows = 0;
    ncols = 0;
    d = nullptr;
}

HwaUtil::Mat_Demo::Mat_Demo(const int nr, const int nc, const HwaUtil::Mat_Demo::MatrixType initType) {
    nrows = nr;
    ncols = nc;
    d = new double[nr * nc];
    switch (initType) {
        case MatrixType::Zero:
            zero();
            break;
        case MatrixType::Identity:
            zero();
            for (int i = 0; i < std::min(ncols,nrows); i++) {
                d[i * ncols + i] = 1;
            }
            break;
        case MatrixType::Random:
            for (int i = 0; i < nrows * ncols; i++) {
                d[i] = rand() / (double) RAND_MAX;
            }
            break;
        case MatrixType::User:
            break;
    }
}

HwaUtil::Mat_Demo::~Mat_Demo() {
    delete[] d;
}

int HwaUtil::Mat_Demo::nr() const {
    return nrows;
}

int HwaUtil::Mat_Demo::nc() const {
    return ncols;
}

double HwaUtil::Mat_Demo::mmax() const {
    double max = d[0];
    for (int i = 0; i < nrows * ncols; i++) {
        if (d[i] > max) {
            max = d[i];
        }
    }
    return max;
}

double HwaUtil::Mat_Demo::mmin() const {
    double min = d[0];
    for (int i = 0; i < nrows * ncols; i++) {
        if (d[i] < min) {
            min = d[i];
        }
    }
    return min;
}

void HwaUtil::Mat_Demo::zero() {
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] = 0;
    }
}

void HwaUtil::Mat_Demo::zero(const int nr, const int nc) {
    nrows = nr;
    ncols = nc;
    delete[] d;
    d = new double[nr * nc];
    zero();
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator*=(const double &a) {
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] *= a;
    }
    return *this;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator+=(const HwaUtil::Mat_Demo &a) {
    if (nrows != a.nrows || ncols != a.ncols) {
        std::cout << "Error: Mat_Demo::operator+=: matrix size mismatch" << std::endl;
        exit(1);
    }
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] += a.d[i];
    }
    return *this;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator-=(const HwaUtil::Mat_Demo &a) {
    if (nrows != a.nrows || ncols != a.ncols) {
        std::cout << "Error: Mat_Demo::operator-=: matrix size mismatch" << std::endl;
        exit(1);
    }
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] -= a.d[i];
    }
    return *this;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator+(const HwaUtil::Mat_Demo &a) {
    if (nrows != a.nrows || ncols != a.ncols) {
        std::cout << "Error: Mat_Demo::operator+: matrix size mismatch" << std::endl;
        exit(1);
    }
    HwaUtil::Mat_Demo *m = new HwaUtil::Mat_Demo(nrows, ncols, MatrixType::Zero);
    for (int i = 0; i < nrows * ncols; i++) {
        m->d[i] = d[i] + a.d[i];
    }
    return *m;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator-(const HwaUtil::Mat_Demo &a) {
    if (nrows != a.nrows || ncols != a.ncols) {
        std::cout << "Error: Mat_Demo::operator-: matrix size mismatch" << std::endl;
        exit(1);
    }
    HwaUtil::Mat_Demo *m = new HwaUtil::Mat_Demo(nrows, ncols, MatrixType::Zero);
    for (int i = 0; i < nrows * ncols; i++) {
        m->d[i] = d[i] - a.d[i];
    }
    return *m;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator=(const HwaUtil::Mat_Demo &a) {
    if (this == &a) {
        return *this;
    }
    nrows = a.nrows;
    ncols = a.ncols;
    delete[] d;
    d = new double[nrows * ncols];
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] = a.d[i];
    }
    return *this;
}

double &HwaUtil::Mat_Demo::operator()(const int i, const int j) {
    if (i < 0 || i >= nrows || j < 0 || j >= ncols) {
        std::cout << "Error: Mat_Demo::operator(): index out of range" << std::endl;
        exit(1);
    }
    return d[i * ncols + j];
}

const double &HwaUtil::Mat_Demo::operator()(const int i, const int j) const {
    if (i < 0 || i >= nrows || j < 0 || j >= ncols) {
        std::cout << "Error: Mat_Demo::operator(): index out of range" << std::endl;
        exit(1);
    }
    return d[i * ncols + j];
}

std::ostream &HwaUtil::operator<<(std::ostream &os, const HwaUtil::Mat_Demo &m) {
    for (int i = 0; i < m.nrows; i++) {
        for (int j = 0; j < m.ncols; j++) {
            os << m(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}

