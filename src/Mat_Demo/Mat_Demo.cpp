//
// Created by Chiahwa Young on 2023/3/17.
//

#include "Mat_Demo.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include "Timer/Timer.h"

//default constructor.
HwaUtil::Mat_Demo::Mat_Demo() {
    nrows = 0;
    ncols = 0;
    d = nullptr;
}

//constructor with initialization.
HwaUtil::Mat_Demo::Mat_Demo(const int nr, const int nc, const HwaUtil::Mat_Demo::MatrixType initType) {
    Timer::tick("HwaUtil::Mat_Demo","()");
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
    Timer::tock("HwaUtil::Mat_Demo","()");
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
    Timer::tick("HwaUtil::Mat_Demo","mmax");
    double max = d[0];
    for (int i = 0; i < nrows * ncols; i++) {
        if (d[i] > max) {
            max = d[i];
        }
    }
    Timer::tock("HwaUtil::Mat_Demo","mmax");
    return max;
}

double HwaUtil::Mat_Demo::mmin() const {
    Timer::tick("HwaUtil::Mat_Demo","mmin");
    double min = d[0];
    for (int i = 0; i < nrows * ncols; i++) {
        if (d[i] < min) {
            min = d[i];
        }
    }
    Timer::tock("HwaUtil::Mat_Demo","mmin");
    return min;
}

void HwaUtil::Mat_Demo::zero() {
    Timer::tick("HwaUtil::Mat_Demo","zero");
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] = 0;
    }
    Timer::tock("HwaUtil::Mat_Demo","zero");
}

void HwaUtil::Mat_Demo::zero(const int nr, const int nc) {
    Timer::tick("HwaUtil::Mat_Demo","zero");
    nrows = nr;
    ncols = nc;
    delete[] d;
    d = new double[nr * nc];
    zero();
    Timer::tock("HwaUtil::Mat_Demo","zero");
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator*=(const double &a) {
    Timer::tick("HwaUtil::Mat_Demo","operator*=");
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] *= a;
    }
    Timer::tock("HwaUtil::Mat_Demo","operator*=");
    return *this;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator+=(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo","operator+=");
    if (nrows != a.nrows || ncols != a.ncols) {
        Timer::tock("HwaUtil::Mat_Demo","operator+=");
        throw std::runtime_error("Error: Mat_Demo::operator+=: matrix size mismatch");
    }
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] += a.d[i];
    }
    Timer::tock("HwaUtil::Mat_Demo","operator+=");
    return *this;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator-=(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo","operator-=");
    if (nrows != a.nrows || ncols != a.ncols) {
        Timer::tock("HwaUtil::Mat_Demo","operator-=");
        throw std::runtime_error("Error: Mat_Demo::operator-=: matrix size mismatch");
    }
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] -= a.d[i];
    }
    Timer::tock("HwaUtil::Mat_Demo","operator-=");
    return *this;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator+(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo","operator+");
    if (nrows != a.nrows || ncols != a.ncols) {
        Timer::tock("HwaUtil::Mat_Demo","operator+");
        throw std::runtime_error("Error: Mat_Demo::operator+: matrix size mismatch");
    }
    auto *m = new HwaUtil::Mat_Demo(nrows, ncols, MatrixType::Zero);
    for (int i = 0; i < nrows * ncols; i++) {
        m->d[i] = d[i] + a.d[i];
    }
    Timer::tock("HwaUtil::Mat_Demo","operator+");
    return *m;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator-(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo","operator-");
    if (nrows != a.nrows || ncols != a.ncols) {
        Timer::tock("HwaUtil::Mat_Demo","operator-");
        throw std::runtime_error("Error: Mat_Demo::operator-: matrix size mismatch");
    }
    auto m = new HwaUtil::Mat_Demo(nrows, ncols, MatrixType::Zero);
    for (int i = 0; i < nrows * ncols; i++) {
        m->d[i] = d[i] - a.d[i];
    }
    Timer::tock("HwaUtil::Mat_Demo","operator-");
    return *m;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator=(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo","operator=");
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
    Timer::tock("HwaUtil::Mat_Demo","operator=");
    return *this;
}

double &HwaUtil::Mat_Demo::operator()(const int i, const int j) {
    if (i < 0 || i >= nrows || j < 0 || j >= ncols) {
        throw std::runtime_error("Error: Mat_Demo::operator(): index out of range");
    }
    return d[i * ncols + j];
}

const double &HwaUtil::Mat_Demo::operator()(const int i, const int j) const {
    if (i < 0 || i >= nrows || j < 0 || j >= ncols) {
        throw std::runtime_error("Error: Mat_Demo::operator(): index out of range");
    }
    return d[i * ncols + j];
}

std::ostream &HwaUtil::operator<<(std::ostream &os, const HwaUtil::Mat_Demo &m) {
    Timer::tick("HwaUtil::(root)","operator<<(Mat_Demo)");
    for (int i = 0; i < m.nrows; i++) {
        for (int j = 0; j < m.ncols; j++) {
            os << m(i, j) << " ";
        }
        os << std::endl;
    }
    Timer::tock("HwaUtil::(root)","operator<<(Mat_Demo)");
    return os;
}

