//
// Created by Chiahwa Young on 2023/3/17.
//

#include "Mat_Demo.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cstring>

#if defined(__APPLE__) || defined(__MACOSX)

#include <Accelerate/Accelerate.h>

#elif defined(__linux__)
#include <cblas.h>
#include <lapacke.h>
#endif

#include "Timer/Timer.h"

//default constructor.
HwaUtil::Mat_Demo::Mat_Demo() {
    nrows = 0;
    ncols = 0;
    d = nullptr;
}

//constructor with initialization.
HwaUtil::Mat_Demo::Mat_Demo(int nr, int nc, MatrixType initType) {
    Timer::tick("HwaUtil::Mat_Demo", "()");
    nrows = nr;
    ncols = nc;
    d = new double[nr * nc];
    switch (initType) {
        case MatrixType::Zero:
            zero();
            break;
        case MatrixType::Identity:
            zero();
            for (int i = 0; i < std::min(ncols, nrows); i++) {
                d[i * ncols + i] = 1;
            }
            break;
        case MatrixType::Random:
            for (int i = 0; i < nrows * ncols; i++) {
                d[i] = rand() / (double) RAND_MAX;
            }
            break;

    }
    Timer::tock("HwaUtil::Mat_Demo", "()");
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
    Timer::tick("HwaUtil::Mat_Demo", "mmax");
    double max = d[0];
    for (int i = 0; i < nrows * ncols; i++) {
        if (d[i] > max) {
            max = d[i];
        }
    }
    Timer::tock("HwaUtil::Mat_Demo", "mmax");
    return max;
}

double HwaUtil::Mat_Demo::mmin() const {
    Timer::tick("HwaUtil::Mat_Demo", "mmin");
    double min = d[0];
    for (int i = 0; i < nrows * ncols; i++) {
        if (d[i] < min) {
            min = d[i];
        }
    }
    Timer::tock("HwaUtil::Mat_Demo", "mmin");
    return min;
}

void HwaUtil::Mat_Demo::zero() {
    Timer::tick("HwaUtil::Mat_Demo", "zero");
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] = 0;
    }
    Timer::tock("HwaUtil::Mat_Demo", "zero");
}

void HwaUtil::Mat_Demo::zero(const int nr, const int nc) {
    Timer::tick("HwaUtil::Mat_Demo", "zero");
    nrows = nr;
    ncols = nc;
    delete[] d;
    d = new double[nr * nc];
    zero();
    Timer::tock("HwaUtil::Mat_Demo", "zero");
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator*=(const double &a) {
    Timer::tick("HwaUtil::Mat_Demo", "operator*=");
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] *= a;
    }
    Timer::tock("HwaUtil::Mat_Demo", "operator*=");
    return *this;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator+=(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo", "operator+=");
    if (nrows != a.nrows || ncols != a.ncols) {
        Timer::tock("HwaUtil::Mat_Demo", "operator+=");
        throw std::runtime_error("Error: Mat_Demo::operator+=: matrix size mismatch");
    }
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] += a.d[i];
    }
    Timer::tock("HwaUtil::Mat_Demo", "operator+=");
    return *this;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator-=(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo", "operator-=");
    if (nrows != a.nrows || ncols != a.ncols) {
        Timer::tock("HwaUtil::Mat_Demo", "operator-=");
        throw std::runtime_error("Error: Mat_Demo::operator-=: matrix size mismatch");
    }
    for (int i = 0; i < nrows * ncols; i++) {
        d[i] -= a.d[i];
    }
    Timer::tock("HwaUtil::Mat_Demo", "operator-=");
    return *this;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator+(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo", "operator+");
    if (nrows != a.nrows || ncols != a.ncols) {
        Timer::tock("HwaUtil::Mat_Demo", "operator+");
        throw std::runtime_error("Error: Mat_Demo::operator+: matrix size mismatch");
    }
    auto *m = new HwaUtil::Mat_Demo(nrows, ncols, MatrixType::Zero);
    for (int i = 0; i < nrows * ncols; i++) {
        m->d[i] = d[i] + a.d[i];
    }
    Timer::tock("HwaUtil::Mat_Demo", "operator+");
    return *m;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator-(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo", "operator-");
    if (nrows != a.nrows || ncols != a.ncols) {
        Timer::tock("HwaUtil::Mat_Demo", "operator-");
        throw std::runtime_error("Error: Mat_Demo::operator-: matrix size mismatch");
    }
    auto m = new HwaUtil::Mat_Demo(nrows, ncols, MatrixType::Zero);
    for (int i = 0; i < nrows * ncols; i++) {
        m->d[i] = d[i] - a.d[i];
    }
    Timer::tock("HwaUtil::Mat_Demo", "operator-");
    return *m;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator=(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo", "operator=");
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
    Timer::tock("HwaUtil::Mat_Demo", "operator=");
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
    Timer::tick("HwaUtil::(root)", "operator<<(Mat_Demo)");
    os << "# rows of the given matrix" << std::endl
       << "nrows: " << m.nrows << std::endl;
    os << "# columns of the given matrix" << std::endl
       << "ncols: " << m.ncols << std::endl;
    os << "# data type of the given matrix" << std::endl
       << "type: " << "double" << std::endl; //TODO:support other types
    os << "# element value" << std::endl
       << "value: " << std::endl;

    for (int i = 0; i < m.nrows; i++) {
        for (int j = 0; j < m.ncols; j++) {
            os << m(i, j) << (j == m.ncols - 1 ? "" : ", ");
        }
        os << std::endl;
    }
    Timer::tock("HwaUtil::(root)", "operator<<(Mat_Demo)");
    return os;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::operator*(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo", "operator*");
    if (ncols != a.nrows) {
        Timer::tock("HwaUtil::Mat_Demo", "operator*");
        throw std::runtime_error("Error: Mat_Demo::operator*: matrix size mismatch");
    }
    auto *m = new HwaUtil::Mat_Demo(nrows, a.ncols, MatrixType::Zero);
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < a.ncols; j++) {
            for (int k = 0; k < ncols; k++) {
                (*m)(i, j) += (*this)(i, k) * a(k, j);
            }
        }
    }
    Timer::tock("HwaUtil::Mat_Demo", "operator*");
    return *m;
}

int HwaUtil::Mat_Demo::is_real_symm() const {
    Timer::tick("HwaUtil::Mat_Demo", "is_real_symm");
    if (nrows != ncols) {
        Timer::tock("HwaUtil::Mat_Demo", "is_real_symm");
        return 0;
    }
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < i; j++) {
            if (std::abs((*this)(i, j) - (*this)(j, i)) > 1e-10) {
                Timer::tock("HwaUtil::Mat_Demo", "is_real_symm");
                return 0;
            }
        }
    }
    Timer::tock("HwaUtil::Mat_Demo", "is_real_symm");
    return 1;
}

HwaUtil::Mat_Demo &HwaUtil::Mat_Demo::blas_mult(const HwaUtil::Mat_Demo &a) {
    Timer::tick("HwaUtil::Mat_Demo", "blas_mult");
    if (ncols != a.nrows) {
        Timer::tock("HwaUtil::Mat_Demo", "blas_mult");
        throw std::runtime_error("Error: Mat_Demo::blas_mult: matrix size mismatch");
    }
    auto *m = new HwaUtil::Mat_Demo(nrows, a.ncols, MatrixType::Zero);
    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                nrows, a.ncols, ncols, 1.0, d, ncols, a.d, a.ncols,
                0.0, m->d, a.ncols);
    Timer::tock("HwaUtil::Mat_Demo", "blas_mult");
    return *m;
}

int HwaUtil::Mat_Demo::lapack_eig(double *eigval, double *eigvec) {
    Timer::tick("HwaUtil::Mat_Demo", "lapack_eig");
    if (!is_real_symm()) {
        Timer::tock("HwaUtil::Mat_Demo", "lapack_eig");
        std::cerr << "Error: Mat_Demo::lapack_eig: matrix is not real symmetric" << std::endl;
        return 0;
    }
    int info = 0;
    int lwork = 3 * nrows - 1;
    double *work = new double[lwork];
#if defined(__APPLE__) || defined(__MACOSX)
    dsyev_("V", "U", &nrows, d, &nrows, eigval, work, &lwork, &info);
#else if defined(__linux__)
    info=LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nrows, d, nrows, eigval);
#endif
    if (info != 0) {
        Timer::tock("HwaUtil::Mat_Demo", "lapack_eig");
        return 0;
    }
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < nrows; j++) {
            eigvec[i * nrows + j] = d[i * nrows + j];
        }
    }
    Timer::tock("HwaUtil::Mat_Demo", "lapack_eig");
    return 1;
}

HwaUtil::Mat_Demo::Mat_Demo(std::istream &is) {
    Timer::tick("HwaUtil::Mat_Demo", "(istream)");
    constexpr auto max_size = std::numeric_limits<std::streamsize>::max();
    std::string line;
    bool ncolget = false, nrowget = false, typeget = false, vallabelget = false;
    while (std::getline(is, line)) {
        std::stringstream buf(line);
        getline(buf, line, '#');
        std::stringstream ssline(line);

        std::string name;
        ssline >> name;
        if (ssline.fail()) {
            continue;
        }
        transform(name.begin(), name.end(), name.begin(), ::tolower);
        name = name.substr(0, name.find(':'));
        if (name == "value") {
            vallabelget = true;
            if (ncolget && nrowget && typeget) {
                break;
            } else {
                Timer::tock("HwaUtil::Mat_Demo", "(istream)");
                throw std::runtime_error("Missing argument!");
            }
        }
        std::string val;
        ssline >> val;
        if (ssline.fail()) {
            Timer::tock("HwaUtil::Mat_Demo", "(istream)");
            throw std::runtime_error("Value of argument " + name + " missing!");
        }
        transform(val.begin(), val.end(), val.begin(), ::tolower);
        if (name == "nrows") {
            nrows = std::stoi(val);
            nrowget = true;
        } else if (name == "ncols") {
            ncols = std::stoi(val);
            ncolget = true;
        } else if (name == "type") {
            //TODO: construct according to specified type
            typeget = true;
        } else {
            Timer::tock("HwaUtil::Mat_Demo", "(istream)");
            throw std::runtime_error("Unknown argument: " + name);
        }
    }
    if (!vallabelget) {
        Timer::tock("HwaUtil::Mat_Demo", "(istream)");
        throw std::runtime_error("Missing \"value\" label!");
    }
    d = new double[nrows * ncols]; //assert: no comment in the middle of the matrix

    for (int i = 0; i < nrows * ncols; i++) {
        while (!isdigit(is.peek()) && is.peek() != '-' && is.peek() != '.') {
            is.ignore();
        }
        is >> d[i];
        if (is.fail()) {
            Timer::tock("HwaUtil::Mat_Demo", "(istream)");
            throw std::runtime_error("Error reading matrix value!");
        }
    }
    Timer::tock("HwaUtil::Mat_Demo", "(istream)");
}

const double *HwaUtil::Mat_Demo::get_ptr(int i, int j) {
    return d + i * ncols + j;
}

HwaUtil::Mat_Demo::Mat_Demo(int nr, int nc, double *data) {
    Timer::tick("HwaUtil::Mat_Demo", "(nr,nc,data)");
    nrows = nr;
    ncols = nc;
    d = new double[nrows * ncols];
    if (data != nullptr) {
        memcpy(d, data, nrows * ncols * sizeof(double));
    }
    Timer::tock("HwaUtil::Mat_Demo", "(nr,nc,data)");
}

HwaUtil::Mat_Demo &HwaUtil::mat_add(const HwaUtil::Mat_Demo &a, const HwaUtil::Mat_Demo &b, double alpha, double beta) {
    Timer::tick("HwaUtil::Mat_Demo", "mat_add");
    if (a.nrows != b.nrows || a.ncols != b.ncols) {
        Timer::tock("HwaUtil::Mat_Demo", "mat_add");
        throw std::runtime_error("Error: Mat_Demo::mat_add: matrix size mismatch");
    }
    auto *m = new HwaUtil::Mat_Demo(a.nrows, a.ncols, Mat_Demo::MatrixType::Zero);
    for (int i = 0; i < a.nrows * a.ncols; i++) {
        m->d[i] = alpha * a.d[i] + beta * b.d[i];
    }
    Timer::tock("HwaUtil::Mat_Demo", "mat_add");
    return *m;
}

