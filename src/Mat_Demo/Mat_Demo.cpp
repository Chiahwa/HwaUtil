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
#include <cmath>

#if defined(__APPLE__) || defined(__MACOSX)

#include <Accelerate/Accelerate.h>

#elif defined(__linux__)

#include <cblas.h>
#include <lapacke.h>

#endif

#ifdef __SCALAPACK__

#include <mpi.h>

extern "C" {
// PDSYEVD function declaration
void pdsyevd_(char *jobz, char *uplo, int *n, double *a, int *ia, int *ja, int *desca, double *w, double *z, int *iz,
              int *jz, int *descz, double *work, int *lwork, int *iwork, int *liwork, int *info);

// CBLACS functions
void Cblacs_pinfo(int *myrank, int *nprocs);
void Cblacs_get(int context, int request, int *value);
void Cblacs_gridinit(int *context, char *order, int np_row, int np_col);
void Cblacs_gridinfo(int context, int *np_row, int *np_col, int *my_row, int *my_col);
void Cblacs_gridexit(int context);
void Cblacs_exit(int status);
void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt,
               int *lld, int *info);
int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
void pdgemr2d_(int *m, int *n, double *A, int *ia, int *ja, int *desca, double *B, int *ib, int *jb,
               int *descb, int *context);
}

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
    memcpy(eigvec, d, sizeof(double) * nrows * nrows);
#if defined(__APPLE__) || defined(__MACOSX)
    dsyev_("V", "U", &nrows, d, &nrows, eigval, work, &lwork, &info);
#else if defined(__linux__) || defined(__linux)
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nrows, eigvec, nrows, eigval);
#endif
    if (info != 0) {
        Timer::tock("HwaUtil::Mat_Demo", "lapack_eig");
        std::cerr << "Error: Mat_Demo::lapack_eig: LAPACK dsyev failed" << std::endl;
        return info;
    }
    /*for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < nrows; j++) {
            eigvec[i * nrows + j] = d[i * nrows + j];
        }
    }*/
    Timer::tock("HwaUtil::Mat_Demo", "lapack_eig");
    return 0;
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

#ifdef __SCALAPACK__

int HwaUtil::Mat_Demo::scalapack_eig(double *eigval, double *eigvec) {
    Timer::tick("HwaUtil::Mat_Demo", "scalapack_eig");
    // Constants
    int ZERO = 0;
    int ONE = 1;
    int NEGONE = -1;

    int my_rank;          // 当前进程的rank
    int num_procs;        // 进程总数
    int blacs_context;    // CBLACS上下文

    // 获取当前进程的rank和进程总数
    Cblacs_pinfo(&my_rank, &num_procs);
    //std::cout<<my_rank<<": init blacs"<<std::endl;

    Cblacs_get(-1, 0, &blacs_context);
    //std::cout<<blacs_context<<std::endl;


    // 设置进程网格
    int nprow = (int) sqrt(num_procs);
    int npcol = num_procs / nprow;
    while (nprow * npcol != num_procs) {
        nprow++;
        npcol = num_procs / nprow;
    }

    char order = 'R';
    Cblacs_gridinit(&blacs_context, &order, nprow, npcol);
    int my_row, my_col;
    Cblacs_gridinfo(blacs_context, &nprow, &npcol, &my_row, &my_col);

    // 设置当前进程的局部矩阵
    int block_size = 2; //TODO: make it a parameter
    int nrow_loc = numroc_(&nrows, &block_size, &my_row, &ZERO, &nprow);
    int ncol_loc = numroc_(&ncols, &block_size, &my_col, &ZERO, &npcol);

    int lld_loc = std::max(1, nrow_loc);

    //std::cout<<"Process "<<my_rank<<", nrow_loc: "<<nrow_loc<<std::endl;

    double *a_loc = new double[nrow_loc * ncol_loc]; // 存储局部矩阵
    for (int i = 0; i < nrow_loc; i++) {
        int i_glb = idx_l2g(i, block_size, my_row, nprow);
        for (int j = 0; j < ncol_loc; j++) {
            int j_lglb = idx_l2g(j, block_size, my_col, npcol);
            a_loc[i + j * nrow_loc] = d[i_glb + j_lglb * nrows];
        }
    }

    double *z_loc = new double[nrows * ncols]; // 存储局部特征向量

    // 获得desc
    int desc_a[9];
    int info;
    descinit_(desc_a, &nrows, &ncols, &block_size, &block_size,
              &ZERO, &ZERO, &blacs_context, &lld_loc, &info);
    if (info != 0) {
        Timer::tock("HwaUtil::Mat_Demo", "scalapack_eig");
        throw std::runtime_error("Error: Mat_Demo::scalapack_eig: descinit_ failed");
    }

    int desc_z[9];
    descinit_(desc_z, &nrows, &ncols, &block_size, &block_size,
              &ZERO, &ZERO, &blacs_context, &lld_loc, &info);
    if (info != 0) {
        Timer::tock("HwaUtil::Mat_Demo", "scalapack_eig");
        throw std::runtime_error("Error: Mat_Demo::scalapack_eig: descinit_ failed");
    }

    //std::cout<<"got desc"<<std::endl;


    // 使用pdsyevd_计算特征值和特征向量
    int lwork = -1;
    int liwork = -1;
    double *work = new double[1];
    int *iwork = new int[1];

    char jobz = 'V'; // 计算特征值和特征向量
    char uplo = 'U'; // 上三角存储
    // double *w = new double[nrow_loc];

    pdsyevd_(&jobz, &uplo, &nrows, a_loc, &ONE, &ONE, desc_a,
             eigval, z_loc, &ONE, &ONE, desc_z,
             work, &lwork, iwork, &liwork,
             &info);
    if (info != 0) {
        Timer::tock("HwaUtil::Mat_Demo", "scalapack_eig");
        throw std::runtime_error("Error: Mat_Demo::scalapack_eig: first trial of pdsyevd_ failed");
    }

    lwork = (int) work[0];
    liwork = iwork[0];
    delete[] work;
    delete[] iwork;
    work = new double[lwork];
    iwork = new int[liwork];
    pdsyevd_(&jobz, &uplo, &nrows, a_loc, &ONE, &ONE, desc_a,
             eigval, z_loc, &ONE, &ONE, desc_z,
             work, &lwork, iwork, &liwork,
             &info);
    if (info != 0) {
        Timer::tock("HwaUtil::Mat_Demo", "scalapack_eig");
        throw std::runtime_error("Error: Mat_Demo::scalapack_eig: second trial of pdsyevd_ failed");
    }

    // 将局部的特征值和特征向量收集到主进程
    const int TAG_Z = 0;
    const int TAG_NROW_LOC = 1;
    const int TAG_NCOL_LOC = 2;
    const int TAG_MY_ROW = 3;
    const int TAG_MY_COL = 4;

    if (my_rank != 0) {
        MPI_Send(z_loc, nrows * ncols, MPI_DOUBLE, 0, TAG_Z, MPI_COMM_WORLD);
        MPI_Send(&nrow_loc, 1, MPI_INT, 0, TAG_NROW_LOC, MPI_COMM_WORLD);
        MPI_Send(&ncol_loc, 1, MPI_INT, 0, TAG_NCOL_LOC, MPI_COMM_WORLD);
        MPI_Send(&my_row, 1, MPI_INT, 0, TAG_MY_ROW, MPI_COMM_WORLD);
        MPI_Send(&my_col, 1, MPI_INT, 0, TAG_MY_COL, MPI_COMM_WORLD);
    }
    if (my_rank == 0) {
        double *z_buf = new double[nrows * ncols];
        int nrow_loc_buf, ncol_loc_buf, my_row_buf, my_col_buf;
        for (int i = 0; i < num_procs; ++i) {
            if (i == 0) {
                memcpy(z_buf, z_loc, nrows * ncols * sizeof(double));
                nrow_loc_buf = nrow_loc;
                ncol_loc_buf = ncol_loc;
                my_row_buf = my_row;
                my_col_buf = my_col;
            } else {
                MPI_Recv(z_buf, nrows * ncols, MPI_DOUBLE, i, TAG_Z, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&nrow_loc_buf, 1, MPI_INT, i, TAG_NROW_LOC, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&ncol_loc_buf, 1, MPI_INT, i, TAG_NCOL_LOC, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&my_row_buf, 1, MPI_INT, i, TAG_MY_ROW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&my_col_buf, 1, MPI_INT, i, TAG_MY_COL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            for (int j = 0; j < nrow_loc_buf; ++j) {
                int i_glb = idx_l2g(j, block_size, my_row_buf, nprow);
                for (int k = 0; k < ncol_loc_buf; ++k) {
                    int j_glb = idx_l2g(k, block_size, my_col_buf, npcol);
                    eigvec[i_glb + j_glb * nrows] = z_buf[j + k * nrow_loc_buf];
                }
            }
        }
        delete[] z_buf;
    }

    MPI_Barrier(MPI_COMM_WORLD); // TODO: MPI 的调用是否具有通用性？

    delete[] a_loc;
    delete[] z_loc;
    delete[] work;
    delete[] iwork;


    Cblacs_gridexit(blacs_context);

    Timer::tock("HwaUtil::Mat_Demo", "scalapack_eig");
    return 0;
}

// 将局部矩阵的行列号转换为全局矩阵的行列号
int HwaUtil::Mat_Demo::idx_l2g(int idx_l, int block_size, int i_proc, int n_proc) {
    return (idx_l / block_size * n_proc * block_size // 完整的周期
            + i_proc * block_size                 // 周期内的进程块偏移
            + idx_l % block_size    // 进程块内的偏移
    );
}

#endif