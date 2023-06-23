/*#include <iostream>
#include <cstddef>
#include <cmath>
#include <iomanip>
#include <mpi.h>

// Declare Scalapack and CBLACS functions using extern "C"
extern "C" {
void Cblacs_pinfo(int *, int *);
void Cblacs_get(int, int, int *);
void Cblacs_gridinit(int *, const char *, int, int);
void Cblacs_gridinfo(int, int *, int *, int *, int *);
void Cblacs_gridexit(int);
void pdgemr2d(int, int, const double *, int, double *, int, int, int);
void pdsyevx_(char *, char *, char *, int *, double *, int *, int *, double *, double *, double *, int *, double *, double *,
         int *, int *, int *, int *);
int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
}*/


/*
// Function to initialize Scalapack context
void init_scalapack(int &context, int &num_procs, int &my_rank) {
    int ictxt, nprow, npcol, myrow, mycol;

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Create a 2D process grid
    Cblacs_pinfo(&my_rank, &num_procs);
    Cblacs_get(0, 0, &ictxt);
    Cblacs_gridinit(&ictxt, "Row-major", 1, 1);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    // Set the Scalapack context
    context = ictxt;
}

// Function to finalize Scalapack context
void finalize_scalapack(int &context) {
    int iam;

    // Get the process ID
    Cblacs_pinfo(&iam, NULL);

    // Finalize the process grid
    Cblacs_gridexit(context);

    // Finalize the MPI environment
    MPI_Finalize();
}

// Function to generate a symmetric matrix
void generate_matrix(double *matrix, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            matrix[i * n + j] = matrix[j * n + i] = std::rand() % 10;
        }
    }
}

int main() {
    int context, num_procs, my_rank;
    int n = 4;  // Size of the matrix
    int info;

    // Initialize Scalapack context
    init_scalapack(context, num_procs, my_rank);

    // Determine the size of the local matrices
    int nb = 2;  // Block size
    int n_loc = numroc_(&n, &nb, &my_rank, &my_rank, &num_procs);

    // Allocate memory for the local matrices
    double *matrix = new double[n_loc * n];

    // Generate the symmetric matrix
    generate_matrix(matrix, n_loc);

    if(my_rank == 0) {
        // Print the matrix
        std::cout << "Matrix:" << std::endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n_loc; ++j) {
                std::cout << matrix[i * n_loc + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // Call dsyevx to compute the eigenvalues and eigenvectors
    char jobz = 'V';  // Compute both eigenvalues and eigenvectors
    char range = 'A'; // Compute all eigenvalues in the interval (vl, vu)
    char uplo = 'U';  // Upper triangle of input matrix is stored
    double vl = 0.0;  // Lower bound of the interval
    double vu = 10.0; // Upper bound of the interval
    int il = 1;       // Index of the smallest eigenvalue to be returned
    int iu = n;       // Index of the largest eigenvalue to be returned
    double abstol = 0.0;  // Absolute tolerance for eigenvalues
    int m;  // Number of eigenvalues found
    double *eigenvalues = new double[n];
    double *eigenvectors = new double[n_loc * n];


    if (my_rank == 0) {
        std::cout << "Reached point before pdsyevx" << std::endl;
    }

    // Call dsyevx
    pdsyevx_(&jobz, &range, &uplo, &n, matrix, &il, &iu, &vl, &vu, &abstol,
             &m, eigenvalues, eigenvectors, &n, &n_loc, &context, &info);

    // Print the eigenvalues and eigenvectors
    if (my_rank == 0) {
        std::cout << "Eigenvalues: ";
        for (int i = 0; i < m; ++i) {
            std::cout << eigenvalues[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Eigenvectors:" << std::endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n_loc; ++j) {
                std::cout << eigenvectors[i * n_loc + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // Free the allocated memory
    delete[] matrix;
    delete[] eigenvalues;
    delete[] eigenvectors;

    // Finalize Scalapack context
    finalize_scalapack(context);

    return 0;
}*/


/*
void printMatrix(const double* matrix, int n, int nprow, int npcol, int myrow, int mycol) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int desc[9];
    int info;
    double* localMatrix;
    int localRows, localCols, localLld;

    if (myrow == 0 && mycol == 0) {
        std::cout << "Matrix:" << std::endl;
    }

    descinit_(desc, &n, &n, &n, &n, &myrow, &mycol, &comm, &nprow, &npcol, &info);
    pdescinit_(desc + 1, &n, &n, &n, &n, &myrow, &mycol, &comm, &nprow, &npcol, &info);

    if (myrow == 0 && mycol == 0) {
        localMatrix = new double[n * n];
    } else {
        localMatrix = new double[1];
    }

    pglue2p_(matrix, desc, localMatrix, desc + 1, &n, &n, &info);
    localRows = numroc_(&n, &n, &myrow, &ZERO, &nprow);
    localCols = numroc_(&n, &n, &mycol, &ZERO, &npcol);
    localLld = std::max(1, localRows);

    if (myrow == 0 && mycol == 0) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cout << std::setw(10) << localMatrix[j * localLld + i] << " ";
            }
            std::cout << std::endl;
        }
        delete[] localMatrix;
    }

    delete[] localMatrix;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int n = 4; // 矩阵的维度
    int nprow = 2; // 行进程数
    int npcol = 2; // 列进程数
    int myrow, mycol; // 当前进程的行和列索引
    int info, lwork, liwork;
    int m, nz, isuppz[2 * n];
    double abstol = 0.0;
    double* matrix;
    double* eigenvalues;
    double* eigenvectors;
    double* work;
    int* iwork;

    // 获取当前进程的行和列索引
    Cblacs_pinfo(&myrow, &mycol);

    // 创建网格通信环境
    Cblacs_get(0, 0, &GRID);
    Cblacs_gridinit(&GRID, "Row-major", nprow, npcol);

    // 分配矩阵
    int localRows = numroc_(&n, &n, &myrow, &ZERO, &nprow);
    int localCols = numroc_(&n, &n, &mycol, &ZERO, &npcol);
    int localLld = std::max(1, localRows);
    matrix = new double[localLld * localCols];

    // 初始化矩阵（示例：对角线为1）
    for (int i = 0; i < localCols; ++i) {
        for (int j = 0; j < localRows; ++j) {
            if (i == j) {
                matrix[i * localLld + j] = 1.0;
            } else {
                matrix[i * localLld + j] = 0.0;
            }
        }
    }

    // 打印原始矩阵
    printMatrix(matrix, n, nprow, npcol, myrow, mycol);

    // 获取工作空间大小
    lwork = -1;
    liwork = -1;
    eigenvalues = new double[n];
    eigenvectors = new double[localLld * localCols];
    dsyevx_("V", "A", "L", &n, matrix, &ONE, &ONE, desc, &ZERO, &ZERO, &ZERO, &ZERO, &abstol, &m, &nz, eigenvalues, eigenvectors, &ONE, &ONE, desc, work, &lwork, iwork, &liwork, isuppz, &info);

    lwork = work[0];
    liwork = iwork[0];
    delete[] work;
    delete[] iwork;
    work = new double[lwork];
    iwork = new int[liwork];

    // 计算特征值和特征向量
    dsyevx_("V", "A", "L", &n, matrix, &ONE, &ONE, desc, &ZERO, &ZERO, &ZERO, &ZERO, &abstol, &m, &nz, eigenvalues, eigenvectors, &ONE, &ONE, desc, work, &lwork, iwork, &liwork, isuppz, &info);

    // 打印特征值和特征向量
    if (myrow == 0 && mycol == 0) {
        std::cout << "Eigenvalues:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << eigenvalues[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Eigenvectors:" << std::endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cout << std::setw(10) << eigenvectors[j * localLld + i] << " ";
            }
            std::cout << std::endl;
        }
    }

    // 释放资源
    delete[] matrix;
    delete[] eigenvalues;
    delete[] eigenvectors;
    delete[] work;
    delete[] iwork;

    // 终止网格通信环境
    Cblacs_gridexit(GRID);

    MPI_Finalize();

    return 0;
}*/






/*#include <iostream>
#include <vector>
#include <cstddef>
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/eigenvalue_problem/syevx.hpp>
#include <mkl_blacs.h>

extern "C" {
// CBLACS functions
void Cblacs_pinfo(int* myrank, int* nprocs);
void Cblacs_get(int context, int request, int* value);
void Cblacs_gridinit(int* context, const char* order, int np_row, int np_col);
void Cblacs_gridinfo(int context, int* np_row, int* np_col, int* my_row, int* my_col);
void Cblacs_gridexit(int context);
}

int main() {
    // Initialize CBLACS
    int num_procs, my_rank;
    Cblacs_pinfo(&my_rank, &num_procs);

    int blacs_context;
    Cblacs_get(-1, 0, &blacs_context);
    Cblacs_gridinit(&blacs_context, "Row-major", 1, num_procs);

    // Matrix size and process grid information
    int n = 6; // Size of the matrix
    int np_row, np_col, my_row, my_col;
    Cblacs_gridinfo(blacs_context, &np_row, &np_col, &my_row, &my_col);

    // Create a block-cyclic matrix descriptor
    scalapack::BlockCyclicMatrix<double> A(n, n, blacs_context, np_row, np_col);

    // Initialize the matrix A with your own data
    // For simplicity, let's assume each process initializes a local block of A
    for (int i = 0; i < A.mlocal(); ++i) {
        for (int j = 0; j < A.nlocal(); ++j) {
            // Assign some values to the local block of A
            A.data(i, j) = my_row * A.mlocal() + i + 1;
        }
    }

    // Perform the eigenvalue computation
    std::vector<double> eigenvalues(n);
    std::vector<double> eigenvectors(n * n);

    scalapack::syevx(scalapack::Job::EigenvaluesVectors, scalapack::Range::All,
                     scalapack::Uplo::Upper, A, 0.0, 0.0, 1, n, 1, n, eigenvalues.data(),
                     eigenvectors.data(), n, blacs_context);

    // Print the eigenvalues and eigenvectors
    if (my_rank == 0) {
        std::cout << "Eigenvalues:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << eigenvalues[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Eigenvectors:" << std::endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cout << eigenvectors[i * n + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // Finalize CBLACS
    Cblacs_gridexit(blacs_context);

    return 0;
}*/

///这段代码使用了scalapackpp11库（https://github.com/dzyla/scalapackpp）来简化与Scalapack和CBLACS的交互。你可以在代码中修改矩阵的大小（变量`n`）和进程网格的大小（变量`np_row`和`np_col`）来适应你的情况。






#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>

extern "C" {
// Scalapack functions
void blacs_pinfo_(int* myrank, int* nprocs);
void blacs_get_(int* context, int* request, int* value);
void blacs_gridinit_(int* context, const char* order, int* np_row, int* np_col);
void blacs_gridinfo_(int* context, int* np_row, int* np_col, int* my_row, int* my_col);
void blacs_gridexit_(int* context);
void pdsyevx_(const char* jobz, const char* range, const char* uplo, int* n, double* A, int* ia,
              int* ja, int* desca, double* vl, double* vu, int* il, int* iu, double* abstol,
              int* m, double* w, double* Z, int* iz, int* jz, int* descz, double* work, int* lwork,
              int* iwork, int* liwork, int* ifail, int* info);
void descinit_(int* desc, int* m, int* n, int* mb, int* nb, int* irsrc, int* icsrc, int* ictxt,
               int* lld, int* info);
int numroc_(int* n, int* nb, int* iproc, int* isrcproc, int* nprocs);
void pdgemr2d_(int* m, int* n, double* A, int* ia, int* ja, int* desca, double* B, int* ib, int* jb,
               int* descb, int* context);
}

int main() {
    // Initialize MPI
    MPI_Init(nullptr, nullptr);

    // Get MPI rank and size
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Initialize BLACS
    int num_procs, my_rank;
    blacs_pinfo_(&my_rank, &num_procs);

    int blacs_context;
    blacs_get_(&blacs_context, nullptr, nullptr);

    int np_row, np_col, my_row, my_col;
    blacs_gridinit_(&blacs_context, "Row-major", &np_row, &np_col);
    blacs_gridinfo_(&blacs_context, &np_row, &np_col, &my_row, &my_col);

    // Matrix size and local block size
    int n = 6; // Size of the matrix
    int nb = 2; // Local block size

    // Create a local block of the matrix
    int mlocal = numroc_(&n, &nb, &my_row, &my_col, &np_row);

    std::vector<double> local_block(mlocal * n, 0.0);

    // Initialize the local block of the matrix with your own data
    for (int i = 0; i < mlocal; ++i) {
        for (int j = 0; j < n; ++j) {
            // Assign some values to the local block
            local_block[i * n + j] = (my_row * mlocal + i + 1) * std::pow(10, my_col);
        }
    }

    // Create a descriptor for the distributed matrix
    int desc[9];
    int info; // Variable to store the return value of Scalapack functions
    descinit_(desc, &n, &n, &nb, &nb, &my_row, &my_col, &blacs_context, &mlocal, &info);

    // Perform the eigenvalue computation
    char jobz = 'V'; // Compute both eigenvalues and eigenvectors
    char range = 'A'; // All eigenvalues are computed
    char uplo = 'U'; // Upper triangle of the matrix is stored
    int ia = 1; // Starting row index of the submatrix
    int ja = 1; // Starting column index of the submatrix
    double vl = 0.0; // Lower bound of eigenvalues
    double vu = 0.0; // Upper bound of eigenvalues
    int il = 0; // Lower bound of eigenvalues index
    int iu = 0; // Upper bound of eigenvalues index
    double abstol = 0.0; // Absolute tolerance
    int m; // Number of eigenvalues found
    std::vector<double> eigenvalues(n);
    std::vector<double> eigenvectors(mlocal * n);

    int iz = 1; // Starting column index of the eigenvector matrix
    int jz = 1; // Starting row index of the eigenvector matrix
    int descz[9]; // Descriptor for the eigenvector matrix
    descinit_(descz, &n, &n, &nb, &nb, &my_row, &my_col, &blacs_context, &mlocal, &info);

    int lwork = -1; // Query the optimal workspace size
    std::vector<double> work(1);
    std::vector<int> iwork(5 * n);
    std::vector<int> ifail(n);
    int liwork = 5 * n;



    pdsyevx_(&jobz, &range, &uplo, &n, local_block.data(), &ia, &ja, desc, &vl, &vu, &il, &iu, &abstol,
             &m, eigenvalues.data(), eigenvectors.data(), &iz, &jz, descz, work.data(), &lwork,
             iwork.data(), &liwork, ifail.data(), &info);

    // Obtain the optimal workspace size
    lwork = static_cast<int>(work[0]);
    work.resize(lwork);

    // Perform the eigenvalue computation with the correct workspace size
    pdsyevx_(&jobz, &range, &uplo, &n, local_block.data(), &ia, &ja, desc, &vl, &vu, &il, &iu, &abstol,
             &m, eigenvalues.data(), eigenvectors.data(), &iz, &jz, descz, work.data(), &lwork,
             iwork.data(), &liwork, ifail.data(), &info);

    // Gather the distributed eigenvectors to the root process
    std::vector<double> gathered_eigenvectors(n * n, 0.0);
    MPI_Gather(eigenvectors.data(), mlocal * n, MPI_DOUBLE, gathered_eigenvectors.data(),
               mlocal * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Print the eigenvalues and eigenvectors on the root process
    if (world_rank == 0) {
        std::cout << "Eigenvalues:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << eigenvalues[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Eigenvectors:" << std::endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cout << gathered_eigenvectors[i * n + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // Release BLACS resources
    blacs_gridexit_(&blacs_context);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}

/*
int main() {
    // Initialize BLACS
    int num_procs, my_rank;
    blacs_pinfo_(&my_rank, &num_procs);

    int blacs_context;
    blacs_get_(&blacs_context, nullptr, nullptr);

    int np_row, np_col, my_row, my_col;
    blacs_gridinit_(&blacs_context, "Row-major", &np_row, &np_col);
    blacs_gridinfo_(&blacs_context, &np_row, &np_col, &my_row, &my_col);

    // Matrix size and local block size
    int n = 6; // Size of the matrix
    int nb = 2; // Local block size

    // Create a local block of the matrix
    int mlocal = numroc_(&n, &nb, &my_row, &my_col, &np_row);

    std::vector<double> local_block(mlocal * n, 0.0);

    // Initialize the local block of the matrix with your own data
    for (int i = 0; i < mlocal; ++i) {
        for (int j = 0; j < n; ++j) {
            // Assign some values to the local block
            local_block[i * n + j] = (my_row * mlocal + i + 1) * std::pow(10, my_col);
        }
    }

    // Create a descriptor for the distributed matrix
    int desc[9];
    int info; // Variable to store the return value of Scalapack functions
    descinit_(desc, &n, &n, &nb, &nb, &my_row, &my_col, &blacs_context, &mlocal, &info);

    // Perform the eigenvalue computation
    char jobz = 'V'; // Compute both eigenvalues and eigenvectors
    char range = 'A'; // All eigenvalues are computed
    char uplo = 'U'; // Upper triangle of the matrix is stored
    int ia = 1; // Starting row index of the submatrix
    int ja = 1; // Starting column index of the submatrix
    double vl = 0.0; // Lower bound of eigenvalues
    double vu = 0.0; // Upper bound of eigenvalues
    int il = 0; // Lower bound of eigenvalues index
    int iu = 0; // Upper bound of eigenvalues index
    double abstol = 0.0; // Absolute tolerance
    int m; // Number of eigenvalues found
    std::vector<double> eigenvalues(n);
    std::vector<double> eigenvectors(mlocal * n);

    int iz = 1; // Starting column index of the eigenvector matrix
    int jz = 1; // Starting row index of the eigenvector matrix
    int descz[9]; // Descriptor for the eigenvector matrix
    descinit_(descz, &n, &n, &nb, &nb, &my_row, &my_col, &blacs_context, &mlocal, &info);

    int lwork = -1; // Query the optimal workspace size
    std::vector<double> work(1);
    std::vector<int> iwork(5 * n);
    std::vector<int> ifail(n);
    int liwork = 5 * n;



    pdsyevx_(&jobz, &range, &uplo, &n, local_block.data(), &ia, &ja, desc, &vl, &vu, &il, &iu, &abstol,
             &m, eigenvalues.data(), eigenvectors.data(), &iz, &jz, descz, work.data(), &lwork,
             iwork.data(), &liwork, ifail.data(), &info);

    // Obtain the optimal workspace size
    lwork = static_cast<int>(work[0]);
    work.resize(lwork);

    // Perform the eigenvalue computation with the correct workspace size
    pdsyevx_(&jobz, &range, &uplo, &n, local_block.data(), &ia, &ja, desc, &vl, &vu, &il, &iu, &abstol,
             &m, eigenvalues.data(), eigenvectors.data(), &iz, &jz, descz, work.data(), &lwork,
             iwork.data(), &liwork, ifail.data(), &info);

    // Create a global matrix for storing the distributed eigenvectors
    std::vector<double> gathered_eigenvectors(n * n, 0.0);

    // Gather the distributed eigenvectors to the root process
    pdgemr2d_(&n, &n, eigenvectors.data(), &iz, &jz, descz, gathered_eigenvectors.data(), &ia, &ja,
              desc, &blacs_context);

    // Print the eigenvalues and eigenvectors on the root process
    if (my_rank == 0) {
        std::cout << "Eigenvalues:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << eigenvalues[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Eigenvectors:" << std::endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cout << gathered_eigenvectors[i * n + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // Release BLACS resources
    blacs_gridexit_(&blacs_context);

    return 0;
}*/
