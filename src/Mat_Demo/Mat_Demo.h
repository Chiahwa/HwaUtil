//
// Created by Chiahwa Young on 2023/3/17.
//

#ifndef HWAUTIL_MAT_DEMO_H
#define HWAUTIL_MAT_DEMO_H

#include <iostream>


namespace HwaUtil {
    class Mat_Demo {
    public:
        //Way to initialize a matrix.
        enum class MatrixType {
            Zero,
            Identity,
            Random
        };

        Mat_Demo();

        Mat_Demo(int nr, int nc, MatrixType initType = MatrixType::Zero);

        Mat_Demo(int nr, int nc, double *data);

        explicit Mat_Demo(std::istream &is);

        ~Mat_Demo();

        //functions to get the number of rows and columns.
        int nr() const;

        int nc() const;

        //functions to get the maximum and minimum values in the matrix.
        double mmax() const;

        double mmin() const;

        //functions to make the matrix all zeros.
        void zero();

        void zero(int nr, int nc);


        //Matrix manipulation.
        Mat_Demo &operator*=(const double &a);

        Mat_Demo &operator+=(const Mat_Demo &a);

        Mat_Demo &operator-=(const Mat_Demo &a);

        Mat_Demo &operator+(const Mat_Demo &a);

        Mat_Demo &operator-(const Mat_Demo &a);

        friend Mat_Demo &mat_add(const Mat_Demo &a, const Mat_Demo &b, double alpha, double beta);

        //Matrix multiplication.
        Mat_Demo &operator*(const Mat_Demo &a); //manual multiplication.
        Mat_Demo &blas_mult(const Mat_Demo &a); //using BLAS.

        int is_real_symm() const; //check if the matrix is real symmetric matrix(1) or not(0).

        int lapack_eig(double *eigval,
                       double *eigvec); //using LAPACK to calculate eigenvalues and eigenvectors, RSM assumed.Return 0 if success.

#ifdef __SCALAPACK__ //if SCALAPACK is available.

        int scalapack_eig(double *eigval,
                          double *eigvec); //using SCALAPACK to calculate eigenvalues and eigenvectors, RSM assumed.Return 0 if success.
#endif

        Mat_Demo &operator=(const Mat_Demo &a);

        //access the matrix elements.
        double &operator()(int i, int j);

        const double &operator()(int i, int j) const;

        const double *get_ptr(int i, int j); //TODO:not safe.


        //print the matrix.
        friend std::ostream &operator<<(std::ostream &os, const Mat_Demo &m);

    private:
        int nrows = 0;
        int ncols = 0;

        //holds data for the matrix.
        double *d = nullptr;

        static int idx_l2g(int idx_l, int block_size, int i_proc, int n_proc);
    };

    // Performs alpha * a + beta * b
    Mat_Demo &mat_add(const HwaUtil::Mat_Demo &a, const HwaUtil::Mat_Demo &b, double alpha, double beta);
}

#endif //HWAUTIL_MAT_DEMO_H
