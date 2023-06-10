//
// Created by Chiahwa Young on 2023/5/31.
//
/* miscellaneous */
#ifndef HWAUTIL_CLASSPROJECT_MISC_H
#define HWAUTIL_CLASSPROJECT_MISC_H

extern "C" {
void Cblacs_pinfo(int*, int*);
void Cblacs_get(int, int, int*);
void Cblacs_gridinit(int*, char const*, int, int);
void Cblacs_gridinfo(int, int*, int*, int*, int*);
void Cblacs_barrier(int , char*);
void Cblacs_gridexit(int);
void Cblacs_exit(int);

int numroc_(int const& n, int const& nb, int const& iproc, int const& isproc, int const& nprocs);
int indxg2p_(int const& glob, int const& nb, int const& iproc, int const& isproc, int const& nprocs);
int indxl2g_(int const& loc, int const& nb, int const& iproc, int const& isproc, int const& nprocs);
void descinit_( int *desc, int const& m, int const& n, int const& mb, int const& nb, int const& irsrc, int const& icsrc, int const& ictxt, int const& lld, int *info);
void pdgesv_( int const& n, int const& nrhs, double *A, int const& ia, int const& ja, int *desca, int* ipiv, double *B, int const& ib, int const& jb, int *descb, int *info);
void pdgemm_( char const *transa, char const *transb, int const& M, int const& N, int const& K, double const& ALPHA,  double * A, int const& IA, int const& JA, int * DESCA, double * B, int const& IB, int const& JB, int * DESCB, double const& BETA, double * C, int const& IC, int const& JC, int * DESCC );
double pdlange_( char const *norm, int const& m, int const& n, double *A, int const& ia, int const& ja, int *desca, double *work);
}


struct Point {
    double x, y, z;
};

struct RadialFunc{
    /* mesh = cutoff / dr + 1 */
    double cutoff; //cutoff radius.
    double dr; //step size.
    int mesh; //number of mesh points.

    //double *r; //mesh points.
    double *v; //values of the radial function.

    /* use cubic spline interpolation to get the value of the radial function at rr */
    double operator()(double rr) const {

        int i = (int) (rr / dr);
        if (i < 0) return v[0];
        if (i >= mesh - 1) return 0;
        double a = (v[i + 1] - v[i]) / dr;
        double b = v[i];
        double c = (v[i + 1] - v[i]) / dr / dr / 2;
        double d = (v[i] - v[i + 1]) / dr / dr / dr;
        double x = rr - i * dr;
        return a * x + b + c * x * x + d * x * x * x;

        /*int n = mesh - 1;
        double h = cutoff / n;
        int i = static_cast<int>(rr / h);
        if (i < 0) return v[0];
        if (i >= mesh - 1) return 0;

        double t = (rr - i * h) / h;
        double t2 = t * t;
        double t3 = t2 * t;

        double f0, f1, d0, d1;

        if (i == 0) {
            f0 = v[0];
            d0 = (v[1] - v[0]) / h;
            f1 = v[1];
            d1 = (v[2] - v[1]) / h;
        } else {
            f0 = v[i - 1];
            d0 = (v[i] - v[i - 1]) / h;
            f1 = v[i];
            d1 = (v[i + 1] - v[i]) / h;
        }

        double a = 2 * (f0 - f1) + (d0 + d1);
        double b = -3 * (f0 - f1) - 2 * d0 - d1;
        double c = d0;
        double d = f0;

        return a * t3 + b * t2 + c * t + d;*/
    }

    /* use linear interpolation to get the value of the radial function at rr */
    [[nodiscard]] double evallinear(double rr) const{
        int i = (int) (rr / dr);
        if (i < 0) return v[0];
        if (i >= mesh - 1) return v[mesh - 1];
        double a = (v[i + 1] - v[i]) / dr;
        double b = v[i];
        double x = rr - i * dr;
        return a * x + b;
    }
};

struct Func3d{
    //double cutoff; //cutoff radius.
    //double dr; //step size.
    int nx,ny,nz; //number of mesh points on each dimension.
    double xrange[2], yrange[2], zrange[2]; //range of mesh points on each dimension.
    double dx,dy,dz; //step size on each dimension.

    //double *x; //mesh points.
    //double *y; //mesh points.
    //double *z; //mesh points.
    double *v; //values of the function.
};

#endif //HWAUTIL_CLASSPROJECT_MISC_H
