//
// Created by Chiahwa Young on 2023/5/31.
//
/* miscellaneous */
#ifndef HWAUTIL_CLASSPROJECT_MISC_H
#define HWAUTIL_CLASSPROJECT_MISC_H

struct Point {
    double x, y, z;
};

struct RadialFunc{
    /* mesh = cutoff / dr + 1 */
    double cutoff; //cutoff radius.
    double dr; //step size.
    int mesh; //number of mesh points.

    double *r; //mesh points.
    double *v; //values of the radial function.
};

struct Func3d{
    //double cutoff; //cutoff radius.
    //double dr; //step size.
    int nx,ny,nz; //number of mesh points on each dimension.
    int xrange[2], yrange[2], zrange[2]; //range of mesh points on each dimension.
    int dx,dy,dz; //step size on each dimension.

    double *x; //mesh points.
    double *y; //mesh points.
    double *z; //mesh points.
    double *v; //values of the function.
};

#endif //HWAUTIL_CLASSPROJECT_MISC_H
