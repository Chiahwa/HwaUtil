//
// Created by Chiahwa Young on 2023/3/17.
//

#ifndef HWAUTIL_MAT_DEMO_H
#define HWAUTIL_MAT_DEMO_H


class Mat_Demo {
public:
    Mat_Demo();
    Mat_Demo(
            const int nr,
            const int nc,
            const bool flag_zero = true);
    ~Mat_Demo();

    int nr() const;
    int nc() const;

    double max() const;
    double min() const;

    void zero();
    void zero(const int nr, const int nc);

    void operator *=(const double a);


private:
    int nrows = 0;
    int ncols = 0;

    //holds data for the matrix.
    double *d = nullptr;
};


#endif //HWAUTIL_MAT_DEMO_H
