//
// Created by Chiahwa Young on 2023/5/31.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include <cmath>

#ifdef __OPENMP__
#include <omp.h>
#endif

#include "ArgumentReader/ArgumentReader.h"
#include "Mat_Demo/Mat_Demo.h"
#include "misc.h"
#include "Timer/Timer.h"

#ifdef __MPI__

#include <mpi.h>

#endif

using namespace std;
using namespace HwaUtil;

int proc_rank = 0, n_proc = 1;

/* global arguments */
HwaUtil::ArgumentReader arguments;
string working_dir;
char *working_dir_cstr;
int working_dir_len;

double lx, ly, lz;
int nx, ny, nz;
int n_points;
int n_f;
Point *points;
RadialFunc *f;
Func3d V;

double *h;
double *H;

void read_point(string points_path, Point *&pts, int n = 2);

void read_f(const string &path, RadialFunc *&ff, int n = 1);

void read_V(string path, Func3d &funcV);

void read_args(const string &input_file_path);

void broadcast();

void compute_H();

int main(int argc, char **argv) {

#ifdef __MPI__
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    if (proc_rank == 0) {
        cout << "MPI initialized. # of processes: " << n_proc << endl;
    }
#else
    cout<<"MPI not enabled."<<endl;
#endif
    Timer::init();
    Timer::tick("HwaUtil::(ClassProject)", "main");
    if (proc_rank == 0) {
        cout << "Hello World from process 0!\nCurrent working dir:" << endl;
        working_dir = std::filesystem::current_path().string() + "/";
        cout << working_dir << endl;

        /* (1) read arguments */
        string input_file_path; //path for argument file
        if (argc != 2) {
            cout << "Please input the path of argument file: " << endl;
            getline(cin, input_file_path);
        } else input_file_path = argv[1];

        //change working dir to the dir of input file
        working_dir = input_file_path.substr(0, input_file_path.find_last_of('/')) + "/";
        read_args(input_file_path);

        //cout << working_dir << endl;

        /* (2) read points */
        n_points = stoi(arguments.GetArgV("n_points"));
        string points_path = arguments.GetArgV("points_path");
        read_point(points_path, points, n_points);

        /* (3) read radial distribution function */
        string distribution_path = arguments.GetArgV("distribution_path");
        read_f(distribution_path, f);

        /* (4) read function V of space points */
        string v_path = arguments.GetArgV("v_path");
        read_V(v_path, V);
        nx = V.nx;
        ny = V.ny;
        nz = V.nz;

    }

    /* (5) broadcast arguments */
#ifdef __MPI__

    MPI_Barrier(MPI_COMM_WORLD);
    broadcast();
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    /* (6) calculate integral */
    if (proc_rank == 0) cout << "Start calculating integral..." << endl;

    compute_H();

    /* (7) reduce */
#ifdef __MPI__
    MPI_Barrier(MPI_COMM_WORLD);
    if (proc_rank == 0)
        H = new double[n_points * n_points];
    MPI_Reduce(h, H, n_points * n_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    H=h;
#endif
    if (proc_rank == 0) {
        cout << "H:" << endl;
        for (int i = 0; i < n_points; i++) {
            for (int j = 0; j < n_points; j++) {
                cout << H[i * n_points + j] << " ";
            }
            cout << endl;
        }
    }


    /* (8) diagonalize */
    if (proc_rank == 0) cout << "Diagonalizing..." << endl;

    if (proc_rank == 0) {
        if (arguments.GetArgV("diago_lib") == "lapack") {
            auto *mat = new Mat_Demo(n_points, n_points, H);
            auto *eig = new double[n_points];
            auto *eigv = new double[n_points * n_points];
            mat->lapack_eig(eig,
                            eigv);

            filesystem::create_directory(working_dir + "output");
            ofstream eig_file(working_dir + "output/eigenvalues.log");
            ofstream eigv_file(working_dir + "output/eigenvectors.log");

            eig_file << "Eigenvalues:" << endl;
            for (int i = 0; i < n_points; i++) {
                eig_file << eig[i] << " ";
            }
            eig_file << endl;
            eigv_file << "Eigenvectors:" << endl;
            for (int i = 0; i < n_points; i++) {
                for (int j = 0; j < n_points; j++) {
                    eigv_file << eigv[i * n_points + j] << " ";
                }
                eigv_file << endl;
            }
            eig_file.close();
            eigv_file.close();
            delete[] eig;
            delete[] eigv;
            delete mat;
        }
    }

    Timer::tock("HwaUtil::(ClassProject)", "main");

#ifdef __MPI__
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    string log_file_name = working_dir+"output/processor_" + to_string(proc_rank) + "_timer.log";
    ofstream log_file(log_file_name);
    HwaUtil::Timer::print_time_usage(log_file);
    log_file.close();

#ifdef __MPI__
    MPI_Finalize();
#endif
    return 0;
}


void read_point(string points_path, Point *&pts, int n) {
    Timer::tick("HwaUtil::(ClassProject)", "read_point");
    if (points_path[0] != '/')
        points_path = working_dir + points_path;
    cout << "Reading points file: " << points_path << "..." << endl;
    ifstream fs(points_path);
    if (fs.fail())
        throw invalid_argument("Invalid points file path");
    pts = new Point[n];
    for (int i = 0; i < n; i++) {
        fs.ignore(256, '(');
        fs >> pts[i].x;
        fs.ignore(256, ',');
        fs >> pts[i].y;
        fs.ignore(256, ',');
        fs >> pts[i].z;
        fs.ignore(256, ')');
    }
    fs.close();
    Timer::tock("HwaUtil::(ClassProject)", "read_point");
}

void read_f(const string &path, RadialFunc *&ff, int n) {
    Timer::tick("HwaUtil::(ClassProject)", "read_f");
    stringstream ss(path);
    n_f = n;
    ff = new RadialFunc[n];
    for (int i = 0; i < n; i++) {
        string pathi;
        if (i < n - 1)
            getline(ss, pathi, '|');
        else
            getline(ss, pathi);
        if (pathi[0] != '/')
            pathi = working_dir + pathi;
        cout << "Reading distribution functions file: " << pathi << "..." << endl;
        ifstream fs(pathi);
        if (fs.fail()) {
            Timer::tock("HwaUtil::(ClassProject)", "read_f");
            throw invalid_argument("Invalid distribution functions file path");
        }

        ArgumentReader ar;
        {
            ar.AddArg("cutoff");
            ar.AddArg("dr");
            ar.AddArg("mesh");
            ar.AddArg("l");
            ar.SetDataLabel("f");
        }
        ar.ReadArgs(fs);
        int l = stoi(ar.GetArgV("l"));
        ff[l - 1].cutoff = stod(ar.GetArgV("cutoff"));
        ff[l - 1].dr = stod(ar.GetArgV("dr"));
        ff[l - 1].mesh = stoi(ar.GetArgV("mesh"));

        //ff[l - 1].r = new double[ff[l - 1].mesh];
        ff[l - 1].v = new double[ff[l - 1].mesh];
        for (int j = 0; j < ff[l - 1].mesh; j++) {
            fs >> ff[l - 1].v[j];
            //    ff[l - 1].r[j] = ff[l - 1].dr * j;
            fs.ignore(256, ',');
        }
    }
    Timer::tock("HwaUtil::(ClassProject)", "read_f");
}

void read_V(string path, Func3d &funcV) {
    Timer::tick("HwaUtil::(ClassProject)", "read_V");
    if (path[0] != '/')
        path = working_dir + path;
    cout << "Reading function V file: " << path << "..." << endl;
    ifstream fs(path);
    if (fs.fail()) {
        Timer::tock("HwaUtil::(ClassProject)", "read_V");
        throw invalid_argument("Invalid V file path");
    }

    ArgumentReader ar;
    {
        ar.AddArg("nx");
        ar.AddArg("ny");
        ar.AddArg("nz");
        ar.SetDataLabel("V");
    }
    ar.ReadArgs(fs);

    funcV.nx = stoi(ar.GetArgV("nx"));
    funcV.ny = stoi(ar.GetArgV("ny"));
    funcV.nz = stoi(ar.GetArgV("nz"));
    funcV.xrange[0] = 0;
    funcV.xrange[1] = lx;
    funcV.yrange[0] = 0;
    funcV.yrange[1] = ly;
    funcV.zrange[0] = 0;
    funcV.zrange[1] = lz;
    funcV.dx = (funcV.xrange[1] - funcV.xrange[0]) / funcV.nx;
    funcV.dy = (funcV.yrange[1] - funcV.yrange[0]) / funcV.ny;
    funcV.dz = (funcV.zrange[1] - funcV.zrange[0]) / funcV.nz;
    funcV.v = new double[funcV.nx * funcV.ny * funcV.nz];
    //TODO: optimize. use 1d array instead of 3d array
    for (int i = 0; i < funcV.nx * funcV.ny * funcV.nz; i++)
        fs >> funcV.v[i];

    fs.close();
    Timer::tock("HwaUtil::(ClassProject)", "read_V");
}

void read_args(const string &input_file_path) {
    Timer::tick("HwaUtil::(ClassProject)", "read_args");
    ifstream fs(input_file_path);
    if (fs.fail()) {
        Timer::tock("HwaUtil::(ClassProject)", "read_args");
        throw invalid_argument("Invalid input file path");
    }
    cout << "Reading argument file: " << input_file_path << "..." << endl;
    {//configure argument reader
        arguments.AddArg("lx");
        arguments.AddArg("ly");
        arguments.AddArg("lz");
        arguments.AddArg("isHexahedral");
        arguments.AddArg("thetaxy");
        arguments.AddArg("thetayz");
        arguments.AddArg("thetaxz");
        arguments.AddArg("support_SH");
        arguments.AddArg("diago_lib");
        arguments.AddArg("support_Periodic_Boundary");
        arguments.AddArg("multi_parallel_strategies");
        arguments.AddArg("n_points");
        arguments.AddArg("points_path");
        arguments.AddArg("v_path");
        arguments.AddArg("distribution_path");
    }
    arguments.ReadArgs(fs);
    fs.close();

    lx = stoi(arguments.GetArgV("lx"));
    ly = stoi(arguments.GetArgV("ly"));
    lz = stoi(arguments.GetArgV("lz"));
    Timer::tock("HwaUtil::(ClassProject)", "read_args");
}

#ifdef __MPI__

void broadcast() {
    Timer::tick("HwaUtil::(ClassProject)", "broadcast");
    if (proc_rank == 0) {
        cout << "Transforming arguments..." << endl;
        working_dir_len = working_dir.length();
        working_dir_cstr = new char[working_dir_len + 1];
        strcpy(working_dir_cstr, working_dir.c_str());
    }
    MPI_Bcast(&n_points, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ly, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_f, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&working_dir_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (proc_rank != 0) {
        points = new Point[n_points];
        f = new RadialFunc[n_f];
        V.v = new double[nx * ny * nz];
        working_dir_cstr = new char[working_dir_len + 1];
    }
    MPI_Bcast(working_dir_cstr, working_dir_len + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (proc_rank != 0)
        working_dir = string(working_dir_cstr);

    if (proc_rank == 0) cout << "Broadcasting points..." << endl;
    for (int i = 0; i < n_points; i++) {
        MPI_Bcast(&points[i].x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&points[i].y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&points[i].z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (proc_rank == 0) cout << "Broadcasting function V..." << endl;
    MPI_Bcast(V.xrange, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(V.yrange, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(V.zrange, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(V.v, nx * ny * nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&V.dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&V.dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&V.dz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    //if(proc_rank==1) cout<<V.xrange[1]<<endl;
    if (proc_rank == 0) cout << "Broadcasting radial distribution functions..." << endl;

    for (int i = 0; i < n_f; i++) {
        MPI_Bcast(&f[i].cutoff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&f[i].mesh, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&f[i].dr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (proc_rank != 0) {
            f[i].v = new double[f[i].mesh];
        }
        MPI_Bcast(f[i].v, f[i].mesh, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    Timer::tock("HwaUtil::(ClassProject)", "broadcast");
}

#endif

void compute_H() {
    Timer::tick("HwaUtil::(ClassProject)", "compute_H");
    h = new double[n_points * n_points];
    int nx_start = nx / n_proc * proc_rank;
    int nx_end = nx / n_proc * (proc_rank + 1) > nx ? nx : nx / n_proc * (proc_rank + 1);
#ifdef __OPENMP__
    omp_set_num_threads(2);
#endif
    double fff[2];
    for (int i = nx_start; i < nx_end; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                for (int l = 0; l < n_points; l++) {
                    for (int m = l; m < n_points; m++) {
#ifdef __OPENMP__
#pragma omp parallel
                        {
                            int id=omp_get_thread_num();
                            if(id==0)   fff[id] = (*f)(sqrt(pow(points[l].x-i*lx/nx,2)+pow(points[l].y-j*ly/ny,2)+pow(points[l].z-k*lz/nz,2)));
                            else        fff[id] = (*f)(sqrt(pow(points[m].x-i*lx/nx,2)+pow(points[m].y-j*ly/ny,2)+pow(points[m].z-k*lz/nz,2)));
                        }

#else
                        fff[0] = (*f)(sqrt(pow(points[l].x - i * lx / nx, 2) + pow(points[l].y - j * ly / ny, 2) +
                                           pow(points[l].z - k * lz / nz, 2)));
                        if (fff[0] == 0) continue;
                        fff[1] = (*f)(sqrt(pow(points[m].x - i * lx / nx, 2) + pow(points[m].y - j * ly / ny, 2) +
                                           pow(points[m].z - k * lz / nz, 2)));
#endif
                        h[l * n_points + m] += V.v[i * ny * nz + j * nz + k] * fff[0] * fff[1] * V.dx * V.dy * V.dz;
                    }
                }
            }
        }
    }

    for (int i = 0; i < n_points; i++) {
        for (int j = i; j < n_points; j++) {
            h[j * n_points + i] = h[i * n_points + j];
        }
    }
    Timer::tock("HwaUtil::(ClassProject)", "compute_H");
}