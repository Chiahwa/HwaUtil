//
// Created by Chiahwa Young on 2023/5/31.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include "ArgumentReader/ArgumentReader.h"
#include "misc.h"

using namespace std;
using namespace HwaUtil;

/* global arguments */
HwaUtil::ArgumentReader arguments;
string working_dir;
int lx, ly, lz;
int nx, ny, nz;
int n_points;
Point *points;
RadialFunc *f;
Func3d V;


void read_point(string points_path, Point *pts, int n = 2) {
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
}

void read_f(const string &path, RadialFunc * &ff, int n = 1) {
    stringstream ss(path);

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
        if (fs.fail())
            throw invalid_argument("Invalid distribution functions file path");

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

        ff[l - 1].r = new double[ff[l - 1].mesh];
        ff[l - 1].v = new double[ff[l - 1].mesh];
        for (int j = 0; j < ff[l - 1].mesh; j++) {
            fs >> ff[l - 1].v[j];
            ff[l - 1].r[j] = ff[l - 1].dr * j;
            fs.ignore(256, ',');
        }
    }
}

void read_V(string path, Func3d &funcV) {
    if (path[0] != '/')
        path = working_dir + path;
    cout << "Reading function V file: " << path << "..." << endl;
    ifstream fs(path);
    if (fs.fail())
        throw invalid_argument("Invalid V file path");

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
    funcV.x = new double[funcV.nx];
    funcV.y = new double[funcV.ny];
    funcV.z = new double[funcV.nz];
    for (int i = 0; i < funcV.nx; i++)
        funcV.x[i] = funcV.xrange[0] + funcV.dx * i;
    for (int i = 0; i < funcV.ny; i++)
        funcV.y[i] = funcV.yrange[0] + funcV.dy * i;
    for (int i = 0; i < funcV.nz; i++)
        funcV.z[i] = funcV.zrange[0] + funcV.dz * i;

}

void read_args(const string &input_file_path) {
    ifstream fs(input_file_path);
    if (fs.fail())
        throw invalid_argument("Invalid input file path");
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
}

int main(int argc, char **argv) {
    cout << "Hello World!\nCurrent working dir:" << endl;
    working_dir = std::filesystem::current_path().string() + "/";
    cout << working_dir << endl;

    /* (1) read arguments */
    string input_file_path; //path for argument file
    if (argc != 2) {
        cout << "Please input the path of argument file: " << endl;
        getline(cin, input_file_path);
    } else input_file_path = argv[1];
    /*
    if (input_file_path[0] != '/') //relative path
        working_dir = working_dir + input_file_path.substr(0, input_file_path.find_last_of('/')) + "/";
    else    //absolute path
     */

    //change working dir to the dir of input file
    working_dir = input_file_path.substr(0, input_file_path.find_last_of('/')) + "/";
    read_args(input_file_path);

    cout<<working_dir<<endl;

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

    cout<<f->evallinear(20)<<endl;
    cout<<(*f)(20)<<endl;

    return 0;
}