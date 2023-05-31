//
// Created by Chiahwa Young on 2023/5/31.
//

#include <iostream>
#include <fstream>
#include "ArgumentReader/ArgumentReader.h"

using namespace std;

/* global arguments */
int lx,ly,lz;
int nx,ny,nz;

int main(int argc, char **argv) {
    cout << "Hello World!" << endl;
    //path for argument file
    string input_file_path="input/INPUT.txt";
    if(argc != 2)
        getline(cin, input_file_path);
    else input_file_path = argv[1];

    ifstream fs(input_file_path);
    if(fs.fail())
        throw invalid_argument("Invalid input file path");
    cout<<"Reading argument file:"<<input_file_path<<"..."<<endl;

    HwaUtil::ArgumentReader argReader;
    argReader.AddArg("lx");
    argReader.AddArg("ly");
    argReader.AddArg("lz");
    argReader.AddArg("isHexahedral");
    argReader.AddArg("thetaxy");
    argReader.AddArg("thetayz");
    argReader.AddArg("thetaxz");
    argReader.AddArg("support_SH");
    argReader.AddArg("diago_lib");
    argReader.AddArg("support_Periodic_Boundary");
    argReader.AddArg("multi_parallel_strategies");
    argReader.AddArg("points_path");
    argReader.AddArg("v_path");
    argReader.AddArg("distribution_path");
    argReader.ReadArgs(fs);
    fs.close();
    return 0;
}