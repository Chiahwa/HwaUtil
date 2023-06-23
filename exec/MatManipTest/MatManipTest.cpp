//
// Created by Chiahwa Young on 2023/4/19.
//
#include <iostream>
#include <fstream>
#include <filesystem>
#include "Mat_Demo/Mat_Demo.h"
#include "ArgumentReader/ArgumentReader.h"
#include "Timer/Timer.h"

using namespace std;
using namespace HwaUtil;

int main(int argc, char *argv[]) {
    ArgumentReader ar;
    ar.AddArg("calculation");
    ar.AddArg("result_print");
    ar.AddArg("timer_print");
    ar.AddArg("input_type");
    ar.AddArg("data_path");

    //path for argument file
    string input_file_path;
    if (argc != 2)
        getline(cin, input_file_path);
    else input_file_path = argv[1];

    ifstream fs(input_file_path);
    if (fs.fail())
        throw invalid_argument("Invalid input file path");
    cout << "Reading argument file:" << input_file_path << "..." << endl;
    ar.ReadArgs(fs);
    fs.close();

    //get calculation type
    auto calculation_str = ar.GetArgV("calculation");
    auto input_type_str = ar.GetArgV("input_type");
    auto result_print_str = ar.GetArgV("result_print");
    auto timer_print_str = ar.GetArgV("timer_print");
    if (calculation_str == "matmul") {
        cout << "Calculation: matmul" << endl;
        //read two matrices
        Mat_Demo *m1, *m2;
        if (input_type_str == "file") {
            auto data_path_str = ar.GetArgV("data_path");
            auto path1 = data_path_str.substr(0, data_path_str.find('|'));
            auto path2 = data_path_str.substr(data_path_str.find('|') + 1);

            fstream fs1(path1);
            fstream fs2(path2);
            m1 = new Mat_Demo(fs1);
            m2 = new Mat_Demo(fs2);
            fs1.close();
            fs2.close();
        } else {
            throw invalid_argument("Invalid input type");
        }

        // do matrix multiplication
        cout << "Performing BLAS matrix multiplication..." << endl;
        auto m3_blas = m1->blas_mult(*m2);
        cout << "Performing manual matrix multiplication..." << endl;
        auto m3_manual = (*m1) * (*m2);

        //print result
        if (result_print_str == "1") {
            cout << R"(Writing result to files "mult_manual.txt" and "mult_manual.txt"...)" << endl;
            fstream fs_manual("mult_manual.txt", ios::out);
            fs_manual << m3_manual;
            fs_manual.close();
            fstream fs_blas("mult_blas.txt", ios::out);
            fs_blas << m3_blas;
            fs_blas.close();
        }

        delete m1, delete m2;
    } else if (calculation_str == "rsmdiago") {
        cout << "Calculation: RSMdiago" << endl;
        //read matrix
        Mat_Demo *m1;
        if (input_type_str == "file") {
            auto data_path_str = ar.GetArgV("data_path");
            fstream fs(data_path_str);
            m1 = new Mat_Demo(fs);
            fs.close();
        } else {
            throw invalid_argument("Invalid input type");
        }

        // do RSMdiago
        cout << "Performing RSMdiago..." << endl;
        int n = m1->nc();
        double *eigval = new double[n];
        double *eigvec = new double[n * n];
        if (m1->lapack_eig(eigval, eigvec) == 0) {
            if (result_print_str == "1") {
                cout << R"(Writing eigenvalues and eigenvectors to files "eig.txt" and "eigvec.txt"...)" << endl;
                fstream fs_eigval("eig.txt", ios::out);
                fstream fs_eigvec("eigvec.txt", ios::out);
                for (int i = 0; i < n; i++)
                    fs_eigval << eigval[i] << endl;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++)
                        fs_eigvec << eigvec[i * n + j] << " ";
                    fs_eigvec << endl;
                }
            }
        } else
            cout << "RSMdiago didn't finish" << endl;

    } else {
        throw invalid_argument("Invalid calculation");
    }
    if (timer_print_str == "1")
        HwaUtil::Timer::print_time_usage();

    return 0;
}