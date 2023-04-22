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

int main() {
    ArgumentReader ar;
    ar.AddArg("calculation");
    ar.AddArg("result_print");
    ar.AddArg("timer_print");
    ar.AddArg("input_type");
    ar.AddArg("data_path");

    ifstream fs("input.txt");
    ar.ReadArgs(fs);
    fs.close();
    //cout<< ar << endl;
    auto calculation_str = ar.GetArgV("calculation");
    if (calculation_str == "matmul") {
        Mat_Demo *m1, *m2;
        auto input_type_str = ar.GetArgV("input_type");
        if (input_type_str == "file") {
            auto data_path_str = ar.GetArgV("data_path");
            auto path1 = data_path_str.substr(0, data_path_str.find('|'));
            auto path2 = data_path_str.substr(data_path_str.find('|') + 1);
            cout<<std::filesystem::current_path()<<endl;
            fstream fs1(path1);
            std::cout<<fs1.fail()<<std::endl;
            fstream fs2(path2);
            m1 = new Mat_Demo(fs1);
            m2 = new Mat_Demo(fs2);
            fs1.close();
            fs2.close();
        } else {
            throw invalid_argument("Invalid input type");
        }
        auto m3_manual = (*m1) * (*m2);
        auto m3_blas = m1->blas_mult(*m2);

        auto result_print_str = ar.GetArgV("result_print");
        if(result_print_str == "1"){
            fstream fs_manual("mult_manual.txt", ios::out);
            fs_manual << m3_manual;
            fs_manual.close();
            fstream fs_blas("mult_blas.txt", ios::out);
            fs_blas << m3_blas;
            fs_blas.close();
        }

        delete m1, delete m2;
    } else if (calculation_str == "RSMdiago") {

    }

    auto timer_print_str = ar.GetArgV("timer_print");
    if(timer_print_str == "1")
        HwaUtil::Timer::print_time_usage();

    return 0;
}