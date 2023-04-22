//
// Created by Chiahwa Young on 2023/4/19.
//
#include <iostream>
#include <fstream>
#include "Mat_Demo/Mat_Demo.h"
#include "ArgumentReader/ArgumentReader.h"

using namespace std;
using namespace HwaUtil;
int main()
{
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
    auto calculation_str=ar.GetArgV("calculation");
    if(calculation_str=="matmul"){
        auto input_type_str=ar.GetArgV("input_type");
        if(input_type_str=="file"){
            auto data_path_str=ar.GetArgV("data_path");
            auto path1=data_path_str.substr(0, data_path_str.find('|'));
            auto path2=data_path_str.substr(data_path_str.find('|')+1);

        }else{
            throw invalid_argument("Invalid input type");
        }
    }
    else if(calculation_str=="RSMdiago"){

    }
    return 0;
}