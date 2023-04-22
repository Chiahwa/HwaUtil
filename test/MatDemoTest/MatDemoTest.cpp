//
// Created by Chiahwa Young on 2023/3/22.
//

#include <iostream>
#include <fstream>
#include <Mat_Demo/Mat_Demo.h>
#include <ArgumentReader/ArgumentReader.h>
#include <string>

using namespace std;
using namespace HwaUtil;

int main() {
    //获取参数
    ArgumentReader ar;
    ar.AddArg("nrows");
    ar.AddArg("ncols");
    ar.AddArg("matrix_type");
    ar.AddArg("matrix_print");
    ar.AddArg("calculation");
    ifstream fs("input.txt");
    ar.ReadArgs(fs);
    fs.close();

    //获取矩阵初始化方式
    auto typestr = ar.GetArgV("matrix_type");
    Mat_Demo::MatrixType matrixType;
    if (typestr == "zero") {
        matrixType = Mat_Demo::MatrixType::Zero;
    } else if (typestr == "identity") {
        matrixType = Mat_Demo::MatrixType::Identity;
    } else if (typestr == "random") {
        matrixType = Mat_Demo::MatrixType::Random;
    } else {
        cout << "Invalid matrix type" << endl;
        return 1;
    }
    //行数列数
    int ncols= stoi(ar.GetArgV("ncols"));
    int nrows= stoi(ar.GetArgV("nrows"));

    Mat_Demo m1(nrows, ncols, matrixType);

    //根据给定的操作进行计算
    cout<<"calculation: ";
    string  cal = ar.GetArgV("calculation");
    cout<<cal<<endl;
    if (cal == "max") {
        cout << m1.mmax() << endl;
    } else if (cal == "min") {
        cout << m1.mmin() << endl;
    } else {
        cout << "Invalid calculation" << endl;
        return 1;
    }

    //按需打印矩阵内容
    string print = ar.GetArgV("matrix_print");
    if (print == "1") {
        cout << "Matrix printed:" << endl;
        cout << m1 << endl;
    } else if (print == "0") {
        cout << "Matrix not printed" << endl;
    } else {
        cout << "Invalid argument matrix_print" << endl;
        return 1;
    }
    return 0;
}