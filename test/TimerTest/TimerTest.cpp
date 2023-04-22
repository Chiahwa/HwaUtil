//
// Created by Chiahwa Young on 2023/4/8.
//
#include <iostream>
#include <string>
#include "Timer/Timer.h"
#include "Mat_Demo/Mat_Demo.h"
#include "ArgumentReader/ArgumentReader.h"
#include <fstream>

using namespace std;
using namespace HwaUtil;

long long fib(int n) {
    HwaUtil::Timer::tick("HwaUtil::(root)", "fib");
    long long ans;
    if (n == 0) ans = 0;
    else if (n == 1) ans = 1;
    else ans = fib(n - 1) + fib(n - 2);
    HwaUtil::Timer::tock("HwaUtil::(root)", "fib");
    return ans;
}

int main() {

    HwaUtil::Timer::tick("HwaUtil::(root)", "main");

    std::cout << "\033[7m" << "Now calculating the first 15 Fibonacci numbers:" << "\033[0m" << std::endl;
    for (int i = 0; i < 10; i++)
        std::cout << fib(i) << std::endl;

    std::cout << "\033[7m" << "Now doing MatDemoTest" << "\033[0m" << std::endl;

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
    } else if (typestr == "user") {
        matrixType = Mat_Demo::MatrixType::User;
    } else {
        cout << "Invalid matrix type" << endl;
        return 1;
    }
    //行数列数
    int ncols = stoi(ar.GetArgV("ncols"));
    int nrows = stoi(ar.GetArgV("nrows"));

    Mat_Demo m1(nrows, ncols, matrixType, cin);

    //根据给定的操作进行计算
    cout << "calculation: ";
    string cal = ar.GetArgV("calculation");
    cout << cal << endl;
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

    HwaUtil::Timer::tock("HwaUtil::(root)", "main");
    HwaUtil::Timer::print_time_usage(std::cout);
    return 0;
}