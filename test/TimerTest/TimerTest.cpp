//
// Created by Chiahwa Young on 2023/4/8.
//
#include <iostream>
#include <string>
#include "Timer/Timer.h"
#include "Mat_Demo/Mat_Demo.h"
#include "ArgumentReader/ArgumentReader.h"
#include <fstream>

#ifdef __MPI__

#include <mpi.h>

#endif

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
//TODO: Change back to serial version. This program doesn't need MPI implementation.
int main(int argc, char *argv[]) {

    int rank = 0, size = 0;
    std::ostream *os = &std::cout;
#ifdef __MPI__
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == 0)
        cout << "MPI enabled" << endl;
#endif
    HwaUtil::Timer::tick("HwaUtil::(root)", "TimerTest");
    HwaUtil::Timer::tick("HwaUtil::(root)", "main");

    /*
    std::cout << "\033[7m" << "Now calculating the first 15 Fibonacci numbers:" << "\033[0m" << std::endl;
    for (int i = 0; i < 10; i++)
        std::cout << fib(i) << std::endl;
    */
    // std::cout << "\033[7m" << "Now doing MatDemoTest" << "\033[0m" << std::endl;

    //获取参数
    ArgumentReader ar;
    if (rank == 0) {
        ar.AddArg("nrows");
        ar.AddArg("ncols");
        ar.AddArg("matrix_type");
        ar.AddArg("matrix_print");
        ar.AddArg("calculation");
        ar.AddArg("print_mpi_log");

        ifstream fs("input.txt");
        ar.ReadArgs(fs);
        fs.close();
#ifdef __MPI__
        auto print_mpi_log_str = ar.GetArgV("print_mpi_log");
        if (print_mpi_log_str == "1") {
            os = new std::ofstream ("processor_" + to_string(rank) + ".log");
        } else if (print_mpi_log_str == "0")
            os = nullptr;
        else {
            HwaUtil::Timer::tock("HwaUtil::(root)", "main");
            throw invalid_argument("Invalid print_mpi_log argument");
        }
#endif
    }
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
        HwaUtil::Timer::tock("HwaUtil::(root)", "main");
        throw invalid_argument("Invalid matrix type");
    }
    //行数列数
    int ncols = stoi(ar.GetArgV("ncols"));
    int nrows = stoi(ar.GetArgV("nrows"));

    Mat_Demo m1(nrows, ncols, matrixType);

    //根据给定的操作进行计算

    string cal = ar.GetArgV("calculation");
    double ans=0.0;
    if(cal=="max") {
        ans=m1.mmax();
    } else if(cal=="min") {
        ans=m1.mmin();
    } else {
        HwaUtil::Timer::tock("HwaUtil::(root)", "main");
        throw invalid_argument("Invalid calculation");
    }
    if(os) {
        (*os) << "calculation: " << cal << endl
        <<ans<<endl;
    }

    //按需打印矩阵内容
    string print = ar.GetArgV("matrix_print");
    if(os) {
        if (print == "1") {
            (*os) << "Matrix printed:" << endl;
            (*os) << m1 << endl;
        } else if (print == "0") {
            (*os) << "Matrix not printed" << endl;
        } else {
            HwaUtil::Timer::tock("HwaUtil::(root)", "main");
            throw invalid_argument("Invalid argument matrix_print");
        }
    }
    HwaUtil::Timer::tock("HwaUtil::(root)", "main");
    if(os)
        HwaUtil::Timer::print_time_usage(*os);
#ifdef __MPI__
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    HwaUtil::Timer::tock("HwaUtil::(root)","TimerTest");
    if(rank==0){
        cout<<"Time Elapsed: "<<((double)(HwaUtil::Timer::func_time("HwaUtil::(root)","TimerTest").count())
                                *std::chrono::nanoseconds::period::num/std::chrono::nanoseconds::period::den)<<endl;
    }
#ifdef __MPI__
    if(os) {
        ((fstream *) os)->close();
        delete os;
    }
    MPI_Finalize();
#endif
    return 0;
}