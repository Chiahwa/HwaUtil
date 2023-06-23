
#ifdef __MPI__

#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>


// 函数声明
extern "C" {
// Scalapack函数接口
void pdsyevx_(char *jobz, char *range, char *uplo, int *n, double *a, int *ia, int *ja, int *desca, double *vl, double *vu,
         int *il, int *iu, double *abstol, int *m, double *w, double *z, int *iz, int *jz, int *descz, double *work,
         int *lwork, int *iwork, int *liwork, int *ifail, int *info);

// CBLACS函数接口
void Cblacs_pinfo(int *mypnum, int *nprocs);
void Cblacs_get(int context, int request, int *value);
void Cblacs_gridinit(int *context, char *order, int np_row, int np_col);
void Cblacs_gridinfo(int context, int *np_row, int *np_col, int *my_row, int *my_col);
void Cblacs_gridexit(int context);
void Cblacs_exit(int status);
void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt,
               int *lld, int *info);
int numroc_(int* n, int* nb, int* iproc, int* isrcproc, int* nprocs);
void pdgemr2d_(int* m, int* n, double* A, int* ia, int* ja, int* desca, double* B, int* ib, int* jb,
               int* descb, int* context);
}

int main(int argc, char **argv) {
    int n;                // 矩阵维度
    int my_rank;          // 当前进程的排名
    int num_procs;        // 进程总数
    int blacs_context;    // CBLACS上下文

    // 初始化MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // 初始化CBLACS
    Cblacs_pinfo(&my_rank, &num_procs);
    Cblacs_get(-1, 0, &blacs_context);

    if (my_rank == 0) {
        std::cout << "BLACS context: " << blacs_context << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // 按需设置矩阵维度n
    if (my_rank == 0) {
        n = 3 /* 设置矩阵维度 */;
        std::cout << "Matrix dimension: " << n << std::endl;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::cout << "Matrix dimension: " << n << std::endl;

    // 设置Scalapack相关参数
    char jobz = 'V';           // 计算特征向量
    char range = 'A';          // 计算所有特征值和特征向量
    char uplo = 'L';           // 输入矩阵为下三角
    int ia = 0;                // 局部矩阵的起始行索引
    int ja = 0;                // 局部矩阵的起始列索引
    int desca[9];              // 矩阵A的描述符
    int descz[9];              // 特征向量矩阵Z的描述符
    int info;                  // 返回的信息代码
    int lwork;                 // 工作数组长度
    int liwork;                // 工作整型数组长度
    double vl = 0.0;           // 范围下界
    double vu = 0.0;           // 范围上界
    int il = 0;                // 返回的特征值起始索引
    int iu = 0;                // 返回的特征值结束索引
    double abstol = 0.0;       // 绝对容差
    int m;                     // 返回的特征值个数
    std::vector<double> w(n);  // 特征值数组
    std::vector<double> z(n);  // 特征向量数组

    // 初始化矩阵A和描述符desca
    // std::vector<double> a(/* 初始化矩阵A */);
    std::vector<double> a = {
            1.0, 2.0, 3.0,
            2.0, 4.0, 5.0,
            3.0, 5.0, 6.0
    };

    // 初始化CBLACS网格
    int np_row = num_procs;       // 网格中的行数
    int np_col = 1;       // 网格中的列数
    int my_row;       // 当前进程所在的行索引
    int my_col;       // 当前进程所在的列索引
    //char order = 'Row'; // 列主序分布
    Cblacs_gridinit(&blacs_context, "Row", np_row, np_col);
    Cblacs_gridinfo(blacs_context, &np_row, &np_col, &my_row, &my_col);

    if (my_rank == 0) {
        std::cout << "Grid size: " << np_row << " x " << np_col << std::endl;

    }

    std::cout<<"My rank: "<< my_rank << ", myrow: " << my_row << ", mycol: " << my_col << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);


    // 设置描述符desca
    int info_1;
    descinit_(desca, &n, &n, &num_procs, &n, &ia, &ja, &blacs_context, &n, &info_1);

    // 设置描述符descz
    int info_2;
    descinit_(descz, &n, &n, &num_procs, &n, &ia, &ja, &blacs_context, &n, &info_2);

    // 计算工作数组长度
    // 设置初始值
    lwork = -1;
    liwork = -1;

    // 调用pdsyevx函数获取所需工作数组的最优长度
    pdsyevx_(&jobz, &range, &uplo, &n, a.data(), &ia, &ja, desca, &vl, &vu, &il, &iu, &abstol, &m, w.data(), z.data(),
             &ia, &ja, descz, nullptr, &lwork, nullptr, &liwork, nullptr, &info);

    // 获取所需的工作数组的最优长度
    lwork = static_cast<int>(w[0]);
    liwork = static_cast<int>(z[0]);

    // 根据返回的长度重新分配工作数组
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    // 调用Scalapack的pdsyevx计算特征值和特征向量
    pdsyevx_(&jobz, &range, &uplo, &n, a.data(), &ia, &ja, desca, &vl, &vu, &il, &iu, &abstol, &m, w.data(),
             z.data(), &ia, &ja, descz, nullptr, &lwork, nullptr, &liwork, nullptr, &info);

    if (info == 0) {
        // 打印特征值和特征向量
        if (my_rank == 0) {
            std::cout << "特征值：" << std::endl;
            for (int i = 0; i < m; ++i) {
                std::cout << w[i] << " ";
            }
            std::cout << std::endl;

            std::cout << "特征向量：" << std::endl;
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    std::cout << z[i * n + j] << " ";
                }
                std::cout << std::endl;
            }
        }
    } else {
        std::cerr << "pdsyevx调用失败，错误代码：" << info << std::endl;
    }


    // 释放资源
    Cblacs_gridexit(blacs_context);
    // Cblacs_exit(0);
    MPI_Finalize();

    return 0;
}


#else

#include <iostream>

int main() {
    std::cout << "Error: MPI not supported." << std::endl;
}

#endif
