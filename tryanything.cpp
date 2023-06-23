
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

// 由本地矩阵索引计算对应的全局索引
int local_to_global(int local_index, int block_size, int process_coord, int process_num) {
    return (local_index / block_size * process_num + process_coord ) * block_size + local_index % block_size;
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
        n = 4 /* 设置矩阵维度 */;
        std::cout << "Matrix dimension: " << n << std::endl;
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 设置Scalapack相关参数
    char jobz = 'V';           // 计算特征向量
    char range = 'A';          // 计算所有特征值和特征向量
    char uplo = 'L';           // 输入矩阵为下三角
    int ia = 1;                // 局部矩阵的起始行索引
    int ja = 1;                // 局部矩阵的起始列索引
    int mb = 2;                // 局部矩阵的行数
    int nb = 2;                // 局部矩阵的列数
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
    std::vector<double> a = {
            1.0, 2.0, 3.0, 4.0,
            2.0, 4.0, 5.0, 7.0,
            3.0, 5.0, 6.0, 8.0,
            4.0, 7.0, 8.0, 9.0
    };

    // 初始化CBLACS网格
    int np_row = 2;       // 网格中的行数
    int np_col = 2;       // 网格中的列数
    int my_row;       // 当前进程所在的行索引
    int my_col;       // 当前进程所在的列索引
    char order[9] = "Row"; // 行主序分布
    Cblacs_gridinit(&blacs_context, order, np_row, np_col);
    Cblacs_gridinfo(blacs_context, &np_row, &np_col, &my_row, &my_col);
    if (my_rank == 0) {
        std::cout << "Grid size: " << np_row << " x " << np_col << std::endl;
    }
    std::cout<<"My rank: "<< my_rank << ", myrow: " << my_row << ", mycol: " << my_col << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // 计算当前进程的子矩阵的行数和列数
    int m_loc = numroc_(&n, &mb, &my_row, &ia, &np_row);
    int n_loc = numroc_(&n, &nb, &my_col, &ja, &np_col);
    std::cout<<"My rank: "<< my_rank << ", m_loc: " << m_loc << ", n_loc: " << n_loc << std::endl;

    // 为局部矩阵分配内存
    std::vector<double> a_loc(m_loc * n_loc);
    std::vector<double> z_loc(m_loc * n_loc); //TODO: N * N?

    // 将矩阵A分发到各个进程
    for (int i = 0; i < m_loc; i++) {
        for (int j = 0; j < n_loc; j++) {
            int glob_i = local_to_global(i, mb, my_row, np_row);
            int glob_j = local_to_global(j, nb, my_col, np_col);

            a_loc[i + j * m_loc] = a[glob_i + glob_j * n];
        }
    }

    // 设置描述符desca
    int info_1;
    descinit_(desca, &n, &n, &mb, &nb, &ia, &ja, &blacs_context, &m_loc, &info_1);
    std::cout<< "desca: " << std::endl;
    for (int i : desca) {
        std::cout<< i << " ";
    }
    std::cout<< std::endl;

    // 设置描述符descz
    int info_2;
    descinit_(descz, &n, &n, &mb, &nb, &ia, &ja, &blacs_context, &m_loc, &info_2);
    std::cout<< "descz: " << std::endl;
    for (int i : descz) {
        std::cout<< i << " ";
    }
    std::cout<< std::endl;

    // 计算工作数组长度
    // 设置初始值
    lwork = -1;
    liwork = -1;

    // 调用pdsyevx函数获取所需工作数组的最优长度
    pdsyevx_(&jobz, &range, &uplo, &n, a_loc.data(), &ia, &ja, desca, &vl, &vu, &il, &iu, &abstol, &m, w.data(), z_loc.data(),
             &ia, &ja, descz, nullptr, &lwork, nullptr, &liwork, nullptr, &info);

    MPI_Barrier(MPI_COMM_WORLD);
    // 获取所需的工作数组的最优长度
    lwork = static_cast<int>(w[0]);
    liwork = static_cast<int>(z[0]);

    if (my_rank == 0) {
        std::cout << "lwork: " << lwork << std::endl;
        std::cout << "liwork: " << liwork << std::endl;
    }

    // 根据返回的长度重新分配工作数组
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    // 调用Scalapack的pdsyevx计算特征值和特征向量
    pdsyevx_(&jobz, &range, &uplo, &n, a_loc.data(), &ia, &ja, desca, &vl, &vu, &il, &iu, &abstol, &m, w.data(),
             z_loc.data(), &ia, &ja, descz, nullptr, &lwork, nullptr, &liwork, nullptr, &info);

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
