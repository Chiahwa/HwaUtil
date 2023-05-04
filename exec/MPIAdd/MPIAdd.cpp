//
// Created by Chiahwa Young on 2023/4/30.
//
#include <iostream>
#include <fstream>
#include <filesystem>
#include "Mat_Demo/Mat_Demo.h"
#include "ArgumentReader/ArgumentReader.h"
#include "Timer/Timer.h"

#ifdef __MPI__

#include <mpi.h>

#endif

using namespace std;
using namespace HwaUtil;

#ifdef __MPI__


Mat_Demo *global_m1 = nullptr;
Mat_Demo *global_m2 = nullptr;
Mat_Demo *receive_m1 = nullptr;
Mat_Demo *receive_m2 = nullptr;
Mat_Demo *blockresult = nullptr;
Mat_Demo *globalresult = nullptr;
int nrow_start_current, nrow_end_current, ncol_start_current, ncol_end_current, nrow_cur, ncol_cur;
int print_mpi_log, timer_print;
string output_to_file_str;
double alpha, beta;
int mpi_rank = 0, mpi_size = 0;

void transmit_matrix(Mat_Demo *global_m, Mat_Demo *&receive_m) {
    HwaUtil::Timer::tick("HwaUtil::(MPIAdd)", "transmit_matrix");
    double *data_send = nullptr;
    double *data_recv = nullptr;
    //transmit matrix
    if (mpi_rank == 0) { //send matrix & rank 0 receive matrix from itself
        cout << "Transmitting matrix..." << endl;
        int nrow = global_m->nr();
        int ncol = global_m->nc();
        int nrow_per_proc = nrow / mpi_size;
        int n_one_more = nrow % mpi_size;
        //int nrow_last_proc = nrow - nrow_per_proc * (mpi_size - 1);
        for (int i = 0; i < mpi_size; i++) {
            //prepare data to send
            int nrow_this_proc = nrow_per_proc + (i < n_one_more ? 1 : 0);
            int nrow_start = nrow_per_proc * i + (i < n_one_more ? i : n_one_more);
            int nrow_end = nrow_start + nrow_this_proc;

            int ncol_start = 0;
            int ncol_end = ncol;
            int ncol_this_proc = ncol_end - ncol_start;
            int n_send = nrow_this_proc * ncol_this_proc;
            data_send = new double[n_send];

            for (int j = nrow_start; j < nrow_end; j++) {
                for (int k = ncol_start; k < ncol_end; k++) {
                    data_send[(j - nrow_start) * ncol + k] = (*global_m)(j, k);
                }
            }
            if (i == 0) {
                nrow_start_current = nrow_start;
                nrow_end_current = nrow_end;
                ncol_start_current = ncol_start;
                ncol_end_current = ncol_end;
                data_recv = data_send;
            } else {
                MPI_Send(&nrow_start, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
                MPI_Send(&nrow_end, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
                MPI_Send(&ncol_start, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
                MPI_Send(&ncol_end, 1, MPI_INT, i, 4, MPI_COMM_WORLD);
                MPI_Send(data_send, n_send, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

                //MPI_Send(&print_mpi_log, 1, MPI_INT, i, 5, MPI_COMM_WORLD);
                //MPI_Send(&timer_print, 1, MPI_INT, i, 6, MPI_COMM_WORLD);
                delete[] data_send;
            }
        }
    } else { //other processes receive matrix
        MPI_Recv(&nrow_start_current, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&nrow_end_current, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&ncol_start_current, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&ncol_end_current, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //MPI_Recv(&print_mpi_log, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //MPI_Recv(&timer_print, 1, MPI_INT, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int n_recv = (nrow_end_current - nrow_start_current) * (ncol_end_current - ncol_start_current);
        data_recv = new double[n_recv];
        MPI_Recv(data_recv, n_recv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
    nrow_cur = nrow_end_current - nrow_start_current;
    ncol_cur = ncol_end_current - ncol_start_current;
    receive_m = new Mat_Demo(nrow_cur, ncol_cur, data_recv);
    delete [] data_recv;
    HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "transmit_matrix");
}

void read_data(int argc, char **argv) {
    HwaUtil::Timer::tick("HwaUtil::(MPIAdd)", "read_data");

    ArgumentReader ar;
    ar.AddArg("calculation");
    ar.AddArg("input_type");
    ar.AddArg("matrix_1");
    ar.AddArg("matrix_2");
    ar.AddArg("alpha");
    ar.AddArg("beta");

    ar.AddArg("timer_print");
    ar.AddArg("output_to_file");
    ar.AddArg("print_mpi_log");

    //path for argument file
    string input_file_path;
    if (argc != 2)
        getline(cin, input_file_path);
    else
        input_file_path = argv[1];
    ifstream fs(input_file_path);
    if (fs.fail())
        throw invalid_argument("Invalid input file path");
    cout << "Reading argument file: " << input_file_path << "..." << endl;
    ar.ReadArgs(fs);
    fs.close();
    string calculation_str = ar.GetArgV("calculation");
    string input_type_str = ar.GetArgV("input_type");
    string m1_path_str = ar.GetArgV("matrix_1");
    string m2_path_str = ar.GetArgV("matrix_2");
    string alpha_str=ar.GetArgV("alpha");
    string beta_str=ar.GetArgV("beta");

    string timer_print_str = ar.GetArgV("timer_print");
    output_to_file_str = ar.GetArgV("output_to_file");
    string print_mpi_log_str = ar.GetArgV("print_mpi_log");


    print_mpi_log = (print_mpi_log_str == "1") ? 1 : 0;
    timer_print = (timer_print_str == "1") ? 1 : 0;

    alpha = stod(alpha_str);
    beta = stod(beta_str);

    if (calculation_str == "matadd_mpi") {
        cout << "Calculation: matadd_mpi" << endl;
        //read matrix

        if (input_type_str == "file") {
            //auto data_path_str = ar.GetArgV("data_path");
            fstream fss(m1_path_str);
            if (fss.fail()) {
                HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "read_data");
                throw invalid_argument(m1_path_str + ": Invalid data file path");
            }
            cout << "Reading matrix from file: " << m1_path_str << "..." << endl;
            global_m1 = new Mat_Demo(fss);
            fss.close();
            fss.open(m2_path_str);
            if (fss.fail()) {
                HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "read_data");
                throw invalid_argument(m2_path_str + ": Invalid data file path");
            }
            cout << "Reading matrix from file: " << m2_path_str << "..." << endl;
            global_m2 = new Mat_Demo(fss);
            fss.close();
        } else {
            //cerr<<"Invalid input type"<<endl;
            HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "read_data");
            throw invalid_argument(input_type_str + ": Invalid input type");
            //MPI_Abort(MPI_COMM_WORLD, 1);
        }
    } else {
        //cerr << "Invalid calculation type"<<endl;
        HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "read_data");
        throw invalid_argument("Invalid calculation type");
        //MPI_Abort(MPI_COMM_WORLD, 1);
    }

    HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "read_data");
}

void print_read_log() {
    HwaUtil::Timer::tick("HwaUtil::(MPIAdd)", "print_read_log");
    if (mpi_rank == 0) {
        filesystem::remove_all("output");
        filesystem::create_directory("output");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    string log_file_name = "output/processor_" + to_string(mpi_rank) + ".log";
    ofstream log_file(log_file_name);
    log_file << "Processor ID: " << mpi_rank << endl;
    log_file << "Block Matrix ID: " << mpi_rank << endl;
    log_file << "Block Matrix Size: " << nrow_end_current - nrow_start_current << " x "
             << ncol_end_current - ncol_start_current << endl;
    log_file << "Start Position: (" << nrow_start_current << ", " << ncol_start_current << ")" << endl;
    log_file << "End Position: (" << nrow_end_current-1 << ", " << ncol_end_current-1 << ")" << endl;
    log_file << "Block Matrix Elements:" << endl;
    log_file << "Matrix_1--------------" << endl;
    for (int i=0;i<receive_m1->nr();++i){
        for (int j=0;j<receive_m1->nc();++j){
            log_file<<(*receive_m1)(i,j)<<", ";
        }
        log_file<<endl;
    }
    log_file << "Matrix_2--------------" << endl;
    for (int i=0;i<receive_m2->nr();++i){
        for (int j=0;j<receive_m2->nc();++j){
            log_file<<(*receive_m2)(i,j)<<", ";
        }
        log_file<<endl;
    }
    log_file.close();
    HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "print_read_log");
}

void bcast_args()
{
    HwaUtil::Timer::tick("HwaUtil::(MPIAdd)", "bcast_args");
    int rootproc=0;
    MPI_Bcast(&print_mpi_log, 1, MPI_INT, rootproc, MPI_COMM_WORLD);
    MPI_Bcast(&timer_print, 1, MPI_INT, rootproc, MPI_COMM_WORLD);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, rootproc, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, MPI_DOUBLE, rootproc, MPI_COMM_WORLD);
    HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "bcast_args");
}

void transmit_back() {
    HwaUtil::Timer::tick("HwaUtil::(MPIAdd)", "transmit_back");
    if (mpi_rank == 0) {
        cout << "Transmitting back result..." << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int rootproc = 0;
    const double *send_buf = blockresult->get_ptr(0,0);
    double *recv_buf;
    int send_size = nrow_cur * ncol_cur;
    int *recv_size=new int[mpi_size];
    int *recv_displs=new int[mpi_size];
    MPI_Gather(&send_size, 1, MPI_INT, recv_size, 1, MPI_INT, rootproc, MPI_COMM_WORLD);
    if(mpi_rank==rootproc){
        int total_size = 0;
        for (int i = 0; i < mpi_size; ++i) {
            recv_displs[i] = total_size;
            total_size += recv_size[i];
        }
        if(total_size!=global_m1->nr()*global_m1->nc()||total_size!=global_m2->nr()*global_m2->nc()){
            HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "transmit_back");
            throw runtime_error("Matrix size error");
        }
        recv_buf = new double[total_size];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(send_buf, send_size, MPI_DOUBLE, recv_buf, recv_size, recv_displs, MPI_DOUBLE, rootproc, MPI_COMM_WORLD);

    if(mpi_rank==rootproc)
        globalresult = new Mat_Demo(global_m1->nr(), global_m1->nc(), recv_buf);
    HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "transmit_back");
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (mpi_rank == 0) {
        cout << "MPI enabled" << endl;
        cout << "# of processors: " << mpi_size << endl;
    }

    HwaUtil::Timer::init();
    HwaUtil::Timer::tick("HwaUtil::(MPIAdd)", "main");
    // read argument file and matrix data.
    if (mpi_rank == 0)
    {
        read_data(argc, argv);
        if(global_m1->nr() != global_m2->nr() || global_m1->nc() != global_m2->nc())
        {
            throw invalid_argument("Matrix size not match");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // broadcast arguments
    bcast_args();
    MPI_Barrier(MPI_COMM_WORLD);

    // divide and distribute
    transmit_matrix(global_m1, receive_m1);
    transmit_matrix(global_m2, receive_m2);
    MPI_Barrier(MPI_COMM_WORLD);
    if (print_mpi_log) {
        print_read_log();
    }

    // calculate
    if (mpi_rank == 0) {
        cout << "Calculating..." << endl;
    }
    blockresult = &HwaUtil::mat_add(*receive_m1, *receive_m2, alpha, beta);
    MPI_Barrier(MPI_COMM_WORLD);

    // transmit back to processor 0
    if (mpi_rank == 0) {
        cout << "Transmitting result back..." << endl;
    }
    transmit_back();
    MPI_Barrier(MPI_COMM_WORLD);

    // print result
    if(mpi_rank==0){
        if(output_to_file_str != "0"){
            ofstream output_file(output_to_file_str);
            if(!output_file.is_open()){
                throw runtime_error("Cannot open output file");
            }
            cout << "Writing result to file: " <<output_to_file_str<<"..." << endl;
            output_file<<(*globalresult);
            output_file.close();
        }
    }

    HwaUtil::Timer::tock("HwaUtil::(MPIAdd)", "main");
    if (timer_print) {
        if (mpi_rank == 0)
            cout << "Time elapsed: " << (double) (HwaUtil::Timer::func_time("HwaUtil::(MPIAdd)", "main")) << "s"
                 << endl;
        if (mpi_rank == 0) {
            filesystem::create_directory("output");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        string log_file_name = "output/processor_" + to_string(mpi_rank) + "_timer.log";
        ofstream log_file(log_file_name);
        HwaUtil::Timer::print_time_usage(log_file);
        log_file.close();
    }

    //delete[] data_recv;
    delete global_m1;
    delete global_m2;
    delete receive_m1;
    delete receive_m2;
    delete globalresult;
    delete blockresult;
    MPI_Finalize();
    return 0;
}

#else
int main(int argc, char *argv[]) {
    throw std::runtime_error("MPI is not enabled at this compilation.");
}
#endif
