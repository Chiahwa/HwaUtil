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


Mat_Demo *global_mat = nullptr;
double *data_recv;
int nrow_start_current, nrow_end_current, ncol_start_current, ncol_end_current;
int print_mpi_log, timer_print;
int mpi_rank = 0, mpi_size = 0;

void transmit_matrix() {
    HwaUtil::Timer::tick("HwaUtil::(MatTransmit)", "transmit_matrix");
    //transmit matrix
    if (mpi_rank == 0) { //send matrix & rank 0 receive matrix from itself
        cout << "Transmitting matrix..." << endl;
        int nrow = global_mat->nr();
        int ncol = global_mat->nc();
        int nrow_per_proc = nrow / mpi_size;
        int nrow_last_proc = nrow - nrow_per_proc * (mpi_size - 1);
        for (int i = 0; i < mpi_size; i++) {
            //prepare data to send
            int nrow_this_proc = (i == mpi_size - 1) ? nrow_last_proc : nrow_per_proc;
            int nrow_start = i * nrow_per_proc;
            int nrow_end = nrow_start + nrow_this_proc;

            int ncol_start = 0;
            int ncol_end = ncol;
            int ncol_this_proc = ncol_end - ncol_start;
            int n_send = nrow_this_proc * ncol_this_proc;
            auto *data_send = new double[n_send];

            for (int j = nrow_start; j < nrow_end; j++) {
                for (int k = ncol_start; k < ncol_end; k++) {
                    data_send[(j - nrow_start) * ncol + k] = (*global_mat)(j, k);
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

                MPI_Send(&print_mpi_log, 1, MPI_INT, i, 5, MPI_COMM_WORLD);
                MPI_Send(&timer_print, 1, MPI_INT, i, 6, MPI_COMM_WORLD);
                delete[] data_send;
            }
        }
    } else { //other processes receive matrix
        MPI_Recv(&nrow_start_current, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&nrow_end_current, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&ncol_start_current, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&ncol_end_current, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&print_mpi_log, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&timer_print, 1, MPI_INT, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int nrow = nrow_end_current - nrow_start_current;
        int ncol = ncol_end_current - ncol_start_current;
        int n_recv = nrow * ncol;
        data_recv = new double[n_recv];
        MPI_Recv(data_recv, n_recv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    HwaUtil::Timer::tock("HwaUtil::(MatTransmit)", "transmit_matrix");
}

void read_data(int argc, char **argv) {
    HwaUtil::Timer::tick("HwaUtil::(MatTransmit)", "read_data");
    string print_mpi_log_str;
    string calculation_str;
    string input_type_str;
    string timer_print_str;
    string data_path_str;
    if (mpi_rank == 0) {
        ArgumentReader ar;
        ar.AddArg("calculation");
        ar.AddArg("timer_print");
        ar.AddArg("input_type");
        ar.AddArg("data_path");
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
        cout << "Reading argument file:" << input_file_path << "..." << endl;
        ar.ReadArgs(fs);
        fs.close();
        calculation_str = ar.GetArgV("calculation");
        input_type_str = ar.GetArgV("input_type");
        timer_print_str = ar.GetArgV("timer_print");
        print_mpi_log_str = ar.GetArgV("print_mpi_log");
        data_path_str = ar.GetArgV("data_path");

        print_mpi_log = (print_mpi_log_str == "1") ? 1 : 0;
        timer_print = (timer_print_str == "1") ? 1 : 0;

        if (calculation_str == "mpi_transmit") {
            cout << "Calculation: mpi_transmit" << endl;
            //read matrix

            if (input_type_str == "file") {
                //auto data_path_str = ar.GetArgV("data_path");
                fstream fss(data_path_str);
                if (fss.fail()) {
                    HwaUtil::Timer::tock("HwaUtil::(MatTransmit)", "read_data");
                    throw invalid_argument("Invalid data file path");
                }
                cout << "Reading matrix from file:" << data_path_str << "..." << endl;
                global_mat = new Mat_Demo(fss);
                fss.close();
            } else {
                //cerr<<"Invalid input type"<<endl;
                HwaUtil::Timer::tock("HwaUtil::(MatTransmit)", "read_data");
                throw invalid_argument("Invalid input type");
                //MPI_Abort(MPI_COMM_WORLD, 1);
            }
        } else {
            //cerr << "Invalid calculation type"<<endl;
            HwaUtil::Timer::tock("HwaUtil::(MatTransmit)", "read_data");
            throw invalid_argument("Invalid calculation type");
            //MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    HwaUtil::Timer::tock("HwaUtil::(MatTransmit)", "read_data");
}

void print_log() {
    HwaUtil::Timer::tick("HwaUtil::(MatTransmit)", "print_log");
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
    log_file << "End Position: (" << nrow_end_current << ", " << ncol_end_current << ")" << endl;
    log_file << "Block Matrix Elements:" << endl;
    for (int i = 0; i < nrow_end_current - nrow_start_current; i++) {
        for (int j = 0; j < ncol_end_current - ncol_start_current; j++) {
            log_file << data_recv[i * (ncol_end_current - ncol_start_current) + j] << ", ";
        }
        log_file << endl;
    }
    log_file.close();
    HwaUtil::Timer::tock("HwaUtil::(MatTransmit)", "print_log");
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
    HwaUtil::Timer::tick("HwaUtil::(MatTransmit)", "main");
    // read argument file and matrix data.
    read_data(argc, argv);
    MPI_Barrier(MPI_COMM_WORLD);

    transmit_matrix();
    if (print_mpi_log) {
        print_log();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    HwaUtil::Timer::tock("HwaUtil::(MatTransmit)", "main");
    if (timer_print) {
        if (mpi_rank == 0)
            cout << "Time elapsed: " << (double) (HwaUtil::Timer::func_time("HwaUtil::(MatTransmit)", "main")) << "s"
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

    delete[] data_recv;
    delete global_mat;
    MPI_Finalize();
    return 0;
}

#else
int main(int argc, char *argv[]) {
    throw std::runtime_error("MPI is not enabled at this compilation.");
}
#endif
