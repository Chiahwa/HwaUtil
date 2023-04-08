//
// Created by Chia-hwa Young on 2023/4/8.
//

#include "Timer.h"
#include <iostream>
#include <string>
#include <chrono>
#include <unordered_map>
#include <iomanip>

namespace HwaUtil {

    // initialize static variables
    auto Timer::program_start_time = std::chrono::high_resolution_clock::now();
    std::unordered_map<std::string, Timer::FuncTimeInfo> Timer::func_time_info;

    void Timer::tick(const std::string &class_name, const std::string &function_name) {
        std::string func_name = class_name + "|||" + function_name;
        func_time_info[func_name].tick(); //TODO: need I check if func_name exists?
    }

    std::ostream &Timer::print_time_usage(std::ostream &os) {
        auto program_end_time = std::chrono::high_resolution_clock::now();
        auto program_total_time = program_end_time - program_start_time;

        os << "Program total time: " << std::chrono::duration_cast<std::chrono::microseconds>(program_total_time).count() << " microseconds" << std::endl;
        const int w_class_name = 28, w_func_name = 28, w_time = 16, w_calls = 10, w_avg = 14, w_per = 13;

        os <<std::left<<std::setfill('-')
        << "\033[7m"
        << std::setw(w_class_name) << "|CLASS_NAME"
        << std::setw(w_func_name) << "|FUNC_NAME"
        << std::setw(w_time) << "|TIME(Sec)"
        << std::setw(w_calls) << "|CALLS"
        << std::setw(w_avg) << "|AVG(Sec)"
        << std::setw(w_per) << "|PER%"
        << "\033[0m"
        << std::endl;

        os<<std::right<<std::setfill(' ');
        auto prog_total_time = double (program_total_time.count())
                *std::chrono::nanoseconds::period::num/std::chrono::nanoseconds::period::den;
        os << "\033[32m"
           << std::setw(w_class_name) << "--"
           << std::setw(w_func_name) << "TOTAL"
           << std::setw(w_time) << prog_total_time
           << std::setw(w_calls) << 1
           << std::setw(w_avg) << prog_total_time
           << std::setw(w_per-1) << 100.0 << '%'
           << "\033[0m"
           << std::endl;
        for (auto &func_time:func_time_info) {
            auto func_total_time = double (func_time.second.total_time().count())
                    *std::chrono::nanoseconds::period::num/std::chrono::nanoseconds::period::den;
            auto total_count = func_time.second.total_count();
            auto avg_time = func_total_time / total_count;
            auto per_time = (double)func_total_time / prog_total_time;
            os << std::setw(w_class_name) << func_time.first.substr(0, func_time.first.find("|||"))
            << std::setw(w_func_name) << func_time.first.substr(func_time.first.find("|||") + 3)
            << std::setw(w_time) << std::setprecision(6)<<func_total_time
            << std::setw(w_calls) << total_count
            << std::setw(w_avg) << avg_time
            << std::setw(w_per-1) << std::fixed<< std::setprecision(1)<<per_time*100 << '%'
            <<std::resetiosflags(std::ios::fixed)
            << std::endl;
        }

        os <<std::setfill('-')<<std::setw(w_calls+w_avg+w_time+w_func_name+w_class_name+w_per)<<"-"<<std::endl;
        return os;
    }

    void Timer::tock(const std::string &class_name, const std::string &function_name) {
        std::string func_name = class_name + "|||" + function_name;
        func_time_info[func_name].tock();
    }

    void Timer::FuncTimeInfo::tick() {
        if (!running_record.is_running()) {
            running_record.init();
        }
        running_record.tick();
    }
    void Timer::FuncTimeInfo::tock() {
        running_record.tock();
        if(!running_record.is_running()){
            records.push_back(running_record);
        }
    }
    std::chrono::nanoseconds Timer::FuncTimeInfo::total_time() {
        std::chrono::nanoseconds totaltime(0);
        for (auto &record:records) {
            totaltime += record.end_time - record.start_time;
        }
        return totaltime;
    }

    unsigned long long Timer::FuncTimeInfo::total_count() {
        return records.size();
    }

    bool Timer::TimerRecord::is_running() const {
        return n_recursion > 0;
    }

    void Timer::TimerRecord::init() {
        n_recursion = 0;
    }

    void Timer::TimerRecord::tick() {
        if (n_recursion == 0) {
            start_time = std::chrono::high_resolution_clock::now();
        }
        n_recursion++;

    }

    void Timer::TimerRecord::tock() {
        n_recursion--;
        if (n_recursion == 0) {
            end_time = std::chrono::high_resolution_clock::now();
        }
    }
} // HwaUtil