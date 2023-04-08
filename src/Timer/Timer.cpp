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

        os <<std::left<<std::setfill('-')
        << "\033[7m"
        << std::setw(25) << "|CLASS_NAME"
        << std::setw(20) << "|FUNC_NAME"
        << std::setw(16) << "|TIME(Sec)"
        << std::setw(10) << "|CALLS"
        << std::setw(13) << "|AVG(Sec)"
        << std::setw(13) << "|PER%"
        << "\033[0m"
        << std::endl;

        os<<std::right<<std::setfill(' ');
        auto prog_total_time = double (program_total_time.count())
                *std::chrono::nanoseconds::period::num/std::chrono::nanoseconds::period::den;
        os << "\033[32m"
           << std::setw(25) << "--"
           << std::setw(20) << "TOTAL"
           << std::setw(16) << prog_total_time
           << std::setw(10) << 1
           << std::setw(13) << prog_total_time
           << std::setw(12) << 100.0 << '%'
           << "\033[0m"
           << std::endl;
        for (auto &func_time:func_time_info) {
            auto func_total_time = double (func_time.second.total_time().count())
                    *std::chrono::nanoseconds::period::num/std::chrono::nanoseconds::period::den;
            auto total_count = func_time.second.total_count();
            auto avg_time = func_total_time / total_count;
            auto per_time = (double)func_total_time / prog_total_time;
            os << std::setw(25) << func_time.first.substr(0, func_time.first.find("|||"))
            << std::setw(20) << func_time.first.substr(func_time.first.find("|||") + 3)
            << std::setw(16) << func_total_time
            << std::setw(10) << total_count
            << std::setw(13) << avg_time
            << std::setw(12) << std::fixed<< std::setprecision(1)<<per_time*100 << '%'
            << std::endl;
        }

        os <<std::setfill('-')<<std::setw(97)<<"-"<<std::endl;
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