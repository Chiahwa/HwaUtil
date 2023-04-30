//
// Created by Chia-hwa Young on 2023/4/8.
//

#ifndef HWAUTIL_TIMER_H
#define HWAUTIL_TIMER_H
#include <iostream>
#include <string>
#include <chrono>
#include <unordered_map>
#include <vector>
namespace HwaUtil {

    class Timer {
    private:
        // class modelling a record of a separate tick-to-tock period
        class TimerRecord{
        private:
            // number of recursion instances currently running as a CHILD of this instance,
            // i.e. the number of times tick() has been called without a corresponding tock() call.
            unsigned long long n_recursion=0;
        public:
            // time at which top-level tick() was called
            double start_time;
            // time at which top-level tock() was called
            double end_time;

            void init();
            bool is_running() const;
            void tick();
            void tock();
        };

        // class modelling the time usage of a function
        class FuncTimeInfo{
        private:
            //bool is_running;
            std::vector<TimerRecord> records;
            TimerRecord running_record;
        public:
            FuncTimeInfo():records(std::vector<TimerRecord>()){};
            void tick();
            void tock();
            double total_time();
            unsigned long long total_count();
        };

    private:
        // map to store the time usage of each function
        static std::unordered_map<std::string, FuncTimeInfo> func_time_info;
        // time at which the program started
        static double program_start_time;

    public:
        // call this function at the beginning of a function
        static void tick(const std::string &class_name, const std::string &function_name);
        // call this function at the end of a function
        static void tock(const std::string &class_name, const std::string &function_name);
        static std::ostream &print_time_usage(std::ostream &os = std::cout);

        static double func_time(const std::string &class_name, const std::string &function_name);
        static void init();
    };

} // HwaUtil

#endif //HWAUTIL_TIMER_H
