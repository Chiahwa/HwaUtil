//
// Created by Chiahwa Young on 2023/4/8.
//
#include <iostream>
#include <string>
#include "Timer/Timer.h"

long long fib(int n) {
    HwaUtil::Timer::tick("HwaUtil::main", "fib");
    long long ans;
    if (n == 0) ans = 0;
    else if (n == 1) ans = 1;
    else ans = fib(n - 1) + fib(n - 2);
    HwaUtil::Timer::tock("HwaUtil::main", "fib");
    return ans;
}

int main() {

    HwaUtil::Timer::tick("HwaUtil::main", "main");

    std::cout << "Now calculating the first 10 Fibonacci numbers:" << std::endl;
    for (int i = 0; i < 10; i++)
        std::cout << fib(i) << std::endl;

    HwaUtil::Timer::tock("HwaUtil::main", "main");
    HwaUtil::Timer::print_time_usage(std::cout);
    return 0;
}