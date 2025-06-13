// common.h
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <chrono>

#include "Types.h"
// Stałe globalne

constexpr real epsilon = 1.19;
constexpr real sigma = 0.369;
constexpr real kB = 8.31e-3;
constexpr real m = 39.95;

constexpr int DEBUG = 1;
constexpr int Term = 500;
constexpr int Heat = 0;
constexpr int Cool = 0;
constexpr real T0 = 50;
constexpr real DT = 0.001;
constexpr int SAVE = 10;
constexpr const char* Fname = "argon_out.csv"; 

constexpr int N_Heating = 3;
constexpr real HeatingTime = 0.5f; // ps
constexpr real margin_coef = 0.5f;


// Pomiar czasu
using TimePoint = std::chrono::high_resolution_clock::time_point;
using Clock = std::chrono::high_resolution_clock;
using WallTimePoint = std::chrono::time_point<Clock>;

class TimerCounter {
private:
    TimePoint t_start;
    double total_time = 0.0;
    int count = 0;

public:
    void start() { t_start = Clock::now(); }
    void add() {
        auto t_end = Clock::now();
        std::chrono::duration<double> elapsed = t_end - t_start;
        total_time += elapsed.count();
        count++;
    }
    void report(const char* label, double total_reference = -1.0) const {
        double ms = total_time * 1e3;
        printf("%-25s: %8.3f ms", label, ms);
        if (count > 0) printf("   (%d razy, %.3f ms avg)", count, ms / count);
        if (total_reference > 0.0)
            printf("   [%.1f%% czasu całkowitego]", 100.0 * total_time / total_reference);
        printf("\n");
    }
    double getTotal() const { return total_time; }
    int getCount() const { return count; }
};

extern TimerCounter T_total, T_verlet, T_force;
extern TimerCounter T_rebuild, T_ghost, T_transfer_up, T_transfer_down, T_stats, T_io;
extern TimerCounter T_kernel, T_copy2gpu, T_copy2cpu;

