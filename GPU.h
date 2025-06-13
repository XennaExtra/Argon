#pragma once
#include "Types.h"
#include "Common.h"

#include <chrono>

#ifndef TILE
#define TILE 32
#endif

#define CALL_CUDA(expr) \
    do { \
        cudaError_t status = (expr); \
        if (status != cudaSuccess) { \
            printf("%s\n", cudaGetErrorString(status)); \
            exit(1); \
        } \
    } while(0)

extern AtomDataGPU ArgonGPU;
extern AtomDataGPU ArgonCPU;

elj3 ForcesEnergies_GPU(int Natoms, AtomData& Atoms, parameters Par);
void InitAllDataGPU(AtomDataGPU *GA, AtomDataGPU *CA, int SIZ);

using TimePoint = std::chrono::high_resolution_clock::time_point;
using Clock = std::chrono::high_resolution_clock;
using WallTimePoint = std::chrono::time_point<Clock>;

__global__ void ForcesEnergiesKernel_naive(int Natoms, AtomDataGPU atoms, real BoxSize, real Rc, real sig6, real sig12);
__global__ void ForcesEnergiesKernel_tiles(int Natoms, AtomDataGPU atoms, real BoxSize, real Rc, real sig6, real sig12, real *gbuf);
__global__ void ForcesEnergiesKernel_shared(int Natoms, AtomDataGPU atoms, real BoxSize, real Rc, real sig6, real sig12, real *host_buf, real *guest_buf);
__global__ void ForcesEnergiesKernel_newton(int Natoms, AtomDataGPU atoms, real BoxSize, real Rc, real sig6, real sig12, real *host_buf, real *guest_buf);
__global__ void reduce_buffer(int Natoms, AtomDataGPU atoms, real *gbuf, int buf_size);
