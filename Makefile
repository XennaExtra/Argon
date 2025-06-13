# ============================================================================
# Makefile – wersja modularna dla nowej struktury projektu
# Obsługuje CPU, GPU, LC, PBC. Każdy tryb to osobny target.
# Wersje GPU mają opisowe nazwy: naive, tiled, shared, newton
# LC i PBC kontrolowane przez flagi makr kompilacyjnych
# ============================================================================

CXX  = g++
NVCC = nvcc

CXXFLAGS = -O2 -Wall -std=c++11
CUFLAGS  = -O2 -arch=compute_35

COMMON_CPU = Common.cpp
COMMON_GPU = CommonData.cu

BINARIES = \
  ArgonCPU ArgonCPU_TILES ArgonCPU_PBC ArgonCPU_LC \
  ArgonGPU_Naive ArgonGPU_Naive_PBC \
  ArgonGPU_Tiled ArgonGPU_Tiled_PBC \
  ArgonGPU_Shared ArgonGPU_Shared_PBC \
  ArgonGPU_Newton ArgonGPU_Newton_PBC \
  ArgonGPU_LC

# ===============================
# Wersje CPU
# ===============================

ArgonCPU:
	$(CXX) $(CXXFLAGS) -c $(COMMON_CPU) -o Common.o
	$(CXX) $(CXXFLAGS) -DFORCE_CPU -c Argon.cpp -o Argon.o
	$(CXX) $(CXXFLAGS) -DFORCE_CPU -c Forces.cpp -o Forces.o
	$(CXX) $(CXXFLAGS) -DFORCE_CPU -c Main.cpp -o Main.o
	$(CXX) $(CXXFLAGS) -o $@ Common.o Argon.o Forces.o Main.o

ArgonCPU_PBC:
	$(CXX) $(CXXFLAGS) -c $(COMMON_CPU) -o Common.o
	$(CXX) $(CXXFLAGS) -DFORCE_CPU_PBC -DUSE_PBC -c Argon.cpp -o Argon_pbc.o
	$(CXX) $(CXXFLAGS) -DFORCE_CPU_PBC -DUSE_PBC -c Forces.cpp -o Forces_pbc.o
	$(CXX) $(CXXFLAGS) -DFORCE_CPU_PBC -DUSE_PBC -c Main.cpp -o Main_pbc.o
	$(CXX) $(CXXFLAGS) -o $@ Common.o Argon_pbc.o Forces_pbc.o Main_pbc.o

ArgonCPU_LC:
	$(CXX) $(CXXFLAGS) -c $(COMMON_CPU) -o Common.o
	#$(CXX) $(CXXFLAGS) -DFORCE_CPU_LC -c Argon.cpp -o Argon.o
	$(CXX) $(CXXFLAGS) -DFORCE_CPU_LC -c Argon_LC.cpp -o Argon_LC.o
	#$(CXX) $(CXXFLAGS) -DFORCE_CPU_LC -c Forces.cpp -o Forces.o
	$(CXX) $(CXXFLAGS) -DFORCE_CPU_LC -c Forces.cpp -o Forces_LC.o
	#$(CXX) $(CXXFLAGS) -DFORCE_CPU_LC -c Forces_LC.cpp -o Forces_LC.o
	$(CXX) $(CXXFLAGS) -DFORCE_CPU_LC -c Main_LC.cpp -o Main_LC.o
	#$(CXX) $(CXXFLAGS) -o $@ Common.o Argon.o Argon_LC.o Forces.o Forces_LC.o Main_LC.o
	$(CXX) $(CXXFLAGS) -o $@ Common.o Argon_LC.o Forces_LC.o Main_LC.o

# ===============================
# Wersje GPU (bez LC)
# ===============================

ArgonGPU:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -c $(COMMON_GPU) -o CommonData.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -c Argon.cu -o Argon_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU -c Forces.cu -o Forces_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU -c Main.cu -o Main_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -o $@ Argon_gpu.o Forces_gpu.o Main_gpu.o CommonData.o



ArgonGPU_PBC:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_PBC -DUSE_PBC -c Argon.cpp -o Argon_gpu_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_PBC -DUSE_PBC -c Forces.cpp -o Forces_gpu_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_PBC -DUSE_PBC -c GPU.cu -o GPU_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_PBC -DUSE_PBC -c Main.cpp -o Main_gpu_pbc.o
	$(NVCC) $(CUFLAGS) -o $@ Argon_gpu_pbc.o Forces_gpu_pbc.o GPU_pbc.o Main_gpu_pbc.o

ArgonGPU_Tiled:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_TILE -c Argon.cpp -o Argon_tile.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_TILE -c Forces.cpp -o Forces_tile.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_TILE -c GPU.cu -o GPU_tile.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_TILE -c Main.cpp -o Main_tile.o
	$(NVCC) $(CUFLAGS) -o $@ Argon_tile.o Forces_tile.o GPU_tile.o Main_tile.o

ArgonGPU_Tiled_PBC:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_TILE_PBC -DUSE_PBC -c Argon.cpp -o Argon_tile_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_TILE_PBC -DUSE_PBC -c Forces.cpp -o Forces_tile_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_TILE_PBC -DUSE_PBC -c GPU.cu -o GPU_tile_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_TILE_PBC -DUSE_PBC -c Main.cpp -o Main_tile_pbc.o
	$(NVCC) $(CUFLAGS) -o $@ Argon_tile_pbc.o Forces_tile_pbc.o GPU_tile_pbc.o Main_tile_pbc.o

ArgonGPU_Shared:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_SHARED -c Argon.cpp -o Argon_shared.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_SHARED -c Forces.cpp -o Forces_shared.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_SHARED -c GPU.cu -o GPU_shared.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_SHARED -c Main.cpp -o Main_shared.o
	$(NVCC) $(CUFLAGS) -o $@ Argon_shared.o Forces_shared.o GPU_shared.o Main_shared.o

ArgonGPU_Shared_PBC:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_SHARED_PBC -DUSE_PBC -c Argon.cpp -o Argon_shared_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_SHARED_PBC -DUSE_PBC -c Forces.cpp -o Forces_shared_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_SHARED_PBC -DUSE_PBC -c GPU.cu -o GPU_shared_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_SHARED_PBC -DUSE_PBC -c Main.cpp -o Main_shared_pbc.o
	$(NVCC) $(CUFLAGS) -o $@ Argon_shared_pbc.o Forces_shared_pbc.o GPU_shared_pbc.o Main_shared_pbc.o

ArgonGPU_Newton:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_AR -c Argon.cpp -o Argon_ar.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_AR -c Forces.cpp -o Forces_ar.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_AR -c GPU.cu -o GPU_ar.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_AR -c Main.cpp -o Main_ar.o
	$(NVCC) $(CUFLAGS) -o $@ Argon_ar.o Forces_ar.o GPU_ar.o Main_ar.o

ArgonGPU_Newton_PBC:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_AR_PBC -DUSE_PBC -c Argon.cpp -o Argon_ar_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_AR_PBC -DUSE_PBC -c Forces.cpp -o Forces_ar_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_AR_PBC -DUSE_PBC -c GPU.cu -o GPU_ar_pbc.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_AR_PBC -DUSE_PBC -c Main.cpp -o Main_ar_pbc.o
	$(NVCC) $(CUFLAGS) -o $@ Argon_ar_pbc.o Forces_ar_pbc.o GPU_ar_pbc.o Main_ar_pbc.o


# ===============================
# Wersje GPU (LC) – przyszłe
# ===============================


ArgonGPU_LC:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Argon.cpp -o Argon_lc_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Argon_LC.cpp -o Argon_LC_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Forces.cpp -o Forces_lc_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c GPU_LC.cu -o GPU_LC.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Main_LC.cpp -o Main_LC_gpu.o
	$(NVCC) $(CUFLAGS) -o $@ Argon_lc_gpu.o Argon_LC_gpu.o Forces_lc_gpu.o GPU_LC.o Main_LC_gpu.o


ArgonGPU_LC_Shared:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Argon.cpp -o Argon_lc_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Argon_LC.cpp -o Argon_LC_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Forces.cpp -o Forces_lc_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c GPU_LC.cu -o GPU_LC.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Main_LC.cpp -o Main_LC_gpu.o
	$(NVCC) $(CUFLAGS) -o $@ Argon_lc_gpu.o Argon_LC_gpu.o Forces_lc_gpu.o GPU_LC.o Main_LC_gpu.o

ArgonGPU_NL_LC:
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Argon.cpp -o Argon_lc_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Argon_LC.cpp -o Argon_LC_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Forces.cpp -o Forces_lc_gpu.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c GPU_LC.cu -o GPU_LC.o
	$(NVCC) $(CUFLAGS) -DUSE_GPU -DFORCE_GPU_LC_SIMPLE -c Main_LC.cpp -o Main_LC_gpu.o
	$(NVCC) $(CUFLAGS) -o $@ Argon_lc_gpu.o Argon_LC_gpu.o Forces_lc_gpu.o GPU_LC.o Main_LC_gpu.o

clean:
	rm -f *.o $(BINARIES)

.PHONY: all clean
