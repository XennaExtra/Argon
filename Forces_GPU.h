#pragma once
#include "Types.h"
#include "Common.h"
#include "Argon.h"
#include "GPU.h"  // jawny warunek – tylko dla USE_GPU


elj3 ForcesEnergies_GPU(int Natoms, AtomData& Atoms, AtomDataGPU& ArgonCPU, AtomDataGPU& ArgonGPU, parameters Par);

elj3 ForcesEnergies_LC_GPU(GridCells& Grid, parameters Par);



