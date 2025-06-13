#pragma once
#include "Types.h"
#include "Common.h"
#include "Argon.h"

/*
#ifdef USE_GPU // USE_GPU

// ===================================
// Sekcja GPU-aware
// ===================================
#include "GPU.h"  // jawny warunek – tylko dla USE_GPU

elj3 ForcesEnergies_GPU(int Natoms, AtomData& Atoms, AtomDataGPU& ArgonCPU, AtomDataGPU& ArgonGPU, parameters Par);

elj3 ForcesEnergies_LC_GPU(GridCells& Grid, parameters Par);

//inline elj3 ForcesEnergiesWrapper(int Natoms, AtomData& Atoms, parameters Par) {
//    return ForcesEnergies_GPU(Natoms, Atoms, ArgonCPU, ArgonGPU, Par);
}

#elif // NO_GPU
*/ 

// ===================================
// CPU-only sekcja (brak dostępu do GPU.h)
// ===================================

#if !defined(FORCE_CPU_LC) && !defined(FORCE_GPU_LC_SIMPLE) && !defined(FORCE_GPU_LC_AR) && !defined(FORCE_GPU_LC_NL)
//#if !defined(FORCE_GPU_LC_SIMPLE) && !defined(FORCE_GPU_LC_AR) && !defined(FORCE_GPU_LC_NL)

elj3 ForcesEnergies(int Natoms, AtomData& Atoms, real Rc);
elj3 ForcesEnergies_pbc(int Natoms, AtomData& Atoms, real Rc);
elj3 ForcesEnergies_LC(GridCells& Grid, parameters Par);

inline elj3 ForcesEnergiesWrapper(int Natoms, AtomData& Atoms, parameters Par) {
#ifdef FORCE_CPU
    return ForcesEnergies(Natoms, Atoms, Par.Rc);
#elif defined(FORCE_CPU_PBC)
    return ForcesEnergies_pbc(Natoms, Atoms, Par.Rc);
#elif defined(FORCE_CPU_LC)
    return ForcesEnergies_LC(Grid, Par);
#elif defined(USE_GPU)
    return ForcesEnergies_GPU(Natoms, Atoms, Par.Rc);
#else
    #error "Brak obsługi dla zadanej konfiguracji "
#endif
}

#endif


// ============================================================================
// Deklaracje funkcji sił i energii dla wersji CPU (prostej)
// ============================================================================

// bez PBC
//elj3 ForcesEnergies(int Natoms, AtomData& Atoms, real Rc);        

// z PBC
//elj3 ForcesEnergies_pbc(int Natoms, AtomData& Atoms, real Rc);   

// z LC
//elj3 ForcesEnergies_LC(GridCells& Grid, parameters Par);

/*
// Wrapper: wybiera właściwą wersję w zależności od makr kompilacji
//elj3 ForcesEnergiesWrapper(int Natoms, AtomData& Atoms, parameters Par);

static inline elj3 ForcesEnergiesWrapper(int Natoms, AtomData& Atoms, parameters Par) {
#ifdef FORCE_CPU
    return ForcesEnergies(Natoms, Atoms, Par.Rc);
#elif defined(FORCE_CPU_PBC)
    return ForcesEnergies_pbc(Natoms, Atoms, Par.Rc);
#elif defined(FORCE_GPU)
    return ForcesEnergies_GPU(Natoms, Atoms, Par.Rc);
#elif defined(FORCE_GPU_PBC)
    return ForcesEnergies_GPU_pbc(Natoms, Atoms, Par.Rc);
#else
    // fallback domyślny
    return ForcesEnergies(Natoms, Atoms, Par.Rc);
#endif

}

#if defined(FORCE_CPU_LC) || defined(FORCE_GPU_LC_SIMPLE) || \
    defined(FORCE_GPU_LC_AR) || defined(FORCE_GPU_LC_NL)

elj3 ForcesEnergiesWrapper_LC(GridCells& Grid, parameters Par) {
#ifdef FORCE_CPU_LC
    return ForcesEnergies_LC(Grid, Par);
#elif defined(FORCE_GPU_LC_SIMPLE)
    return ForcesEnergies_GPU_LC(Grid, Par);
#elif defined(FORCE_GPU_LC_AR)
    return ForcesEnergies_GPU_LC_ar(Grid, Par);
#elif defined(FORCE_GPU_LC_NL)
    return ForcesEnergies_GPU_LC_nl(Grid, Par);
#else
    #error "Brak zdefiniowanej wersji FORCE_* dla modelu LC"
#endif
}
#endif

*/
