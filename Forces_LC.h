#pragma once
#include "Types.h"
#include "Argon_LC.h"

// ============================================================================
// Deklaracje funkcji sił i energii dla wersji linked-cell (CPU i GPU)
// ============================================================================

// CPU: siły bez zasady Newtona
elj3 ForcesEnergies_LC(GridCells& Grid, parameters Par);

// GPU: różne warianty (implementacje w GPU.cu / GPU_LC.cu)
elj3 ForcesEnergies_GPU_LC(GridCells& Grid, parameters Par);        // wersja podstawowa
elj3 ForcesEnergies_GPU_LC_ar(GridCells& Grid, parameters Par);     // z zasadą akcji-reakcji
elj3 ForcesEnergies_GPU_LC_nl(GridCells& Grid, parameters Par);     // z listą sąsiedztwa

// Wrapper wybierający wersję FORCE_*_LC
// elj3 ForcesEnergiesWrapper_LC(GridCells& Grid, parameters Par);

#if defined(FORCE_CPU_LC) || defined(FORCE_GPU_LC_SIMPLE) || defined(FORCE_GPU_LC_AR) || defined(FORCE_GPU_LC_NL)

static inline elj3 ForcesEnergiesWrapper_LC(GridCells& Grid, parameters Par) {
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

#endif // definicje wrappera LC tylko dla FORCE_*_LC