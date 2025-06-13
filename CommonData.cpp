#include "Types.h"
#include "Common.h"
#include "Argon.h"
#ifdef USE_GPU
#include "GPU.h"
#endif

// Główne dane globalne

AtomData Argon;         // struktura hostowa (CPU)
parameters Pars;        // parametry symulacji

AtomDataGPU ArgonGPU;   // struktura danych na GPU
AtomDataGPU ArgonCPU;   // kopia na CPU dla transferu (host mirror)

void msg(const char* name){
    printf("%s done\n",name);
};