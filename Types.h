#pragma once

#include <chrono>

using TimePoint = std::chrono::high_resolution_clock::time_point;
using Clock = std::chrono::high_resolution_clock;
using WallTimePoint = std::chrono::time_point<Clock>;

using real = float;

//#include "common.h"

struct AtomData {
    real *X;   // aktualne położenia atomów (X)
    real *Y;   // aktualne położenia atomów (Y)
    real *Z;   // aktualne położenia atomów (Z)

    real *Vx;  // prędkości atomów (X)
    real *Vy;  // prędkości atomów (Y)
    real *Vz;  // prędkości atomów (Z)

    real *Fx;  // siły działające na atomy (X)
    real *Fy;  // siły działające na atomy (Y)
    real *Fz;  // siły działające na atomy (Z)

    real *X0;  // referencyjne położenia początkowe – do śledzenia przemieszczeń
    real *Y0;
    real *Z0;

    bool* track;   // czy atom podlega śledzeniu (np. wykluczenie atomów brzegowych)
    int Natoms;    // liczba atomów w strukturze
};

struct AtomDataGPU {
    real *X;  // położenia atomów w kierunku X (GPU)
    real *Y;  // położenia atomów w kierunku Y (GPU)
    real *Z;  // położenia atomów w kierunku Z (GPU)

    real *Fx; // siły działające na atomy w kierunku X (GPU)
    real *Fy; // siły działające na atomy w kierunku Y (GPU)
    real *Fz; // siły działające na atomy w kierunku Z (GPU)

    real *ELJ6;   // składowa (σ/r)^6 energii LJ – suma po atomach (GPU)
    real *ELJ12;  // składowa (σ/r)^12 energii LJ – suma po atomach (GPU)
    real *ELJ;    // pełna energia LJ (ε[(σ/r)^12 - 2(σ/r)^6]) – suma po atomach (GPU)
};

struct parameters {
    real Rc, BoxSize, Buffer;
    int GridSize, CoreSize;
    bool PBC;
};

struct elj3 { real ELJ, ELJ6, ELJ12; };
struct energy3 { real Epot, Ekin, Etot; };

struct DisplacementStats {
    real avg_disp, max_disp;
    bool rebuild_needed;
};

//struct energy3 { real Epot, Ekin, Etot; };

constexpr int MAX_ATOMS_PER_CELL = 96;

struct CellAtomData {
    real X[MAX_ATOMS_PER_CELL], Y[MAX_ATOMS_PER_CELL], Z[MAX_ATOMS_PER_CELL];
    real Vx[MAX_ATOMS_PER_CELL], Vy[MAX_ATOMS_PER_CELL], Vz[MAX_ATOMS_PER_CELL];
    real Fx[MAX_ATOMS_PER_CELL], Fy[MAX_ATOMS_PER_CELL], Fz[MAX_ATOMS_PER_CELL];
    real X0[MAX_ATOMS_PER_CELL], Y0[MAX_ATOMS_PER_CELL], Z0[MAX_ATOMS_PER_CELL];
    int Natoms;
};

struct GridCells {
    int D, G, Ncells;
    real CellSize;
    CellAtomData* Cells;
    int* GlobalIndex;
};

struct CellAtomDataGPU {
    real X[MAX_ATOMS_PER_CELL], Y[MAX_ATOMS_PER_CELL], Z[MAX_ATOMS_PER_CELL];
    real ELJ12[MAX_ATOMS_PER_CELL], ELJ6[MAX_ATOMS_PER_CELL], ELJ[MAX_ATOMS_PER_CELL];
    real Fx[MAX_ATOMS_PER_CELL], Fy[MAX_ATOMS_PER_CELL], Fz[MAX_ATOMS_PER_CELL];
    real X0[MAX_ATOMS_PER_CELL], Y0[MAX_ATOMS_PER_CELL], Z0[MAX_ATOMS_PER_CELL];
    int Natoms;
};


//struct DisplacementStats { real avg_disp, max_disp; bool rebuild_needed; };

void ComputeCenterOfMass(const AtomData& Atoms, int Natoms, real& CMx, real& CMy, real& CMz);
void ShiftPositions(AtomData& Atoms, int Natoms, real dx, real dy, real dz);
void ComputeDisplacementStats(int Natoms, const AtomData& Atoms, real& avg_disp, real& max_disp);


