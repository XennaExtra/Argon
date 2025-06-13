#pragma once
#include "Types.h"





extern GridCells Grid;

struct Shift3D { real sx, sy, sz; };
struct cell3 { int x, y, z; };

inline int cellIndex(cell3 c, int G) { return c.x + G * (c.y + G * c.z); }
inline int ghost_src(int i, int D, int G) { return i + (i == 0) * D - (i == G - 1) * (G - 2); }

// LC functions
void InitCells(GridCells* Grid, int D, real CellSize);
void FreeCells(GridCells* Grid);
void AssignAtomsToCells(const AtomData& Argon, int Natoms, GridCells* Grid);
void CopyShadow(GridCells* Grid, int dst, int src, Shift3D shift);
void FillGhostCells(GridCells* Grid);
void GatherAtomsFromCells(AtomData& Argon, int Natoms, const GridCells* Grid);
void GatherPositionsFromCells(AtomData& Argon, int Natoms, const GridCells* Grid);
void GatherVelocitiesFromCells(AtomData& Argon, int Natoms, const GridCells* Grid);
real MaxDisplacement_LC(const GridCells& Grid);
void UpdateHalfVelocity_LC(GridCells& Grid, real DT2);
DisplacementStats UpdatePositions_LC(GridCells& Grid, real dt, parameters Par);
real UpdateKinetic_LC(const GridCells& Grid);
energy3 Verlet_LC(int Natoms, AtomData& Argon, GridCells& Grid, real dt, parameters Par, elj3* OutEpot);
void RunHeatingPhase_LC(int Natoms, AtomData& Atoms, GridCells& Grid, parameters Par, real Vmax, real t0);
void RunThermalisationPhase_LC(int Natoms, AtomData& Atoms, GridCells& Grid, parameters Par, real dt, int Term_, int save);
void RunProductionPhase_LC(int Natoms, AtomData& Atoms, GridCells& Grid, parameters Par, real dt, int Term_, int save);
void Simulation_LC(int SIZ, real dt, int Steps, int save, AtomData& Atoms, GridCells& Grid, parameters Par, real t0, bool doThermalize);
real Restart(AtomData& Argon, int Natoms, GridCells& Grid, real Vmax);
