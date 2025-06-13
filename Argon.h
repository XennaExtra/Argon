#pragma once

#include "Types.h"

extern AtomData Argon;
extern parameters Pars;
extern energy3 E;

// CPU-naive functions
WallTimePoint TimerStart();
double TimerStop(WallTimePoint start);
void InitAllData(AtomData *Atoms, int SIZ);
int Build(AtomData& Atoms, int SIZ, parameters& Par);
real Start(int Natoms, AtomData& Atoms, real Vmax);
void FreeAll(AtomData Atoms);
energy3 Verlet(int Natoms, AtomData& Atoms, real DT2, parameters Par, elj3* OutEpot);
real UpdateKinetic(int Natoms, AtomData& Atoms);
void UpdateHalfVelocity(int Natoms, AtomData& Atoms, real DT2);
void UpdatePositions(int Natoms, AtomData& Atoms, real DT2);
void ComputeDisplacementStats(int Natoms, const AtomData& Atoms, real& avg_disp, real& max_disp);
void RunHeatingPhase(int Natoms, AtomData& Atoms, parameters Par, real Vmax, real t0);
void RunThermalisationPhase(int Natoms, AtomData& Atoms, parameters Par, real dt, int Term_, int save);
void RunProductionPhase(int Natoms, AtomData& Atoms, parameters Par, real dt, int Term_, int save);
void Simulation(int SIZ, real dt, int Steps, int save, AtomData& Atoms, parameters Par, real t0, bool doThermalize);
void FixPositions(int Natoms, AtomData& Atoms, real Box);
void ComputeCenterOfMass(const AtomData& Atoms, int Natoms, real& CMx, real& CMy, real& CMz);
void ShiftPositions(AtomData& Atoms, int Natoms, real dx, real dy, real dz);
void PrintPhaseTime(const char* phase, double seconds);
void msg(const char* name);