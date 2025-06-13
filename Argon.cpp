#include "Types.h"
#include "Common.h"
#include "Argon.h"

#ifdef USE_GPU
#include "GPU.h"
#endif

#include "Forces.h"

extern int MySize;
extern AtomData Argon;
extern parameters Pars;
extern energy3 E;
extern GridCells Grid;
extern TimerCounter T_total;
extern TimerCounter T_force;
extern TimerCounter T_rebuild;
extern TimerCounter T_ghost;
extern TimerCounter T_transfer_up;
extern TimerCounter T_transfer_down;
extern TimerCounter T_stats;
extern TimerCounter T_io;
extern TimerCounter T_verlet;



// Verlet()
// główna pętla integracyjna (Verlet), wywołuje ForcesEnergiesWrapper()
energy3 Verlet(int Natoms, AtomData& Atoms, real DT2, parameters Par, elj3* OutEpot){
    energy3 locE = {} ;
    elj3 locELJ = {} ;

    T_verlet.start();

    // 1. Pierwsza połowa kroku prędkości (v ← v + 0.5*a*dt)
    UpdateHalfVelocity(Natoms, Atoms, DT2);

    // 2. Pozycje (r ← r + v*dt)
    UpdatePositions(Natoms, Atoms, DT2);

    // 3. Korekta pozycji przy PBC – jeśli atom wyszedł poza pudełko, zawijamy go z powrotem
    if (Par.PBC) {
        FixPositions(Natoms, Atoms, Par.BoxSize);
    }

    // 4. Siły i potencjał (a = f/m)
    T_force.start();
    locELJ = ForcesEnergiesWrapper(Natoms, Atoms, Par);
    T_force.add();
    *OutEpot = locELJ;  // energia przekazana na zewnątrz

    // 5. Druga połowa kroku prędkości
    UpdateHalfVelocity(Natoms, Atoms, DT2);

    // 6. Energetyka
    locE.Epot = locELJ.ELJ;
    locE.Ekin = UpdateKinetic(Natoms, Atoms);
    locE.Etot = locE.Epot + locE.Ekin;

    T_verlet.add();
    return locE;
}

// UpdateKinetic()
// obliczanie energii kinetycznej
real UpdateKinetic(int Natoms, AtomData& Atoms){
    real DEkin;
    real locEkin=0;
    
    for (int i=0;i<Natoms;i++) {
        DEkin = Atoms.Vx[i]*Atoms.Vx[i];
        DEkin += Atoms.Vy[i]*Atoms.Vy[i];
        DEkin += Atoms.Vz[i]*Atoms.Vz[i];
        DEkin = DEkin*m*0.5;
        // computed kinetic energy of the i-th atom
        locEkin+=DEkin;
        // computed kinetic energy of all atoms up to i-th
    }
    return locEkin;
};

// UpdateHalfVelocity()
// aktualizacja prędkości o połowę kroku czasowego
void UpdateHalfVelocity(int Natoms, AtomData& Atoms, real DT2){
    real ax, ay, az;
    real halfDT=DT2*0.5;
    
    for (int i=0;i<Natoms;i++) {
        ax = Atoms.Fx[i]/m;
        ay = Atoms.Fy[i]/m;
        az = Atoms.Fz[i]/m;
        Atoms.Vx[i]=Atoms.Vx[i]+ax*halfDT;
        Atoms.Vy[i]=Atoms.Vy[i]+ay*halfDT;
        Atoms.Vz[i]=Atoms.Vz[i]+az*halfDT;
    }
}

// UpdatePositions()
// aktualizacja położeń atomów
void UpdatePositions(int Natoms, AtomData& Atoms, real DT2){
    for (int i=0;i<Natoms;i++) {
        Atoms.X[i]=Atoms.X[i]+Atoms.Vx[i]*DT2;
        Atoms.Y[i]=Atoms.Y[i]+Atoms.Vy[i]*DT2;
        Atoms.Z[i]=Atoms.Z[i]+Atoms.Vz[i]*DT2;
    }
    // computed new positions at new full step - step #2
};


// RunHeatingPhase()
// podgrzewanie układu przez zadany czas

void RunHeatingPhase(int Natoms, AtomData& locAtoms, parameters Par, real Vmax, real t0) {
    auto start = TimerStart();
    constexpr int N_Heating = 3;
    constexpr real HeatingTime = 0.5f; // ps
    constexpr real HeatingStep = 0.02f;
    const int HeatingSteps = int(HeatingTime / HeatingStep);

    real Ekin_temp = Start(Natoms, locAtoms, Vmax * sqrtf(0.5f));
    printf("Thermalization Start / %d at T ≈ %.2f K (is %f)\n", N_Heating, t0 * sqrt(0.5f),Ekin_temp);

    for (int step_ = 0; step_ < HeatingSteps; step_++) {
        elj3 dummyEpot = {} ;
        energy3 Etmp = Verlet(Natoms, locAtoms, HeatingStep, Par, &dummyEpot);
        real Tinst = 2 * Etmp.Ekin / ((3 * Natoms - 3) * kB);
        if (step_ == HeatingSteps - 1) {
            printf("  -> Final temperature for starting cycle: %.2f K\n", Tinst);
        }
    }
    for (int cycle = 0; cycle < N_Heating; cycle++) {
        printf("Thermalization cycle %d / %d at T ≈ %.2f K\n", cycle + 1, N_Heating, t0);
        real Ekin_temp = Start(Natoms, locAtoms, Vmax);
        for (int step_ = 0; step_ < HeatingSteps; step_++) {
            elj3 dummyEpot = {} ;
            energy3 Etmp = Verlet(Natoms, locAtoms, HeatingStep, Par, &dummyEpot);
            real Tinst = 2 * Etmp.Ekin / ((3 * Natoms - 3) * kB);
            if (step_ == HeatingSteps - 1) {
                printf("  -> Final temperature in cycle %d: %.2f K\n", cycle + 1, Tinst);
            }
        }
    }
    PrintPhaseTime("Heating", TimerStop(start));
}
// RunThermalisationPhase()
// termalizacja układu (np. regulacja temperatury)
void RunThermalisationPhase(int Natoms, AtomData& Atoms, parameters Par, real dt, int Term_, int save) {
    auto start = TimerStart();
    FILE *CSV = fopen("therm_energies.csv", "w");
    FILE *DISP = fopen("therm_disp.csv", "w");
    fprintf(CSV, "  Time,     Etot,      Epot,      Ekin,     Elj12,      Elj6,    Temp\n");
    fprintf(DISP, "  Time, AvgDisp, MaxDisp\n");

    printf("\n\nThermalisation for %d steps\n\n",Term_);
    double Eavg = 0.0, Evar = 0.0, Estdev = 0.0, Exx = 0.0, Ex = 0.0;
    double Emin = 1e20, Emax = -1e20, Estart = 0.0, Eend = 0.0;
    energy3 locE = {};
    elj3 Epot = {};
    real Temp;
    for (int i = 0; i < Natoms; i++) {
        Atoms.X0[i] = Atoms.X[i];
        Atoms.Y0[i] = Atoms.Y[i];
        Atoms.Z0[i] = Atoms.Z[i];
    }
    for (int step = 0; step < Term_; step++) {
        locE = Verlet(Natoms, Atoms, dt, Par, &Epot);
        if (step == 0) Estart = locE.Etot;
        if ((step % save) == 0) {
            real avg_disp = 0.0f, max_disp = 0.0f;
            ComputeDisplacementStats(Natoms, Atoms, avg_disp, max_disp);
            Temp = 2 * locE.Ekin / ((3 * Natoms - 3) * kB);
            fprintf(DISP, "%8.3f, %8.5f, %8.5f \n", step * dt, avg_disp, max_disp);
            fprintf(CSV, "%8.3f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n",step * dt, locE.Etot, locE.Epot, locE.Ekin, Epot.ELJ12, Epot.ELJ6, Temp);
            printf("Step %i, ETOT = %f\n", step, locE.Etot);
        }
        Ex += locE.Etot;
        Exx += locE.Etot * locE.Etot;
        if (locE.Etot < Emin) Emin = locE.Etot;
        if (locE.Etot > Emax) Emax = locE.Etot;
    }
    Eend = locE.Etot;
    Eavg = Ex / Term_;
    Evar = Exx / Term_ - Eavg * Eavg;
    Estdev = sqrt(Evar);
    Temp = 2 * locE.Ekin / ((3 * Natoms - 3) * kB);
    fprintf(CSV, "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n",Term_ * dt, locE.Etot, locE.Epot, locE.Ekin, Epot.ELJ12, Epot.ELJ6, Temp);
    printf("Step %i, ETOT = %lf \n", Term_, locE.Etot);
    printf("Etot = %f +/- %f, Emin=%f  Emax=%f, DE=%f ", Eavg, Estdev, Emin, Emax, Emax - Emin);
    printf("Estart=%f, Eend=%f, DELTA=%f \n", Estart, Eend, Estart - Eend);
    fclose(CSV);
    fclose(DISP);
    PrintPhaseTime("Thermalisation", TimerStop(start));
}

// RunProductionPhase()
// faza produkcyjna – właściwa symulacja
void RunProductionPhase(int Natoms, AtomData& Atoms, parameters Par, real dt, int Term_, int save) {
    auto start = TimerStart();
  
    FILE *CSV = fopen(Fname, "w");
    FILE *DISP = fopen("displacement.csv", "w");
    fprintf(CSV, "  Time,     Etot,      Epot,      Ekin,     Elj12,      Elj6,    Temp\n");
    fprintf(DISP, "  Time, AvgDisp, MaxDisp\n");

    double Eavg = 0.0, Evar = 0.0, Estdev = 0.0, Exx = 0.0, Ex = 0.0;
    double Emin = 1e20, Emax = -1e20, Estart = 0.0, Eend = 0.0;
    energy3 locE = {};
    elj3 Epot = {};
    real Temp;

    printf("\n\nProduction phase for %d steps\n\n",Term_);
    for (int i = 0; i < Natoms; i++) {
        Atoms.X0[i] = Atoms.X[i];
        Atoms.Y0[i] = Atoms.Y[i];
        Atoms.Z0[i] = Atoms.Z[i];
    }
    for (int step = 0; step < Term_; step++) {
        locE = Verlet(Natoms, Atoms, dt, Par, &Epot);
        if (step == 0) Estart = locE.Etot;
        if ((step % save) == 0) {
            real avg_disp = 0.0f, max_disp = 0.0f;
            T_stats.start();
            ComputeDisplacementStats(Natoms, Atoms, avg_disp, max_disp);
            T_stats.add();
            Temp = 2 * locE.Ekin / ((3 * Natoms - 3) * kB);
            T_io.start();
            fprintf(DISP, "%8.3f, %8.5f, %8.5f \n", step * dt, avg_disp, max_disp);
            fprintf(CSV, "%8.3f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n",step * dt, locE.Etot, locE.Epot, locE.Ekin, Epot.ELJ12, Epot.ELJ6, Temp);
            printf("Step %i, ETOT = %f\n", step, locE.Etot);
            T_io.add();
        }
        Ex += locE.Etot;
        Exx += locE.Etot * locE.Etot;
        if (locE.Etot < Emin) Emin = locE.Etot;
        if (locE.Etot > Emax) Emax = locE.Etot;
    }
    Eend = locE.Etot;
    Eavg = Ex / Term_;
    Evar = Exx / Term_ - Eavg * Eavg;
    Estdev = sqrt(Evar);
    Temp = 2 * locE.Ekin / ((3 * Natoms - 3) * kB);
    T_io.start();
    fprintf(CSV, "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n",Term_ * dt, locE.Etot, locE.Epot, locE.Ekin, Epot.ELJ12, Epot.ELJ6, Temp);
    printf("Step %i, ETOT = %lf \n", Term_, locE.Etot);
    printf("Etot = %f +/- %f, Emin=%f  Emax=%f, DE=%f ", Eavg, Estdev, Emin, Emax, Emax - Emin);
    printf("Estart=%f, Eend=%f, DELTA=%f \n", Estart, Eend, Estart - Eend);
    fclose(CSV);
    fclose(DISP);
    T_io.add();
    
    PrintPhaseTime("Production", TimerStop(start));
}

// Simulation()
// pełna procedura: budowa układu, termalizacja, produkcja
void Simulation(int SIZ, real dt, int Steps, int save, AtomData& Atoms, parameters Par, real t0, bool doThermalize) {
    int Natoms = Atoms.Natoms;
    int therm_steps = 100;
    real therm_time_step = 0.02;
    T_total.start();
    real Vmax = sqrt(3 * kB * t0 / m);
    RunHeatingPhase(Natoms, Atoms, Par, Vmax, t0);
    RunThermalisationPhase(Natoms, Atoms, Par, therm_time_step, therm_steps, save);
    RunProductionPhase(Natoms, Atoms, Par, dt, Steps, save);
    T_total.add();
    FreeAll(Atoms);
}

// FixPositions()
// dopasowanie atomów do pudełka (np. po przesunięciu)
void FixPositions(int Natoms, AtomData& Atoms, real Box) {
    for (int i = 0; i < Natoms; i++) {
        // X-koordynata
        if (Atoms.X[i] < 0.0f)
            Atoms.X[i] += Box;
        else if (Atoms.X[i] >= Box)
            Atoms.X[i] -= Box;

        // Y-koordynata
        if (Atoms.Y[i] < 0.0f)
            Atoms.Y[i] += Box;
        else if (Atoms.Y[i] >= Box)
            Atoms.Y[i] -= Box;

        // Z-koordynata
        if (Atoms.Z[i] < 0.0f)
            Atoms.Z[i] += Box;
        else if (Atoms.Z[i] >= Box)
            Atoms.Z[i] -= Box;
    }
}

// PrintPhaseTime()
// pomocniczy wypis statystyk czasowych


