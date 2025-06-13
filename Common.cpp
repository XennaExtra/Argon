#include "Types.h"
#include "Common.h"

int MySize;
AtomData Argon;
parameters Pars;
energy3 E;
GridCells Grid;
TimerCounter T_total;
TimerCounter T_force;
TimerCounter T_rebuild;
TimerCounter T_ghost;
TimerCounter T_transfer_up;
TimerCounter T_transfer_down;
TimerCounter T_stats;
TimerCounter T_io;
TimerCounter T_verlet;

WallTimePoint TimerStart();
double TimerStop(WallTimePoint start);

WallTimePoint TimerStart() {
    return Clock::now();
}

double TimerStop(WallTimePoint start) {
    auto end = Clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return elapsed.count();
}

void PrintPhaseTime(const char* phase, double seconds) {
    printf("--- %s phase took %.6f seconds ---\n", phase, seconds);
}

// InitAllData()
// alokacja i zerowanie tablic AtomData
void InitAllData(AtomData *Atoms, int SIZ){
    int SizOfreal;
    MySize = 4 * SIZ * SIZ * SIZ;
    SizOfreal = sizeof(real);
    //printf("In InitAllData\n");

    Atoms->X = (real*) malloc(MySize * SizOfreal);
    Atoms->Y = (real*) malloc(MySize * SizOfreal);
    Atoms->Z = (real*) malloc(MySize * SizOfreal);
    Atoms->Vx = (real*) malloc(MySize * SizOfreal);
    Atoms->Vy = (real*) malloc(MySize * SizOfreal);
    Atoms->Vz = (real*) malloc(MySize * SizOfreal);
    Atoms->Fx = (real*) malloc(MySize * SizOfreal);
    Atoms->Fy = (real*) malloc(MySize * SizOfreal);
    Atoms->Fz = (real*) malloc(MySize * SizOfreal);
    Atoms->X0 = (real*) malloc(MySize * SizOfreal);
    Atoms->Y0 = (real*) malloc(MySize * SizOfreal);
    Atoms->Z0 = (real*) malloc(MySize * SizOfreal);
    Atoms->track = (bool*) malloc(MySize * sizeof(bool));
    Atoms->Natoms = MySize;

    for (int i = 0; i < MySize; i++) {
        Atoms->Fx[i] = Atoms->Fy[i] = Atoms->Fz[i] = 0.0f;
        Atoms->X0[i] = 0.0f;
        Atoms->Y0[i] = 0.0f;
        Atoms->Z0[i] = 0.0f;
    }
    //printf("Finished InitAllData(), Natoms = %d \n",Atoms->Natoms);
}

// Build()
// budowa początkowego układu (np. siatka FCC), ustalanie parametrów
int Build(AtomData& Atoms, int SIZ, parameters& Par) {
    //printf("In Build()\n");
    fflush(stdout);
    constexpr real T_ref = 5.0;
    constexpr real alpha_box = 0.0015;
    real scale = 1.0f + alpha_box * (T0 - T_ref);
    const real box_fcc = sigma * sqrt(2.0f) * scale / 2;
    Par.BoxSize = SIZ*box_fcc*2;

    printf("Building atoms\n");
    fflush(stdout);
    int i = 0;
    for (int ix = 0; ix < SIZ; ix++) {
        for (int iy = 0; iy < SIZ; iy++) {
            for (int iz = 0; iz < SIZ; iz++) {

                bool is_edge = (ix == 0 || ix == SIZ - 1 ||
                                iy == 0 || iy == SIZ - 1 ||
                                iz == 0 || iz == SIZ - 1);

                // Atom 1
                Atoms.X[i] = box_fcc * (0.5 + ix * 2);
                Atoms.Y[i] = box_fcc * (0.5 + iy * 2);
                Atoms.Z[i] = box_fcc * (0.5 + iz * 2);
                Atoms.track[i] = !is_edge;
                i++;

                // Atom 2
                Atoms.X[i] = box_fcc * (0.5 + ix * 2);
                Atoms.Y[i] = box_fcc * (1.5 + iy * 2);
                Atoms.Z[i] = box_fcc * (1.5 + iz * 2);
                Atoms.track[i] = !is_edge;
                i++;

                // Atom 3
                Atoms.X[i] = box_fcc * (1.5 + ix * 2);
                Atoms.Y[i] = box_fcc * (0.5 + iy * 2);
                Atoms.Z[i] = box_fcc * (1.5 + iz * 2);
                Atoms.track[i] = !is_edge;
                i++;

                // Atom 4
                Atoms.X[i] = box_fcc * (1.5 + ix * 2);
                Atoms.Y[i] = box_fcc * (1.5 + iy * 2);
                Atoms.Z[i] = box_fcc * (0.5 + iz * 2);
                Atoms.track[i] = !is_edge;
                i++;
            }
        }
    }
    //printf("Almost there \n");
    Atoms.Natoms = i;
    for (int i = 0; i < Atoms.Natoms; i++) {
        Atoms.Fx[i] = Atoms.Fy[i] = Atoms.Fz[i] = 0.0;
        Atoms.Vx[i] = Atoms.Vy[i] = Atoms.Vz[i] = 0.0;
        Atoms.X0[i] = Atoms.Y0[i] = Atoms.Z0[i] = 0.0;
    }
    //printf("Done \n");
    return Atoms.Natoms;
}
// Start()
// inicjalizacja prędkości zgodnie z rozkładem Maxwella
real Start(int Natoms, AtomData& Atoms, real Vmax){
    printf("In start\n");
    real CMx, CMy, CMz;
    ComputeCenterOfMass(Atoms, Natoms, CMx, CMy, CMz);
    ShiftPositions(Atoms, Natoms, -CMx, -CMy, -CMz);
    printf("Vmax %f \n",Vmax);
    real MVx = 0, MVy = 0, MVz = 0;
    real Ekin;

    for (int i = 0; i < Natoms; i++) {
        Atoms.Vx[i] = Vmax * (2.0f * rand() / RAND_MAX - 1.0f);
        Atoms.Vy[i] = Vmax * (2.0f * rand() / RAND_MAX - 1.0f);
        Atoms.Vz[i] = Vmax * (2.0f * rand() / RAND_MAX - 1.0f);
        MVx += Atoms.Vx[i];
        MVy += Atoms.Vy[i];
        MVz += Atoms.Vz[i];
    }
    MVx /= Natoms;
    MVy /= Natoms;
    MVz /= Natoms;
    // Usuwanie globalnego ruchu translacyjnego i rotacyjnego – symulacja powinna startować
    // bez przesunięcia środka masy ani rotacji względem osi.
    
    // Usuwanie ruchu środka masy
    for (int i = 0; i < Natoms; i++) {
        Atoms.Vx[i] -= MVx;
        Atoms.Vy[i] -= MVy;
        Atoms.Vz[i] -= MVz;
    }

    // Zerowanie momentu pędu
    real Lx = 0.0, Ly = 0.0, Lz = 0.0;
    real Ix = 0.0, Iy = 0.0, Iz = 0.0;
    for (int i = 0; i < Natoms; i++) {
        real x = Atoms.X[i], y = Atoms.Y[i], z = Atoms.Z[i];
        real vx = Atoms.Vx[i], vy = Atoms.Vy[i], vz = Atoms.Vz[i];
        Lx += m * (y * vz - z * vy);
        Ly += m * (z * vx - x * vz);
        Lz += m * (x * vy - y * vx);
        Ix += m * (y*y + z*z);
        Iy += m * (x*x + z*z);
        Iz += m * (x*x + y*y);
    }

    real wx = Lx / Ix;
    real wy = Ly / Iy;
    real wz = Lz / Iz;

    for (int i = 0; i < Natoms; i++) {
        real x = Atoms.X[i], y = Atoms.Y[i], z = Atoms.Z[i];
        Atoms.Vx[i] -= wy * z - wz * y;
        Atoms.Vy[i] -= wz * x - wx * z;
        Atoms.Vz[i] -= wx * y - wy * x;
    }

    // Przywrócenie położenia układu
    ShiftPositions(Atoms, Natoms, CMx, CMy, CMz);

    Ekin = 0.0;
    for (int i = 0; i < Natoms; i++) {
        Ekin += (Atoms.Vx[i]*Atoms.Vx[i] + Atoms.Vy[i]*Atoms.Vy[i] + Atoms.Vz[i]*Atoms.Vz[i]);
    }
    printf("Start done\n");
    return 0.5 * m * Ekin;
}


// ComputeCenterOfMass()
// obliczanie środka masy układu
void ComputeCenterOfMass(const AtomData& Atoms, int Natoms, real& CMx, real& CMy, real& CMz) {
    CMx = CMy = CMz = 0.0f;
    for (int i = 0; i < Natoms; i++) {
        CMx += Atoms.X[i];
        CMy += Atoms.Y[i];
        CMz += Atoms.Z[i];
    }
    CMx /= Natoms;
    CMy /= Natoms;
    CMz /= Natoms;
}

// ShiftPositions()
// przesunięcie wszystkich atomów o dany wektor
void ShiftPositions(AtomData& Atoms, int Natoms, real dx, real dy, real dz) {
    for (int i = 0; i < Natoms; i++) {
        Atoms.X[i] += dx;
        Atoms.Y[i] += dy;
        Atoms.Z[i] += dz;
    }
}

// ComputeDisplacementStats()
// obliczanie średniego i maksymalnego przemieszczenia atomów
void ComputeDisplacementStats(int Natoms, const AtomData& Atoms, real& avg_disp, real& max_disp) {
    // Uwaga: ta wersja działa tylko na danych z AtomData (tryb prosty PBC)
    // W trybie linked-cell (LC) należy użyć osobnej wersji śledzenia przesunięć,
    // operującej na strukturze GridCells i ograniczonej do realnych komórek.
    // Ta funkcja w LC moze być użyta tylko do gromadzenia statystyki, 
    // z częstotliwością równą odświeżaniu pozycji w AtomData. 
    real sum_disp2 = 0.0f;
    real max_disp2 = 0.0f;
    int tracked = 0;

    for (int i = 0; i < Natoms; i++) {
        if (!Atoms.track[i]) continue;

        real dx = Atoms.X[i] - Atoms.X0[i];
        real dy = Atoms.Y[i] - Atoms.Y0[i];
        real dz = Atoms.Z[i] - Atoms.Z0[i];
        real disp2 = dx*dx + dy*dy + dz*dz;

        sum_disp2 += disp2;
        if (disp2 > max_disp2) max_disp2 = disp2;
        tracked++;
    }

    avg_disp = sqrt(sum_disp2 / tracked);
    max_disp = sqrt(max_disp2);
}

// FreeAll()
// zwalnianie pamięci AtomData
void FreeAll(AtomData Atoms){
    
    free(Atoms.X);
    free(Atoms.Y);
    free(Atoms.Z);
    free(Atoms.Fx);
    free(Atoms.Fy);
    free(Atoms.Fz);
    free(Atoms.Vx);
    free(Atoms.Vy);
    free(Atoms.Vz);
    free(Atoms.X0);
    free(Atoms.Y0);
    free(Atoms.Z0);
    free(Atoms.track);
}