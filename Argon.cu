#include "Types.h"
#include "Common.h"
#include "Argon.h"
#include "GPU.h"
#include "Forces_GPU.h"
//#include "Forces.h"

int MySize;
extern AtomData Argon;
extern AtomDataGPU ArgonGPU;
extern AtomDataGPU ArgonCPU;

extern CellAtomData cellArgon;
extern CellAtomDataGPU cellArgonGPU;
extern CellAtomDataGPU cellArgonCPU;

#define BW 32



extern parameters Pars;
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
TimerCounter T_kernel;
TimerCounter T_copy2gpu; 
TimerCounter T_copy2cpu;



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

void InitAllDataGPU(AtomDataGPU *GA, AtomDataGPU *CA, int SIZ){
    printf("In init_all_data GPU\n");
    int Sreal;    
    MySize=4*SIZ*SIZ*SIZ; 
    Sreal=sizeof(real);
    
    int BufSize = MySize*(MySize + BW-1)/BW;
 
    CALL_CUDA(cudaMalloc((void**)&GA->X , MySize*Sreal));
    CALL_CUDA(cudaMalloc((void**)&GA->Y , MySize*Sreal)); 
    CALL_CUDA(cudaMalloc((void**)&GA->Z , MySize*Sreal)); 
    CALL_CUDA(cudaMalloc((void**)&GA->Fx , MySize*Sreal)); 
    CALL_CUDA(cudaMalloc((void**)&GA->Fy , MySize*Sreal)); 
    CALL_CUDA(cudaMalloc((void**)&GA->Fz , MySize*Sreal)); 
    CALL_CUDA(cudaMalloc((void**)&GA->ELJ6 , MySize*Sreal)); 
    CALL_CUDA(cudaMalloc((void**)&GA->ELJ12 , MySize*Sreal)); 
    CALL_CUDA(cudaMalloc((void**)&GA->ELJ , MySize*Sreal)); 


    // alokowanie pamięci na CPU na tablice energii potencjalnych 

    //real *elj6, *elj12, *elj;
    //elj6=(real*) malloc(MySize*Sreal);
    //elj12=(real*) malloc(MySize*Sreal);
    //elj=(real*) malloc(MySize*Sreal);

    // przypisanie zaalokowanych tablic do struktury ArgonCPU
    CA->ELJ6= (real*) malloc(MySize*Sreal);
    CA->ELJ12= (real*) malloc(MySize*Sreal);
    CA->ELJ= (real*) malloc(MySize*Sreal);
    CA->X = (real*) malloc(MySize*Sreal);
    CA->Y = (real*) malloc(MySize*Sreal);
    CA->Z = (real*) malloc(MySize*Sreal);
    CA->Fx = (real*) malloc(MySize*Sreal);
    CA->Fy = (real*) malloc(MySize*Sreal);
    CA->Fz = (real*) malloc(MySize*Sreal);

    // Alokacja buforów pośrednich 

    //CALL_CUDA(cudaMalloc((void**)&host_buffer, 5*BufSize*Sreal)) 
    //CALL_CUDA(cudaMalloc((void**)&guest_buffer, 5*BufSize*Sreal)) 
    //CALL_CUDA(cudaMalloc((void**)&gHostELJ12 , BufSize*Sreal)) 
    //CALL_CUDA(cudaMalloc((void**)&gHostFx , BufSize*Sreal)) 
    //CALL_CUDA(cudaMalloc((void**)&gHostFy , BufSize*Sreal)) 
    //CALL_CUDA(cudaMalloc((void**)&gHostFz , BufSize*Sreal)) 

    //cHostELJ6=(real*) malloc(BufSize*Sreal);
    //cHostELJ12=(real*) malloc(BufSize*Sreal);
    //cHostFx=(real*) malloc(BufSize*Sreal);
    //cHostFy=(real*) malloc(BufSize*Sreal);
    //cHostFz=(real*) malloc(BufSize*Sreal);
    
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

    printf("sigma = %f,box_fcc = %f\n\n",sigma,box_fcc);

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
    printf("Done \n");
    //for (int i =0;i<10;i++) {
    //    printf("%d X=%f Y=%f Z=%f \n",i,Atoms.X[i],Atoms.Y[i],Atoms.Z[i]);
    //}
    return Atoms.Natoms;
}

// Start()
// inicjalizacja prędkości zgodnie z rozkładem Maxwella
real Start(int Natoms, AtomData& Atoms, real Vmax){
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

    return 0.5 * m * Ekin;
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

// Verlet()
// główna pętla integracyjna (Verlet), wywołuje ForcesEnergiesWrapper()
energy3 Verlet(int Natoms, AtomData& Atoms, real DT2, parameters Par, elj3* OutEpot){
    energy3 locE = {} ;
    elj3 locELJ = {} ;
    //msg("Verlet start\n");
    T_verlet.start();
    //msg("T_verlet.start");
    // 1. Pierwsza połowa kroku prędkości (v ← v + 0.5*a*dt)
    UpdateHalfVelocity(Natoms, Atoms, DT2);
    //msg("UpdateHalfVelocity");
    // 2. Pozycje (r ← r + v*dt)
    UpdatePositions(Natoms, Atoms, DT2);
    //msg("UpdatePositions");
    // 3. Korekta pozycji przy PBC – jeśli atom wyszedł poza pudełko, zawijamy go z powrotem
    if (Par.PBC) {
        FixPositions(Natoms, Atoms, Par.BoxSize);
        msg("FixPositions");
    }

    // 4. Siły i potencjał (a = f/m)
    T_force.start();
    //msg("T_force.start");
    //locELJ = ForcesEnergies_GPU(Natoms, Atoms, Par, );
    //locELJ = ForcesEnergies_GPU(Natoms, Atoms, ArgonGPU, ArgonCPU,Pars); 
    locELJ = ForcesEnergies_GPU(Natoms, Atoms, ArgonGPU, ArgonCPU,Pars); 
    T_force.add();
    //msg("T_force.add");
    *OutEpot = locELJ;  // energia przekazana na zewnątrz

    // 5. Druga połowa kroku prędkości
    UpdateHalfVelocity(Natoms, Atoms, DT2);
    //msg("UpdateHalfVelocity");
    // 6. Energetyka
    locE.Epot = locELJ.ELJ;
    locE.Ekin = UpdateKinetic(Natoms, Atoms);
    locE.Etot = locE.Epot + locE.Ekin;

    T_verlet.add();
    //msg("T_verlet.add");
    //msg("Verlet");
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
    real X, Y, Z;
    real X2, Y2, Z2;
    X=Y=Z=X2=Y2=Z2=0.0f;
    for (int i=0;i<Natoms;i++) {
        Atoms.X[i]=Atoms.X[i]+Atoms.Vx[i]*DT2; X+=Atoms.X[i];X2+=Atoms.X[i]*Atoms.X[i];
        Atoms.Y[i]=Atoms.Y[i]+Atoms.Vy[i]*DT2; Y+=Atoms.Y[i];Y2+=Atoms.Y[i]*Atoms.Y[i];
        Atoms.Z[i]=Atoms.Z[i]+Atoms.Vz[i]*DT2; Z+=Atoms.Z[i];Z2+=Atoms.Z[i]*Atoms.Z[i];
    }
    // computed new positions at new full step - step #2
    X=X/Natoms; X2=sqrt(X2)/Natoms;
    Y=Y/Natoms; Y2=sqrt(Y2)/Natoms;
    Z=Z/Natoms; Z2=sqrt(Z2)/Natoms;
    //printf("X = %f, Y= %f, Z=%f, R = %f\n",X,Y,Z,sqrt(X2+Y2+Z2));
};

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
    real therm_time_step = 0.01;
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

// PrintPhaseTime()
// pomocniczy wypis statystyk czasowych

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
