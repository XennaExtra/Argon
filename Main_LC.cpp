#include "Types.h"
#include "Common.h" 

#ifdef USE_GPU
#include "GPU.h"
#endif

#include "Argon.h"
#include "Argon_LC.h" 
#include "Forces.h"

int main(int argc, char** argv) {
    int size, term;
    real dt, T;
    int save;

    srand(777);

    if (argc != 7) {
        printf("Usage: %s Rozmiar Krok Liczba_kroków Temperatura Grid Save\n", argv[0]);
        printf("Rozmiar - liczba komórek FCC w jednym kierunku\n");
        printf("Krok - krok czasowy symulacji - rozsądne wartości między 0.001 a 0.03\n");
        printf("Domyślnie wpisujemy 0.01\n");
        printf("Liczba kroków w fazie produkcyjnej symulacji 0.01\n");
        printf("Temperatura - temperatura symulacji. Wartości między 70 a 90 (Kelwinów\n");
        printf("Grid - rozmiar siatki obliczeniowej\n");
        printf("Save - częstotliwość zapisu parametrów trajektorii\n");
        exit(1);
    }

    size = atoi(argv[1]);
    dt   = atof(argv[2]);
    term = atoi(argv[3]);
    T    = atof(argv[4]);
    Pars.CoreSize = atoi(argv[5]);
    save = atoi(argv[6]);
    if (T == 0.0) T = T0;

    Pars.Rc = 1.0;        // Promień odcięcia (ustalony na sztywno)
    //size = 24;            // Liczba komórek FCC w każdym wymiarze
    //Pars.CoreSize = 13;
    

    printf("Number of fcc cells: %i\n", size * size * size);
    printf("dt = %lf, term = %i, T = %lf\n", dt, term, T);

    InitAllData(&Argon, size);        // Alokacja i inicjalizacja struktury AtomData
    Build(Argon, size, Pars);         // Budowa sieci FCC i wyliczenie rozmiaru pudła

    Pars.Buffer = (Pars.BoxSize - Pars.CoreSize * Pars.Rc) / Pars.CoreSize;

    printf("Konfiguracja\n\n");
    printf("Rcut %f bufor %f, core grid size %d Natoms %d\n",Pars.Rc, Pars.Buffer, Pars.CoreSize, Argon.Natoms);

    // Inicjalizacja i konfiguracja siatki linked-cell
    InitCells(&Grid, Pars.CoreSize, Pars.Rc + Pars.Buffer);
    AssignAtomsToCells(Argon, Argon.Natoms, &Grid);
    FillGhostCells(&Grid);

    // Właściwa symulacja (fazowanie + produkcja)
    Simulation_LC(size, dt, term, save, Argon, Grid, Pars, T, true);

    // Raport czasu wykonania
    T_total.report("Total");
    T_verlet.report("Krok Verlet", T_total.getTotal());
    T_force.report("Siły", T_total.getTotal());
    T_io.report("Zapis danych", T_total.getTotal());
    T_stats.report("Statystyki", T_total.getTotal());
    T_rebuild.report("Rekalkulacja gridu", T_total.getTotal());
    T_ghost.report("Kopiowanie ghostów", T_total.getTotal());

    FreeCells(&Grid);
}
