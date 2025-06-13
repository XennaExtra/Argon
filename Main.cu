// Plik Main.cpp
// Wersja main() dla modelu prostego (CPU i GPU naive/tiled/shared/newton)
// Używa struktur AtomData i parameters
// Wywołuje: Simulation()
// Obsługuje argumenty wejściowe (opcjonalnie), wypisuje czas itd.

#include "Types.h"
#include "Common.h"
#include "Argon.h"
//#include "Forces.h"
#include "GPU.h"

extern AtomData Argon;
extern AtomDataGPU ArgonGPU;   // struktura danych na GPU
extern AtomDataGPU ArgonCPU;   // kopia na CPU dla transferu (host mirror)

int main(int argc, char** argv) {
    int size, term;
    real dt, T;
    int save;

    srand(777); // Ustalanie ziarna losowego – deterministyczne wyniki

     if (argc != 6) {
        printf("Usage: %s Rozmiar Krok Liczba_kroków Temperatura Save\n", argv[0]);
        printf("Rozmiar - liczba komórek FCC w jednym kierunku\n");
        printf("Krok - krok czasowy symulacji - rozsądne wartości między 0.001 a 0.03\n");
        printf("Domyślnie wpisujemy 0.01\n");
        printf("Liczba kroków symulacji w fazie produkcyjnej\n");
        printf("Temperatura - temperatura symulacji. Wartości między 70 a 90 (Kelwinów\n");
        //printf("Grid - rozmiar siatki obliczeniowej\n");
        printf("Save - częstotliwość zapisu parametrów trajektorii\n");
        exit(1);
    }

    size = atoi(argv[1]);
    dt   = atof(argv[2]);
    term = atoi(argv[3]);
    T    = atof(argv[4]);
    save = atoi(argv[5]);
    if (T == 0.0) T = T0;

    save = atoi(argv[5]);
    if (save < 1) save = SAVE;

    Pars.Rc = 1.0; // Promień odcięcia oddziaływań LJ
    
    printf("Liczba komórek FCC: %i (układ: %i atomów)\n", size * size * size, 4 * size * size * size);
    printf("dt = %.4f, steps = %d, T = %.2f K, zapis co %d kroków\n", dt, term, T, save);


    // Alokacja i inicjalizacja danych atomów
    InitAllData(&Argon, size);
    InitAllDataGPU(&ArgonGPU,&ArgonCPU,size);
    Build(Argon, size, Pars);

    // Bufor na potrzeby linked-cell (niewykorzystywany w wersji naiwnej)
    Pars.Buffer = (Pars.BoxSize - Pars.CoreSize * Pars.Rc) / Pars.CoreSize;

    printf("Konfiguracja symulacji\n");
    printf("Rcut = %.3f, Liczba atomów = %d\n", Pars.Rc, Argon.Natoms);

    // Uruchomienie głównej pętli symulacyjnej (CPU-naive)
    Simulation(size, dt, term, save, Argon, Pars, T, true);

    // Raporty z pomiarów czasu wykonania
    T_total.report("Total");
    T_verlet.report("Krok Verlet", T_total.getTotal());
    T_force.report("Siły", T_total.getTotal());
    T_io.report("Zapis danych", T_total.getTotal());
    T_stats.report("Statystyki", T_total.getTotal());

    return 0;
}

