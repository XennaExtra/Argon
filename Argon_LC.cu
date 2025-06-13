#include "Types.h"
#include "Common.h"
#include "GPU.h"
#include "Argon.h"
#include "Argon_LC.h"
#include "Forces_GPU.h"


// InitCells()
// alokacja siatki komórek
void InitCells(GridCells* Grid, int D, real CellSize) {
    // Inicjalizacja struktury GridCells.
    // Tworzy trójwymiarową siatkę komórek o rozmiarze (D+2)^3,
    // gdzie D to liczba realnych komórek, a dodatkowe warstwy to ghost cells.
    // Każda komórka otrzymuje prealokowaną strukturę do przechowywania danych atomów.
    Grid->D = D;
    Grid->G = D + 2;
    Grid->Ncells = Grid->G * Grid->G * Grid->G;
    Grid->CellSize = CellSize;
    printf("Core size: %d GridSize %D",Grid->D = D,Grid->G);

    // Alokacja: tablica komórek i liczników
    Grid->Cells = (CellAtomData*)malloc(Grid->Ncells * sizeof(CellAtomData));
    Grid->GlobalIndex = (int*)malloc(Grid->Ncells * MAX_ATOMS_PER_CELL * sizeof(int));
}

// FreeCells()
// zwolnienie pamięci siatki komórek
void FreeCells(GridCells* Grid) {
    // Zwolnienie pamięci zaalokowanej na strukturę GridCells.
    // Dotyczy zarówno danych atomów w komórkach, jak i indeksów globalnych.

    if (!Grid) return;

    free(Grid->Cells);
    //free(Grid->Counts);
    free(Grid->GlobalIndex);

    Grid->GlobalIndex = NULL;
    Grid->Cells = NULL;
    //Grid->Counts = NULL;
    Grid->Ncells = 0;
    Grid->D = 0;
    Grid->G = 0;
}

// AssignAtomsToCells()
// przypisanie atomów do komórek na podstawie położeń
void AssignAtomsToCells(const AtomData& Argon, int Natoms, GridCells* Grid) {
    // Przypisuje każdy atom do odpowiedniej komórki siatki GridCells.
    // Pozycje atomów są mapowane do indeksów komórek przy uwzględnieniu PBC (fmod).
    // Zmienna slot określa pozycję atomu w danej komórce (maksymalnie MAX_ATOMS_PER_CELL).
    // GlobalIndex przechowuje numer atomu w głównej strukturze, by umożliwić późniejszy zapis zwrotny.

    int G = Grid->G;
    real CellSize = Grid->CellSize;
    real BoxSize = CellSize*Grid->D;

    // Zerujemy liczniki atomów w komórkach
    for (int cid = 0; cid < Grid->Ncells; cid++) {
    //    Grid->Counts[i] = 0;
          CellAtomData* cell = &Grid->Cells[cid];
          cell->Natoms=0;
    }

    for (int i = 0; i < Natoms; i++) {

        real x = fmod(fmod(Argon.X[i], BoxSize) + BoxSize, BoxSize);
        int cx = (int)(x / CellSize);
        real y = fmod(fmod(Argon.Y[i], BoxSize) + BoxSize, BoxSize);
        int cy = (int)(y / CellSize);
        real z = fmod(fmod(Argon.Z[i], BoxSize) + BoxSize, BoxSize);
        int cz = (int)(z / CellSize);

        // Przesunięcie o +1, aby wejść w zakres [1 ... D]
        cx += 1; cy += 1; cz += 1;

        // Identyfikator komórki
        cell3 c = {cx, cy, cz};
        int cid = cellIndex(c, G);

        //printf(" Cell Index: %3d;  cx %3d cy %3d cz %3d \n",cid,cx,cy,cz);

        CellAtomData* cell = &Grid->Cells[cid];
        int slot = cell->Natoms;
        cell->Natoms++;

        if (cell->Natoms == MAX_ATOMS_PER_CELL) {
            printf("Error: too many atoms in cell %d\n", cid);
            exit(1);
        }

        cell->X[slot]  = Argon.X[i];
        cell->Y[slot]  = Argon.Y[i];
        cell->Z[slot]  = Argon.Z[i];

        cell->Vx[slot] = Argon.Vx[i];
        cell->Vy[slot] = Argon.Vy[i];
        cell->Vz[slot] = Argon.Vz[i];

        cell->X0[slot] = Argon.X[i];
        cell->Y0[slot] = Argon.Y[i];
        cell->Z0[slot] = Argon.Z[i];
        Grid->GlobalIndex[cid * MAX_ATOMS_PER_CELL + slot] = i;

    }
}

// CopyShadow()
// kopiowanie zawartości komórki sąsiadującej (ghost)
void CopyShadow(GridCells* Grid, int dst, int src, Shift3D shift) {
    // Kopiuje pozycje atomów z komórki src (komórka rzeczywista) 
    // do komórki obrazu (dst), z odpowiednim przesunięciem pozycji.
    // Używane przy budowie periodycznych obrazów komórek do implementacji PBC w siatce komórek.

    CellAtomData* src_cell = &Grid->Cells[src];
    CellAtomData* dst_cell = &Grid->Cells[dst];

    int count = src_cell->Natoms;
    dst_cell->Natoms = count;

    real dx = shift.sx;
    real dy = shift.sy;
    real dz = shift.sz;

    for (int k = 0; k < count; k++) {
        dst_cell->X[k] = src_cell->X[k] + dx;
        dst_cell->Y[k] = src_cell->Y[k] + dy;
        dst_cell->Z[k] = src_cell->Z[k] + dz;
    }
}


// FillGhostCells()
// uzupełnienie komórek ghostowych wokół siatki
void FillGhostCells(GridCells* Grid) {
    // Tworzy warstwy obrazów periodycznych  wzdłuż trzech osi siatki GridCells.
    // Kopiuje dane z odpowiednich brzegowych komórek realnych, dodając przesunięcie pozycji.
    // Dzięki temu każda komórka na brzegu ma dostęp do sąsiadów również w obecności PBC.
    int D = Grid->D;
    int G = Grid->G;
    real Box = D * Grid->CellSize;

    // ---------------------------
    // Obrazy w osi X: x = 0 i x = G-1
    // ---------------------------
    for (int x = 0; x < G; x += (G - 1)) {
        int src_x = ghost_src(x, D, G);
        real shift_x = (x == 0) ? -Box : +Box;
        for (int y = 0; y < G; ++y) {
            for (int z = 0; z < G; ++z) {
                int dst = cellIndex({x, y, z}, G);
                int src = cellIndex({src_x, y, z}, G);
                Shift3D shift = {shift_x, 0.0f, 0.0f};
                CopyShadow(Grid, dst, src, shift);
            }
        }
    }

    // ---------------------------
    // Obrazy w osi Y: y = 0 i y = G-1
    // ---------------------------
    for (int y = 0; y < G; y += (G - 1)) {
        int src_y = ghost_src(y, D, G);
        real shift_y = (y == 0) ? -Box : +Box;
        for (int x = 0; x < G; ++x) {
            for (int z = 0; z < G; ++z) {
                int dst = cellIndex({x, y, z}, G);
                int src = cellIndex({x, src_y, z}, G);
                Shift3D shift = {0.0f, shift_y, 0.0f};
                CopyShadow(Grid, dst, src, shift);
            }
        }
    }

    // ---------------------------
    // Obrazy w osi Z: z = 0 i z = G-1
    // ---------------------------
    for (int z = 0; z < G; z += (G - 1)) {
        int src_z = ghost_src(z, D, G);
        real shift_z = (z == 0) ? -Box : +Box;
        for (int x = 0; x < G; ++x) {
            for (int y = 0; y < G; ++y) {
                int dst = cellIndex({x, y, z}, G);
                int src = cellIndex({x, y, src_z}, G);
                Shift3D shift = {0.0f, 0.0f, shift_z};
                CopyShadow(Grid, dst, src, shift);
            }
        }
    }
}

// GatherAtomsFromCells()
// przeniesienie pozycji i prędkości z siatki do AtomData
void GatherAtomsFromCells(AtomData& Argon, int Natoms, const GridCells* Grid) {
    // Zrzuca dane pozycji i prędkości z komórek realnych z powrotem do głównej struktury AtomData.
    // Używane podczas rekonstrukcji siatki i przed propagacją pozycji.
    // Indeks globalny (GlobalIndex) wskazuje, gdzie dane mają trafić.
    int D = Grid->D;
    int G = Grid->G;

    for (int cz = 1; cz <= D; cz++) {
        for (int cy = 1; cy <= D; cy++) {
            for (int cx = 1; cx <= D; cx++) {
                cell3 c = {cx, cy, cz};
                int cid = cellIndex(c, G);
                const CellAtomData* cell = &Grid->Cells[cid];

                for (int k = 0; k < cell->Natoms; k++) {
                    int i = Grid->GlobalIndex[cid * MAX_ATOMS_PER_CELL + k];

                    if (i < 0 || i >= Natoms) {
                        printf("Error: invalid index i=%d in cell (%d,%d,%d), cid=%d, slot %d\n",
                               i, cx, cy, cz, cid, k);
                        exit(1);
                        exit(1);
                    }

                    Argon.X[i]  = cell->X[k];
                    Argon.Y[i]  = cell->Y[k];
                 
                    Argon.Z[i]  = cell->Z[k];

                    Argon.Vx[i] = cell->Vx[k];
                    Argon.Vy[i] = cell->Vy[k];
                    Argon.Vz[i] = cell->Vz[k];

                    Argon.X0[i] = cell->X0[k];
                    Argon.Y0[i] = cell->Y0[k];
                    Argon.Z0[i] = cell->Z0[k];
                }
            }
        }
    }
}

// GatherPositionsFromCells()
// przeniesienie tylko pozycji atomów z siatki do AtomData
void GatherPositionsFromCells(AtomData& Argon, int Natoms, const GridCells* Grid) {
    // Zbiera tylko pozycje z komórek realnych do struktury AtomData.
    // Uywane do śledzenia ewolucji układu.
    int D = Grid->D;
    int G = Grid->G;

    for (int cz = 1; cz <= D; cz++) {
        for (int cy = 1; cy <= D; cy++) {
            for (int cx = 1; cx <= D; cx++) {
                cell3 c = {cx, cy, cz};
                int cid = cellIndex(c, G);
                const CellAtomData* cell = &Grid->Cells[cid];

                for (int k = 0; k < cell->Natoms; k++) {
                    int i = Grid->GlobalIndex[cid * MAX_ATOMS_PER_CELL + k];

                    if (i < 0 || i >= Natoms) {
                        printf("Error: invalid index i=%d in cell (%d,%d,%d), cid=%d, slot %d\n",
                               i, cx, cy, cz, cid, k);
                        exit(1);
                        exit(1);
                    }

                    Argon.X[i]  = cell->X[k];
                    Argon.Y[i]  = cell->Y[k];
                    Argon.Z[i]  = cell->Z[k];

                }
            }
        }
    }
}

// GatherVelocitiesFromCells()
// przeniesienie tylko prędkości atomów z siatki do AtomData
void GatherVelocitiesFromCells(AtomData& Argon, int Natoms, const GridCells* Grid) {
    /// Zbiera tylko prędkości z komórek realnych do struktury AtomData.
    // Na razie nie używane.
    int D = Grid->D;
    int G = Grid->G;

    for (int cz = 1; cz <= D; cz++) {
        for (int cy = 1; cy <= D; cy++) {
            for (int cx = 1; cx <= D; cx++) {
                cell3 c = {cx, cy, cz};
                int cid = cellIndex(c, G);
                const CellAtomData* cell = &Grid->Cells[cid];

                for (int k = 0; k < cell->Natoms; k++) {
                    int i = Grid->GlobalIndex[cid * MAX_ATOMS_PER_CELL + k];

                    if (i < 0 || i >= Natoms) {
                        printf("Error: GatherVelocitiesFromCells – invalid index i=%d in cell (%d,%d,%d) cid=%d, slot=%d\n",
       i, cx, cy, cz, cid, k);
                    }

                    Argon.Vx[i] = cell->Vx[k];
                    Argon.Vy[i] = cell->Vy[k];
                    Argon.Vz[i] = cell->Vz[k];
                }
            }
        }
    }
}

// MaxDisplacement_LC()
// obliczenie maksymalnego przemieszczenia w modelu LC

// UpdateHalfVelocity_LC()
// aktualizacja prędkości o połowę kroku w modelu LC
void UpdateHalfVelocity_LC(GridCells& Grid, real dt) {
    int D = Grid.D;
    int G = Grid.G;
    real dt_half = 0.5*dt;

    for (int host_z = 1; host_z <= D; host_z++) {
        for (int host_y = 1; host_y <= D; host_y++) {
            for (int host_x = 1; host_x <= D; host_x++) {
                cell3 host = {host_x, host_y, host_z};
                int cid = cellIndex(host, G);
                CellAtomData& cell = Grid.Cells[cid];

                for (int i = 0; i < cell.Natoms; ++i) {
                    cell.Vx[i] += dt_half * cell.Fx[i] / m;
                    cell.Vy[i] += dt_half * cell.Fy[i] / m;
                    cell.Vz[i] += dt_half * cell.Fz[i] / m;
                }
            }
        }
    }
}

// UpdatePositions_LC()
// aktualizacja położeń atomów w modelu LC i sprawdzanie przemieszczenia
DisplacementStats UpdatePositions_LC(GridCells& Grid, real dt, parameters Par) {
    int D = Grid.D;
    int G = Grid.G;

    real sum_disp2 = 0.0f;
    real max_disp2 = 0.0f;
    //bool rebuild_needed = false;

    const real threshold2 = margin_coef * Par.Buffer * margin_coef * Par.Buffer;  // (margin_coef  * Buffer)^2
    int total_atoms = 0;

    for (int host_z = 1; host_z <= D; ++host_z) {
        for (int host_y = 1; host_y <= D; ++host_y) {
            for (int host_x = 1; host_x <= D; ++host_x) {
                cell3 host = {host_x, host_y, host_z};
                int cid = cellIndex(host, G);
                if (cid < 0 || cid >= Grid.Ncells) {
                    printf("ERROR: bad cid=%d from (%d %d %d) at G=%d\n", cid, host.x, host.y, host.z, G);
                    exit(1);
                }
                //printf("cid=%d from  (%d %d %d) \n",cid, host.x, host.y, host.z, G);
                CellAtomData& cell = Grid.Cells[cid];
                total_atoms += cell.Natoms;
                //printf("Natoms %d total atoms %d\n",cell.Natoms,total_atoms);

                for (int i = 0; i < cell.Natoms; ++i) {
                    // Propagacja pozycji
                    cell.X[i] += dt * cell.Vx[i];
                    cell.Y[i] += dt * cell.Vy[i];
                    cell.Z[i] += dt * cell.Vz[i];

                    // Różnica względem pozycji referencyjnej
                    real dx0 = cell.X[i] - cell.X0[i];
                    real dy0 = cell.Y[i] - cell.Y0[i];
                    real dz0 = cell.Z[i] - cell.Z0[i];
                    real disp2 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;

                    sum_disp2 += disp2;
                    if (disp2 > max_disp2) max_disp2 = disp2;
                }
            }
        }
    }
    //printf("total_atoms %d ",total_atoms);
    DisplacementStats stats;
    stats.avg_disp = (total_atoms > 0) ? sqrt(sum_disp2 / total_atoms) : 0.0f;
    //printf("Avg disp %f ",stats.avg_disp);
    stats.max_disp = sqrt(max_disp2);
     //printf("Max disp %f ",stats.max_disp);
    stats.rebuild_needed = (max_disp2 > threshold2);
    //printf("Rebuild? %d\n",(int)stats.rebuild_needed );
    return stats;
}

// UpdateKinetic_LC()
// obliczanie energii kinetycznej w modelu LC
real UpdateKinetic_LC(const GridCells& Grid) {

    int D = Grid.D;
    int G = Grid.G;
    real Ekin = 0.0;

    for (int host_z = 1; host_z <= D; host_z++) {
        for (int host_y = 1; host_y <= D; host_y++) {
            for (int host_x = 1; host_x <= D; host_x++) {
                cell3 host = {host_x, host_y, host_z};
                int cid = cellIndex(host, G);
                CellAtomData& cell = Grid.Cells[cid];

                for (int i = 0; i < cell.Natoms; ++i) {
                    real vx = cell.Vx[i];
                    real vy = cell.Vy[i];
                    real vz = cell.Vz[i];
                    Ekin += 0.5 * m * (vx * vx + vy * vy + vz * vz);
                }
            }
        }
    }
    return Ekin;
}

// Verlet_LC()
// integracja ruchu w modelu LC
energy3 Verlet_LC(int Natoms, AtomData& Argon, GridCells& Grid, real dt, parameters Par, elj3* OutEpot) {
    energy3 locE = {};
    elj3 locELJ = {} ;

    T_verlet.start();
    // 1. Pierwsza połowa kroku prędkości (v ← v + 0.5 * a * dt)
    UpdateHalfVelocity_LC(Grid, dt);

    // 2. Pozycje (r ← r + v * dt)
    DisplacementStats stats = UpdatePositions_LC(Grid, dt, Par);

    // 3. Odświeżenie ghostów (po aktualizacji realnych pozycji)
    T_ghost.start();
    FillGhostCells(&Grid);
    T_ghost.add();

    // 4. Sprawdzenie przemieszczeń względem X0 i ewentualna rekonstrukcja
    //real max_disp = MaxDisplacement_LC(Grid);
    // if (max_disp > Par.Buffer*0.5) {
    if (stats.rebuild_needed) {
        //printf("Rebuilding\n");
         T_rebuild.start();
        // a) Zrzucenie danych z komórek do globalnej struktury
        GatherAtomsFromCells(Argon, Natoms, &Grid);
        // b) Ponowne przypisanie do komórek, ustawia nowe X0
        AssignAtomsToCells(Argon, Natoms, &Grid);
        // c) Odświeżenie ghostów po nowym przypisaniu
        FillGhostCells(&Grid);
        T_rebuild.add();
    }
    
    // 5. Obliczenie sił i energii potencjalnej
    T_force.start();
    locELJ = ForcesEnergiesWrapper_LC(Grid,Par);
    T_force.add();
    *OutEpot = locELJ;

    // 6. Druga połowa kroku prędkości
    UpdateHalfVelocity_LC(Grid, dt);

    // 7. Energia kinetyczna i całkowita
    locE.Epot = locELJ.ELJ;
    locE.Ekin = UpdateKinetic_LC(Grid);
    locE.Etot = locE.Epot + locE.Ekin;
    printf("Energies (tot pot kin) %f %f %f \n", locE.Etot,locE.Epot,locE.Ekin);
     T_verlet.add();
    return locE;
}

// RunHeatingPhase_LC()
// faza podgrzewania w modelu LC
void RunHeatingPhase_LC(int Natoms, AtomData& locAtoms, GridCells& Grid, parameters Par, real Vmax, real t0) {
    // Faza "podgrzewania" układu: krótkie symulacje z kolejnymi re-inicjalizacjami prędkości.
    // Pozwala osiągnąć docelową temperaturę w układzie.
    auto start = TimerStart();
    constexpr real HeatingStep = 0.02f;
    const int HeatingSteps = int(HeatingTime / HeatingStep);
    elj3 locELJ = {} ;

    //    real Ekin_temp = Start(Natoms, locAtoms, Vmax * sqrtf(0.5f));
    real Ekin_temp =  Restart(locAtoms, Natoms, Grid, Vmax *sqrt(0.5f));
    printf("Thermalization Start / %d at T ≈ %.2f K (is %f)\n", N_Heating, t0 * sqrt(0.5f),Ekin_temp);

    locELJ = ForcesEnergies_LC(Grid,Par);
    for (int step_ = 0; step_ < HeatingSteps; step_++) {
        elj3 dummyEpot = {};
        
        energy3 Etmp = Verlet_LC(Natoms, locAtoms, Grid, HeatingStep, Par, &dummyEpot);
        real Tinst = 2 * Etmp.Ekin / ((3 * Natoms - 3) * kB);
        if (step_ == HeatingSteps - 1) {
            printf("  -> Final temperature for starting cycle: %.2f K\n", Tinst);
        }
    }
    for (int cycle = 0; cycle < N_Heating; cycle++) {
        printf("Thermalization cycle %d / %d at T ≈ %.2f K\n", cycle + 1, N_Heating, t0);
        real Ekin_temp =  Restart(locAtoms, Natoms, Grid, Vmax *sqrt(0.5f));
        locELJ = ForcesEnergies_LC(Grid,Par);
        for (int step_ = 0; step_ < HeatingSteps; step_++) {
            elj3 dummyEpot = {};
            energy3 Etmp = Verlet_LC(Natoms, locAtoms, Grid, HeatingStep, Par, &dummyEpot);
            real Tinst = 2 * Etmp.Ekin / ((3 * Natoms - 3) * kB);
            if (step_ == HeatingSteps - 1) {
                printf("  -> Final temperature in cycle %d: %.2f K\n", cycle + 1, Tinst);
            }
        }
    }
    PrintPhaseTime("Heating", TimerStop(start));
}

// RunThermalisationPhase_LC()
// faza termalizacji w modelu LC
void RunThermalisationPhase_LC(int Natoms, AtomData& Atoms, GridCells& Grid, parameters Par, real dt, int Term_, int save) {
    // Faza pomocniacza
    // Pozwala osiągnąć równomierny rozkład energii kinetycznej i temperatury w układzie.
    
    auto start = TimerStart();
    FILE *CSV = fopen("therm_energies.csv", "w");
    FILE *DISP = fopen("therm_disp.csv", "w");
    elj3 locELJ = {} ;
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
        locE = Verlet_LC(Natoms, Atoms, Grid, dt, Par, &Epot);
        if (step == 0) Estart = locE.Etot;
        if ((step % save) == 0) {
            real avg_disp = 0.0f, max_disp = 0.0f;
            T_stats.start();
            GatherPositionsFromCells(Atoms, Natoms, &Grid);
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
    fprintf(CSV, "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f\n",Term_ * dt, locE.Etot, locE.Epot, locE.Ekin, Epot.ELJ12, Epot.ELJ6, Temp);
    printf("Step %i, ETOT = %lf \n", Term_, locE.Etot);
    printf("Etot = %f +/- %f, Emin=%f  Emax=%f, DE=%f ", Eavg, Estdev, Emin, Emax, Emax - Emin);
    printf("Estart=%f, Eend=%f, DELTA=%f \n", Estart, Eend, Estart - Eend);
    fclose(CSV);
    fclose(DISP);
    PrintPhaseTime("Thermalisation", TimerStop(start));
}

// RunProductionPhase_LC()
// faza produkcyjna w modelu LC
void RunProductionPhase_LC(int Natoms, AtomData& Atoms, GridCells& Grid, parameters Par, real dt, int Term_, int save) {
    // Główna faza symulacji. 
    // Wykonują kolejne kroki MD i zapisują dane do plików CSV.
    // Wersja LC opiera się na siatce komórek i uwzględnia rekonstrukcje siatki.
    // Zbiera dane do obliczeń statystycznych: energia, temperatura, przesunięcia.
    auto start = TimerStart();
    FILE *CSV = fopen(Fname, "w");
    FILE *DISP = fopen("displacement.csv", "w");
    fprintf(CSV, "  Time,     Etot,      Epot,      Ekin,     Elj12,      Elj6,    Temp\n");
    fprintf(DISP, "  Time, AvgDisp, MaxDisp\n");

    double Eavg = 0.0, Evar = 0.0, Estdev = 0.0, Exx = 0.0, Ex = 0.0;
    double Emin = 1e20, Emax = -1e20, Estart = 0.0, Eend = 0.0;
    energy3 locE = {} ;
    elj3 Epot = {};
    real Temp;
    elj3 locELJ = {} ;
    printf("\n\nProduction phase for %d steps\n\n",Term_);
    for (int i = 0; i < Natoms; i++) {
        Atoms.X0[i] = Atoms.X[i];
        Atoms.Y0[i] = Atoms.Y[i];
        Atoms.Z0[i] = Atoms.Z[i];
    }
    for (int step = 0; step < Term_; step++) {
        locE = Verlet_LC(Natoms, Atoms, Grid, dt, Par, &Epot);
        if (step == 0) Estart = locE.Etot;
        if ((step % save) == 0) {
            real avg_disp = 0.0f, max_disp = 0.0f;
            T_stats.start();
            GatherPositionsFromCells(Atoms, Natoms, &Grid);
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
// Simulation_LC()
// pełna procedura w modelu LC
void Simulation_LC(int SIZ, real dt, int Steps, int save, AtomData& Atoms, GridCells& Grid, parameters Par, real t0, bool doThermalize) {
    int Natoms = Atoms.Natoms;
    int therm_steps = 100;
    real therm_time_step = 0.01;

    real Vmax = sqrt(3 * kB * t0 / m);
    T_total.start();
    RunHeatingPhase_LC(Natoms, Atoms, Grid, Par, Vmax, t0);
    RunThermalisationPhase_LC(Natoms, Atoms, Grid, Par, dt, therm_steps, save);
    RunProductionPhase_LC(Natoms, Atoms, Grid, Par, dt, Steps, save);
    T_total.add();
    FreeAll(Atoms);

}

// Restart()
// przepisanie danych między siatką a AtomData i ponowna inicjalizacja
real Restart(AtomData& Argon, int Natoms, GridCells& Grid, real Vmax) {
    // Reinicjalizuje prędkości atomów i przypisuje je ponownie do komórek.
    // Używane przed rozpoczęciem symulacji lub nowej fazy termalizacji.
    // Daje możliwość nadania losowych prędkości zgodnie z określoną temperaturą (Vmax).

    // Przeniesienie danych między strukturami
    GatherPositionsFromCells(Argon, Natoms, &Grid);
    // Nadanie nowych prędkości
    real Ekin = Start(Natoms, Argon, Vmax);
    // Przeniesienie danych między strukturami
    AssignAtomsToCells(Argon, Natoms, &Grid);
    //printf("Restart(): T ≈ %.2f K (z Ekin = %f)\n", 2 * Ekin / ((3 * Natoms - 3) * kB), Ekin);

    return Ekin;
}