// ForcesEnergies()
// obliczanie sił i energii – prosta wersja, bez PBC

// ForcesEnergies_pbc()
// obliczanie sił i energii z periodycznymi warunkami brzegowymi

// ForcesEnergies_LC()
// obliczanie sił i energii w modelu linked-cell (bez Newtona)

// (opcjonalnie) inne wersje LC jeśli nie są przeniesione do GPU_LC.cu


#include "Types.h"
#include "Common.h" 

#ifdef USE_GPU
#include "GPU.h"
#endif

#include "Argon.h"
#include "Forces.h"
#include "Argon_LC.h"

#ifdef FORCE_CPU_PBC
#warning "FORCE_CPU_PBC aktywne w forces.cpp"
#else
#warning "FORCE_CPU_PBC NIEaktywne w forces.cpp"
#endif

elj3 ForcesEnergies(int Natoms, AtomData& Atoms, real Rc) {
    real r2, r4, r6, r12, DF, DX, DY, DZ;
    real sig2 = sigma * sigma;
    real sig4 = sig2 * sig2;
    real sig6 = sig2 * sig4;
    real sig12 = sig6 * sig6;
    real DFX, DFY, DFZ, DELJ6, DELJ12;
    real Rcut2 = Rc * Rc;

    elj3 locEpot = {0.0f, 0.0f, 0.0f};
    
    // Zerowanie sił przed obliczeniami 
    // – w pętli sumujem wkłady od par, ze względu na sposób realizacji:
    // od razu dodajemy do właściwej zmiennej,  musimy tę zmienną wyzerowac. 
    for (int i = 0; i < Natoms; i++) {
        Atoms.Fx[i] = Atoms.Fy[i] = Atoms.Fz[i] = 0.0f;
    }

    for (int i = 0; i < Natoms; i++) {
        for (int j = i + 1; j < Natoms; j++) {
            //  j > i  (prawo akcji i reakcji)

            DX = Atoms.X[i] - Atoms.X[j];
            DY = Atoms.Y[i] - Atoms.Y[j];
            DZ = Atoms.Z[i] - Atoms.Z[j];
            r2 = DX * DX + DY * DY + DZ * DZ;

            if (r2 < Rcut2) {
                r4 = r2 * r2;
                r6 = r2 * r4;
                r12 = r6 * r6;
                DELJ6 = sig6 / r6;
                DELJ12 = sig12 / r12;

                locEpot.ELJ6 -= DELJ6;
                locEpot.ELJ12 += DELJ12;

                DF = 12.0f * epsilon * (DELJ12 - DELJ6);
                DFX = DF * DX / r2;
                DFY = DF * DY / r2;
                DFZ = DF * DZ / r2;

                Atoms.Fx[i] += DFX;
                Atoms.Fx[j] -= DFX;
                Atoms.Fy[i] += DFY;
                Atoms.Fy[j] -= DFY;
                Atoms.Fz[i] += DFZ;
                Atoms.Fz[j] -= DFZ;
            }
        }
    }

    locEpot.ELJ6 *= 2 * epsilon;
    locEpot.ELJ12 *= epsilon;
    locEpot.ELJ = locEpot.ELJ6 + locEpot.ELJ12;

    return locEpot;
}

//elj3 ForcesEnergiesWrapper(int Natoms, AtomData& Atoms, parameters Par) {
//    return ForcesEnergies(Natoms, Atoms, Par.Rc);
//}

// Minimal image convention – periodyczne warunki brzegowe (PBC)
//
// W przestrzeni periodycznej każdy atom ma nieskończenie wiele obrazów
// rozmieszczonych co Box w każdej osi. Kiedy obliczamy siły między parą atomów,
// musimy wybrać ten obraz j-tego atomu, który jest *najbliżej* i-tego atomu.
//
// Innymi słowy, nie interesuje nas, gdzie *naprawdę* jest j-ty atom,
// tylko która jego kopia (obraz periodyczny) znajduje się najbliżej
// względem pozycji i-tego atomu.
//
// Korekta wektora odległości polega na „zawinięciu” różnicy współrzędnych,
// tak by długość była minimalna (czyli w zakresie [−Box/2, Box/2]).
//
// Formuła:
//     DX -= Box * roundf(DX / Box);
// gwarantuje właśnie tę minimalizację długości wektora odległości,
// a nie śledzenie trajektorii przez granice.
//
// Przykład:
//     Box = 10
//     X[i] = 1.0, X[j] = 9.8 → DX = -8.8 → round(-0.88) = -1 → DX += 10 → DX = 1.2
//     ⇒ najbliższy obraz j-tego atomu jest po lewej stronie, nie po prawej
//
// To *nie* oznacza, że cząstka przeszła przez granicę – to tylko matematyczny wybór
// obrazu j-tego atomu, który daje najkrótszy wektor r_ij.

/*
#ifndef FORCE_GPU
elj3 ForcesEnergiesWrapper(int Natoms, AtomData& Atoms, parameters Par) {
#ifdef FORCE_CPU
    return ForcesEnergies(Natoms, Atoms, Par.Rc);
#elif defined(FORCE_CPU_PBC)
    return ForcesEnergies_pbc(Natoms, Atoms, Par.Rc);
#else
    #error "Brak definicji CPU wrappera"
#endif
}
#endif
*/ 

elj3 ForcesEnergies_pbc(int Natoms, AtomData& Atoms, real Rc) {
    real r2, r4, r6, r12, DF, DX, DY, DZ;
    real sig2 = sigma * sigma;
    real sig4 = sig2 * sig2;
    real sig6 = sig2 * sig4;
    real sig12 = sig6 * sig6;
    real DFX, DFY, DFZ, DELJ6, DELJ12;
    real Rcut2 = Rc * Rc;
    // Rozmiar pudełka symulacyjnego – potrzebny do periodycznych warunków brzegowych
    real Box = Pars.BoxSize;
    
    elj3 locEpot = {0.0f, 0.0f, 0.0f};

    // Zerowanie sił
    for (int i = 0; i < Natoms; i++) {
        Atoms.Fx[i] = Atoms.Fy[i] = Atoms.Fz[i] = 0.0f;
    }

    // Pętla po wszystkich parach i < j
    for (int i = 0; i < Natoms; i++) {
        for (int j = i + 1; j < Natoms; j++) {
            DX = Atoms.X[i] - Atoms.X[j];
            DY = Atoms.Y[i] - Atoms.Y[j];
            DZ = Atoms.Z[i] - Atoms.Z[j];

            // Minimal image convention – patrzymy, która kopia j-tego atomu
            // jest najbliżej i-tego atomu w przestrzeni periodycznej.
            DX -= Box * roundf(DX / Box);
            DY -= Box * roundf(DY / Box);
            DZ -= Box * roundf(DZ / Box);

            //DX -= Box * (int)(DX / Box);
            //DY -= Box * (int)(DY / Box);
            //DZ -= Box * (int) (DZ / Box);
            r2 = DX * DX + DY * DY + DZ * DZ;

            if (r2 < Rcut2) {
                r4 = r2 * r2;
                r6 = r2 * r4;
                r12 = r6 * r6;
                DELJ6 = sig6 / r6;
                DELJ12 = sig12 / r12;

                locEpot.ELJ6 -= DELJ6;
                locEpot.ELJ12 += DELJ12;

                DF = 12.0f * epsilon * (DELJ12 - DELJ6);
                DFX = DF * DX / r2;
                DFY = DF * DY / r2;
                DFZ = DF * DZ / r2;

                Atoms.Fx[i] += DFX;
                Atoms.Fx[j] -= DFX;
                Atoms.Fy[i] += DFY;
                Atoms.Fy[j] -= DFY;
                Atoms.Fz[i] += DFZ;
                Atoms.Fz[j] -= DFZ;
            }
        }
    }

    locEpot.ELJ6 *= 2 * epsilon;
    locEpot.ELJ12 *= epsilon;
    locEpot.ELJ = locEpot.ELJ6 + locEpot.ELJ12;

    return locEpot;
}
/*
elj3 ForcesEnergies_LC(GridCells& Grid,parameters Par) {
    real Rc2 = Par.Rc*Par.Rc;
    real sig2 = sigma * sigma;
    real sig4 = sig2 * sig2;
    real sig6 = sig2 * sig4;
    real sig12 = sig6 * sig6;

    elj3 totalE = {0.0f, 0.0f, 0.0f};
    int D = Grid.D;
    int G = Grid.G;
    
    printf("Rc2 = %f\n",Rc2);
    // Obliczenia sił i energii w podejściu linked-cell, atomy w komórkach znajdujacych się na brzegu
    // oddziałują z atomami z komórek cieni/duchów (ghosts) wygenerowanymi przez rzutowanie 
    // obrazów komórek znadujących się po drugiej stronie pudła symulacyjnego. 
    // Iterujemy po wszystkich 27 komórkach wokół bieżącej (hosta), włącznie z nią samą.
    // Interakcje host-guest: każdy sąsiad

    /*
        Dla każdej komórki host (centralnej), przeglądamy wszystkie 26 sąsiednich komórek,
        tworząc sześcian 3×3×3 wokół hosta (czyli łącznie 27 komórek, w tym sam host).
                   
        W szczególności:
        - Jeśli guest == host, to liczymy wszystkie pary (i < j) w tej samej komórce,
            korzystając z zasady akcji-reakcji (obie cząstki dostają siłę).
        - Jeśli guest ≠ host, to liczymy oddziaływanie jednostronnie: tylko atomy hosta dostają siłę.
            Guest traktowany jest jako pasywne źródło potencjału (typowy model host/guest).


        Z punktu widzenia fizyki, układ jest zamknięty w pudełku z periodycznymi warunkami brzegowymi
        (PBC), które tworzą torus w 3D: jeśli cząstka wychodzi przez jedną ścianę pudełka,
        pojawia się po przeciwnej stronie. To oznacza, że potencjalna odległość między dwiema cząstkami
        musi być liczona względem ich najbliższych obrazów.

        Zamiast dynamicznego przeliczania odległości zgodnie z tzw. minimal image convention
        (czyli szukania obrazu atomu najbliższego drugiemu atomowi), stosujemy podejście z kopiami.
        Rozszerzamy siatkę komórek dodając do niej "cienie" — ghost cells — które zawierają
        kopie atomów z przeciwległych ścian. Dzięki temu wystarczy liczyć odległości wprost,
        bez korekt na PBC.

        Ten trick upraszcza znacznie kod i pozwala uzyskać regularną strukturę danych — bardzo istotną
        przy implementacji na GPU. Powtarzające się operacje i brak potrzeby dynamicznych korekt
        indeksów dają wysoką wydajność i prostotę debugowania.

        Dodatkowo jest to podejście szybsze: zamiast wykonywać dla każdej pary cząstek kosztowne
        obliczenia korekt periodycznych (co rośnie jak N²), wykonujemy raz liniowo (O(N)) kopiowanie
        cząstek do komórek cienia. Pozwala to uniknąć wielokrotnego dzielenia i zaokrąglania,
        co w klasycznym podejściu z minimal image convention pojawia się w każdej iteracji sił.
    

    int count_pairs=0;
    double avg_pairs;
    int sum_atoms=0;

    for (int host_z = 1; host_z <= D; host_z++) {
        for (int host_y = 1; host_y <= D; host_y++) {
            for (int host_x = 1; host_x <= D; host_x++) {
                cell3 host = {host_x, host_y, host_z};
                int host_id = cellIndex(host, G);
                CellAtomData* Host = &Grid.Cells[host_id];
                int Nhost = Host->Natoms;
                sum_atoms += Nhost;
                real* Fx_local = Host->Fx;
                real* Fy_local = Host->Fy;
                real* Fz_local = Host->Fz;

                // Zerowanie sił przed sumowaniem nowych wartości w tej komórce
                for (int i=0;i<Nhost;i++) {
                    Fx_local[i]=Fy_local[i]=Fz_local[i]=0.0;
                }
                elj3 Elocal = {0.0f, 0.0f, 0.0f};


                // Iteracja po sąsiadach (w tym po sobie)
                for (int dz = -1; dz <= 1; dz++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dx = -1; dx <= 1; dx++) {
                            cell3 guest = {host_x + dx, host_y + dy, host_z + dz};
                            int guest_id = cellIndex(guest, G);
                            CellAtomData* Guest = &Grid.Cells[guest_id];
                            int Nguest = Guest->Natoms;
                            
                            // Interakcje wewnątrz jednej komórki (host == guest)
                            if (host_id == guest_id) {
                                for (int i = 0; i < Nhost; i++) {
                                    //sum_atoms++;
                                    for (int j = i + 1; j < Nhost; j++) {
                                        real dx = Host->X[i] - Host->X[j];
                                        real dy = Host->Y[i] - Host->Y[j];
                                        real dz = Host->Z[i] - Host->Z[j];
                                        //count_pairs++;
                                        real r2 = dx * dx + dy * dy + dz * dz;
                                        if (r2 < Rc2) {
                                            count_pairs++;
                                            real r4 = r2 * r2;
                                            real r6 = r2 * r4;
                                            real r12 = r6 * r6;
                                            real DELJ6 = sig6 / r6;
                                            real DELJ12 = sig12 / r12;

                                            real DF = 12.0f * epsilon * (DELJ12 - DELJ6) / r2;
                                            real Fx = DF * dx;
                                            real Fy = DF * dy;
                                            real Fz = DF * dz;

                                            Fx_local[i] += Fx;
                                            Fy_local[i] += Fy;
                                            Fz_local[i] += Fz;

                                            Fx_local[j] -= Fx;
                                            Fy_local[j] -= Fy;
                                            Fz_local[j] -= Fz;

                                            Elocal.ELJ6  -= 2.0f * DELJ6;
                                            Elocal.ELJ12 += 2.0f * DELJ12;
                                        }
                                    }
                                }
                            } else {
                                // Interakcje między komórkami host/guest, tylko atomy hosta aktualizują siły, 
                                // guest traktowany jako pasywne źródło potencjału
                                for (int i = 0; i < Nhost; i++) {
                                    //sum_atoms++;
                                    for (int j = 0; j < Nguest; j++) {
                                        real dx = Host->X[i] - Guest->X[j];
                                        real dy = Host->Y[i] - Guest->Y[j];
                                        real dz = Host->Z[i] - Guest->Z[j];
                                        real r2 = dx * dx + dy * dy + dz * dz;
                                        if (r2 < Rc2) {
                                            count_pairs++;
                                            real r4 = r2 * r2;
                                            real r6 = r2 * r4;
                                            real r12 = r6 * r6;
                                            real DELJ6 = sig6 / r6;
                                            real DELJ12 = sig12 / r12;

                                            real DF = 12.0f * epsilon * (DELJ12 - DELJ6) / r2;
                                            real Fx = DF * dx;
                                            real Fy = DF * dy;
                                            real Fz = DF * dz;

                                            Fx_local[i] += Fx;
                                            Fy_local[i] += Fy;
                                            Fz_local[i] += Fz;

                                            Elocal.ELJ6  -= DELJ6;
                                            Elocal.ELJ12 += DELJ12;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                // Dodanie energii (jeszcze przed podziałem)
                totalE.ELJ6  += Elocal.ELJ6;
                totalE.ELJ12 += Elocal.ELJ12;
            }
        }
    }

    //printf("Total number of atoms %d pairs %d, average_pairs%f \n",sum_atoms,count_pairs,1.0*count_pairs/sum_atoms);
    // Finalizacja energii
    // Energię dzielimy przez 2, bo każda para była liczona dwukrotnie (symetrycznie)
    totalE.ELJ6  *= epsilon * 0.5f;
    totalE.ELJ12 *= epsilon * 0.5f;
    totalE.ELJ    = totalE.ELJ6 + totalE.ELJ12;
    //printf("totalE: (%f %f %f)\n",totalE.ELJ6,totalE.ELJ12,totalE.ELJ);
    return totalE;
}
*/
