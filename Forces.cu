#include "Types.h"
#include "Common.h" 

#ifdef USE_GPU
#include "GPU.h"
#endif

#include "Argon.h"
//#include "Forces.h"
#include "Argon_LC.h"
#include "GPU.h"

#include <stdio.h>
#include <cuda_runtime.h>

/*
elj3 ForcesEnergies_GPU(int Natoms, AtomData& Atoms, AtomDataGPU& gpuAtoms, AtomDataGPU& cpuAtoms, parameters Pars) {
    // Stałe Lennard-Jonesa
    real sig2 = sigma * sigma;
    real sig6 = sig2 * sig2 * sig2;
    real sig12 = sig6 * sig6;

    // Zerowanie sił na CPU
    for (int i = 0; i < Natoms; i++) {
        Atoms.Fx[i] = Atoms.Fy[i] = Atoms.Fz[i] = 0.0;
    }

    // Przesyłanie położeń na GPU
    CALL_CUDA(cudaMemcpy(gpuAtoms.X, Atoms.X, Natoms * sizeof(real), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(gpuAtoms.Y, Atoms.Y, Natoms * sizeof(real), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(gpuAtoms.Z, Atoms.Z, Natoms * sizeof(real), cudaMemcpyHostToDevice));

    // Konfiguracja kernela
    dim3 threadsPerBlock(128, 1, 1);
    dim3 numBlocks((Natoms + 127) / 128, 1, 1);

#ifdef FORCE_GPU
    ForcesEnergiesKernel_naive<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12);
    CALL_CUDA(cudaDeviceSynchronize());

#elif defined(FORCE_GPU_TILE)
    ForcesEnergiesKernel_tiles<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12, gpuAtoms.ELJ);  // przykład
    CALL_CUDA(cudaDeviceSynchronize());

#elif defined(FORCE_GPU_SHARED)
    ForcesEnergiesKernel_shared<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12, nullptr, nullptr);  // bufory TBD
    CALL_CUDA(cudaDeviceSynchronize());

#elif defined(FORCE_GPU_NEWTON)
    ForcesEnergiesKernel_newton<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12, nullptr, nullptr);  // bufory TBD
    CALL_CUDA(cudaDeviceSynchronize());

#else
    #error "Musisz zdefiniować jedną z opcji: FORCE_GPU, FORCE_GPU_TILE, FORCE_GPU_SHARED, FORCE_GPU_NEWTON"
#endif


    // Pobranie wyników z GPU
    CALL_CUDA(cudaMemcpy(cpuAtoms.Fx, gpuAtoms.Fx, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.Fy, gpuAtoms.Fy, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.Fz, gpuAtoms.Fz, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.ELJ6, gpuAtoms.ELJ6, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.ELJ12, gpuAtoms.ELJ12, Natoms * sizeof(real), cudaMemcpyDeviceToHost));

    // Sumowanie energii
    elj3 locE;
    locE.ELJ6 = locE.ELJ12 = locE.ELJ = 0.0f;

    for (int i = 0; i < Natoms; i++) {
        locE.ELJ6  += cpuAtoms.ELJ6[i];
        locE.ELJ12 += cpuAtoms.ELJ12[i];
    }

    locE.ELJ6  *= 2.0f * epsilon;
    locE.ELJ12 *= epsilon;
    locE.ELJ   = locE.ELJ6 + locE.ELJ12;

    return locE;
}

*/

elj3 ForcesEnergies_GPU(int Natoms, AtomData& Atoms, AtomDataGPU& gpuAtoms, AtomDataGPU& cpuAtoms, parameters Pars) {
    //msg("ForcesEnergies_GPU Entered");

    
    // Stałe Lennard-Jonesa
    real sig2 = sigma * sigma;
    real sig6 = sig2 * sig2 * sig2;
    real sig12 = sig6 * sig6;

    // Zerowanie sił na CPU
    for (int i = 0; i < Natoms; i++) {
        Atoms.Fx[i] = Atoms.Fy[i] = Atoms.Fz[i] = 0.0;
    }


    //for (int i =0;i<10;i++) {
    //    printf("%d X=%f Y=%f Z=%f \n",i,Atoms.X[i],Atoms.Y[i],Atoms.Z[i]);
    //}

    // Przesyłanie położeń na GPU
    CALL_CUDA(cudaMemcpy(gpuAtoms.X, Atoms.X, Natoms * sizeof(real), cudaMemcpyHostToDevice));
    //msg("cudaMemcpy");
    CALL_CUDA(cudaMemcpy(gpuAtoms.Y, Atoms.Y, Natoms * sizeof(real), cudaMemcpyHostToDevice));
    //msg("cudaMemcpy");
    CALL_CUDA(cudaMemcpy(gpuAtoms.Z, Atoms.Z, Natoms * sizeof(real), cudaMemcpyHostToDevice));
    //msg("cudaMemcpy");

    for(int i=0;i<Natoms;i++) {
        Atoms.Fx[i] = cpuAtoms.Fx[i];
        Atoms.Fy[i] = cpuAtoms.Fy[i];
        Atoms.Fz[i] = cpuAtoms.Fz[i];
    }

    // Konfiguracja kernela
    dim3 threadsPerBlock(128, 1, 1);
    dim3 numBlocks((Natoms + 127) / 128, 1, 1);

#ifdef FORCE_GPU
    ForcesEnergiesKernel_naive<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12);
    CALL_CUDA(cudaDeviceSynchronize());

#elif defined(FORCE_GPU_TILE)
    ForcesEnergiesKernel_tiles<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12, gpuAtoms.ELJ);  // przykład
    CALL_CUDA(cudaDeviceSynchronize());

#elif defined(FORCE_GPU_SHARED)
    ForcesEnergiesKernel_shared<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12, nullptr, nullptr);  // bufory TBD
    CALL_CUDA(cudaDeviceSynchronize());

#elif defined(FORCE_GPU_NEWTON)
    ForcesEnergiesKernel_newton<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12, nullptr, nullptr);  // bufory TBD
    CALL_CUDA(cudaDeviceSynchronize());

#else
    #error "Musisz zdefiniować jedną z opcji: FORCE_GPU, FORCE_GPU_TILE, FORCE_GPU_SHARED, FORCE_GPU_NEWTON"
#endif

    //msg("Kernel done");
    // Pobranie wyników z GPU
    CALL_CUDA(cudaMemcpy(cpuAtoms.Fx, gpuAtoms.Fx, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.Fy, gpuAtoms.Fy, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.Fz, gpuAtoms.Fz, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.ELJ6, gpuAtoms.ELJ6, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.ELJ12, gpuAtoms.ELJ12, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    //msg("All cudaMemcpy");
    // Sumowanie energii
    elj3 locE;
    locE.ELJ6 = locE.ELJ12 = locE.ELJ = 0.0f;
    real SumFx=0; 

    for (int i = 0; i < Natoms; i++) {
        locE.ELJ6  += cpuAtoms.ELJ6[i];
        locE.ELJ12 += cpuAtoms.ELJ12[i];
        SumFx += cpuAtoms.Fx[i];
        Atoms.Fx[i]=cpuAtoms.Fx[i];
        Atoms.Fy[i]=cpuAtoms.Fy[i];
        Atoms.Fz[i]=cpuAtoms.Fz[i];
    }

    locE.ELJ6  *= 2.0f;
    //locE.ELJ12 *= epsilon;
    locE.ELJ   = locE.ELJ6 + locE.ELJ12;

    //printf("Suma ELJ12=%f ELJ6=%f ELJ=%f Fx=%f\n",locE.ELJ12,locE.ELJ6,locE.ELJ,SumFx);

    return locE;
}



// Dummy wersja kernela: 1 wątek na atom
__global__ void ForcesEnergiesKernel_naive_dummy(int Natoms, AtomDataGPU atoms, float BoxSize, float Rc, float sig6, float sig12) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx == 0) {
        printf("[GPU] Wywołano ForcesEnergiesKernel_naive_dummy\n");
        asm("trap;");
    }
}

// Dummy wersja kernela tiled (bez shared memory)
__global__ void ForcesEnergiesKernel_tiles_dummy(int Natoms, AtomDataGPU atoms, float BoxSize, float Rc, float sig6, float sig12, float* gbuf) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx == 0) {
        printf("[GPU] Wywołano ForcesEnergiesKernel_tiles_dummy\n");
        asm("trap;");
    }
}

// Dummy wersja kernela z shared memory
__global__ void ForcesEnergiesKernel_shared_dummy(int Natoms, AtomDataGPU atoms, float BoxSize, float Rc, float sig6, float sig12, float* host_buf, float* guest_buf) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx == 0) {
        printf("[GPU] Wywołano ForcesEnergiesKernel_shared_dummy\n");
        asm("trap;");
    }
}

// Dummy wersja kernela z zasadą Newtona (akcja–reakcja)
__global__ void ForcesEnergiesKernel_newton_dummy(int Natoms, AtomDataGPU atoms, float BoxSize, float Rc, float sig6, float sig12, float* host_buf, float* guest_buf) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if (idx == 0) {
        printf("[GPU] Wywołano ForcesEnergiesKernel_newton_dummy\n");
        asm("trap;");
    }
}

__global__ void ForcesEnergiesKernel_naive(int Natoms, AtomDataGPU atoms, real BoxSize, real Rc,real sig6, real sig12){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    real r2, r6, r12;
    real DF;
    real DX, DY, DZ;
    //real X, Y, Z;
    real Fx, Fy, Fz;
    //real sig2, sig4, sig6, sig12;
    //real DFX,DFY,DFZ;
    real DELJ6, DELJ12;
    real ELJ6, ELJ12;
    real Rcut2=Rc*Rc;

    //if (idx < 3) {
    //    printf("%d X=%f Y=%f Z=%f \n",idx,atoms.X[idx],atoms.Y[idx],atoms.Z[idx]);
    //}

    if (idx <Natoms) {
        Fx = Fy = Fz =0.0;
        ELJ6 = ELJ12 = 0.0;
        for (int j=0;j<Natoms;j++) if (j != idx) {
            DX=atoms.X[idx]-atoms.X[j];
            DY=atoms.Y[idx]-atoms.Y[j];
            DZ=atoms.Z[idx]-atoms.Z[j];

            #ifdef USE_PBC
                DX -= Box * roundf(DX / Box);
                DY -= Box * roundf(DY / Box);
                DZ -= Box * roundf(DZ / Box);
            #endif

            r2 = DX*DX+DY*DY+DZ*DZ;
            if (r2<Rcut2) {
                r6 = r2*r2*r2;
                r12=r6*r6;
                DELJ6 = sig6/r6;
                DELJ12 = sig12/r12;
                ELJ6 -= DELJ6;
                ELJ12 += DELJ12;
                DF = 12*epsilon*(DELJ12-DELJ6);
                Fx += DF*DX/r2;
                Fy += DF*DY/r2;
                Fz += DF*DZ/r2;
            }
        }
        ELJ6 = ELJ6/2;
        ELJ12 = ELJ12/2;
        atoms.ELJ6[idx] = ELJ6;
        atoms.ELJ12[idx] =ELJ12;
        atoms.ELJ[idx] = epsilon*(ELJ12 - ELJ6);
        atoms.Fx[idx] = Fx;
        atoms.Fy[idx] = Fy;
        atoms.Fz[idx] = Fz;
    }
    //if (idx < 3) {
    //    printf("%d ELJ12=%f ELJ6=%f ELJ=%f \n",idx,atoms.ELJ12[idx],atoms.ELJ6[idx],atoms.ELJ[idx]);
    //}

}

__global__ void ForcesEnergiesKernel_tiles(int Natoms, AtomDataGPU atoms, real BoxSize, real Rc,real sig6, real sig12){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    /*
    real r2, r6, r12;
    real DF;
    real DX, DY, DZ;
    //real X, Y, Z;
    real Fx, Fy, Fz;
    //real sig2, sig4, sig6, sig12;
    //real DFX,DFY,DFZ;
    real DELJ6, DELJ12;
    real ELJ6, ELJ12;
    real Rcut2=Rc*Rc;
    */
    if (idx==0) {
        printf("DEBUG: entering trap in ForcesEnergiesKernel_tiles\n");
        asm("trap;");
    }

}

__global__ void ForcesEnergiesKernel_shared(int Natoms, AtomDataGPU atoms, real BoxSize, real Rc,real sig6, real sig12){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    /*
    real r2, r6, r12;
    real DF;
    real DX, DY, DZ;
    //real X, Y, Z;
    real Fx, Fy, Fz;
    //real sig2, sig4, sig6, sig12;
    //real DFX,DFY,DFZ;
    real DELJ6, DELJ12;
    real ELJ6, ELJ12;
    real Rcut2=Rc*Rc;
    */
    if (idx==0) {
        printf("DEBUG: entering trap in ForcesEnergiesKernel_shared\n");
        asm("trap;");
    }

}
__global__ void ForcesEnergiesKernel_newton(int Natoms, AtomDataGPU atoms, real BoxSize, real Rc,real sig6, real sig12){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    /*
    real r2, r6, r12;
    real DF;
    real DX, DY, DZ;
    //real X, Y, Z;
    real Fx, Fy, Fz;
    //real sig2, sig4, sig6, sig12;
    //real DFX,DFY,DFZ;
    real DELJ6, DELJ12;
    real ELJ6, ELJ12;
    real Rcut2=Rc*Rc;
    */
    if (idx==0) {
        printf("DEBUG: entering trap in ForcesEnergiesKernel_newton\n");
        asm("trap;");
    }

}

\

// Nowy wrapper

void ForcesEnergies_GPU_Simple(int Natoms, AtomData& Atoms, AtomDataGPU& gpuAtoms, AtomDataGPU& cpuAtoms, parameters Pars) {
    //msg("ForcesEnergies_GPU_Simpe Entered");

    // Stałe Lennard-Jonesa
    real sig2 = sigma * sigma;
    real sig6 = sig2 * sig2 * sig2;
    real sig12 = sig6 * sig6;
 
    // Konfiguracja kernela
    dim3 threadsPerBlock(128, 1, 1);
    dim3 numBlocks((Natoms + 127) / 128, 1, 1);

#ifdef FORCE_GPU
    ForcesEnergiesKernel_naive<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12);
    CALL_CUDA(cudaDeviceSynchronize());

#elif defined(FORCE_GPU_TILE)
    ForcesEnergiesKernel_tiles<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12, gpuAtoms.ELJ);  // przykład
    CALL_CUDA(cudaDeviceSynchronize());

#elif defined(FORCE_GPU_SHARED)
    ForcesEnergiesKernel_shared<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12, nullptr, nullptr);  // bufory TBD
    CALL_CUDA(cudaDeviceSynchronize());

#elif defined(FORCE_GPU_NEWTON)
    ForcesEnergiesKernel_newton<<<numBlocks, threadsPerBlock>>>(
        Natoms, gpuAtoms, Pars.BoxSize, Pars.Rc, sig6, sig12, nullptr, nullptr);  // bufory TBD
    CALL_CUDA(cudaDeviceSynchronize());

#else
    #error "Musisz zdefiniować jedną z opcji: FORCE_GPU, FORCE_GPU_TILE, FORCE_GPU_SHARED, FORCE_GPU_NEWTON"
#endif

    //return locE;
}

elj3 ForcesEnergies_GPU2(int Natoms, AtomData& Atoms, AtomDataGPU& gpuAtoms, AtomDataGPU& cpuAtoms, parameters Pars) {
    //msg("ForcesEnergies_GPU Entered");

    


    // Zerowanie sił na CPU
    for (int i = 0; i < Natoms; i++) {
        Atoms.Fx[i] = Atoms.Fy[i] = Atoms.Fz[i] = 0.0;
    }

    // Przesyłanie położeń na GPU
    CALL_CUDA(cudaMemcpy(gpuAtoms.X, Atoms.X, Natoms * sizeof(real), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(gpuAtoms.Y, Atoms.Y, Natoms * sizeof(real), cudaMemcpyHostToDevice));
    CALL_CUDA(cudaMemcpy(gpuAtoms.Z, Atoms.Z, Natoms * sizeof(real), cudaMemcpyHostToDevice));

    for(int i=0;i<Natoms;i++) {
        Atoms.Fx[i] = cpuAtoms.Fx[i];
        Atoms.Fy[i] = cpuAtoms.Fy[i];
        Atoms.Fz[i] = cpuAtoms.Fz[i];
    }

    #ifdef GPU_SIMPLE
    //msg("Kernel done");
       ForcesEnergies_GPU_Simple(Natoms, Atoms, gpuAtoms, cpuAtoms, Pars);
    #endif 
    // Pobranie wyników z GPU
    CALL_CUDA(cudaMemcpy(cpuAtoms.Fx, gpuAtoms.Fx, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.Fy, gpuAtoms.Fy, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.Fz, gpuAtoms.Fz, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.ELJ6, gpuAtoms.ELJ6, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    CALL_CUDA(cudaMemcpy(cpuAtoms.ELJ12, gpuAtoms.ELJ12, Natoms * sizeof(real), cudaMemcpyDeviceToHost));
    //msg("All cudaMemcpy");
    // Sumowanie energii
 
    elj3 locE;
    locE.ELJ6 = locE.ELJ12 = locE.ELJ = 0.0f;
    real SumFx=0; 

    for (int i = 0; i < Natoms; i++) {
        locE.ELJ6  += cpuAtoms.ELJ6[i];
        locE.ELJ12 += cpuAtoms.ELJ12[i];
        SumFx += cpuAtoms.Fx[i];
        Atoms.Fx[i]=cpuAtoms.Fx[i];
        Atoms.Fy[i]=cpuAtoms.Fy[i];
        Atoms.Fz[i]=cpuAtoms.Fz[i];
    }

    locE.ELJ6  *= 2.0f;
    //locE.ELJ12 *= epsilon;
    locE.ELJ   = locE.ELJ6 + locE.ELJ12;

    return locE;
}

