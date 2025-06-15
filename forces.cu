#include <cuda_runtime.h>
#include "forces.h"
#include "parameters.h"
#include "stdio.h"

#define TILE 256
#define FULL_MASK 0xffffffff

__device__ float minimum_image(float cordi, float cell_length) {
    float t = 0.5f * cell_length;
    return cordi + cell_length * ((cordi <= -t) - (cordi > t));
}

__constant__ float device_L;

__global__ void forces_kernel(const float* rx, const float* ry, const float* rz,
                              float* fx, float* fy, float* fz,
                              float* epot, float* pres_vir) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    
    int tid = threadIdx.y * blockDim.x + threadIdx.x;
    int lid  = tid % warpSize;

    int i = idx;
    int j_init = idy * TILE;
       
    float rcut2 = RCUT * RCUT;

    float xi = rx[i];
    float yi = ry[i];
    float zi = rz[i];

    float fxi = 0.0f;
    float fyi = 0.0f;
    float fzi = 0.0f;

    float epoti = 0.0f;
    float pres_viri = 0.0f;

    __shared__ float sh_fx[TILE];
    __shared__ float sh_fy[TILE];
    __shared__ float sh_fz[TILE];

    __shared__ float sh_epot;
    __shared__ float sh_pres_vir;

    if (tid == 0) {
        for (int k = 0; k < TILE; k++) {
            sh_fx[k] = 0.0f;
            sh_fy[k] = 0.0f;
            sh_fz[k] = 0.0f;
        }
        
        sh_epot = 0.0f;
        sh_pres_vir = 0.0f;
    }

    if (i < N && i < j_init + TILE - 1) {
        int j = j_init > i ? 0 : i + 1 - j_init;
        for (; j < TILE && j < N - j_init; j++) {

            float xj = rx[j_init + j];
            float yj = ry[j_init + j];
            float zj = rz[j_init + j];
        
            float xij = minimum_image(xi - xj, device_L);
            float yij = minimum_image(yi - yj, device_L);
            float zij = minimum_image(zi - zj, device_L);

            float rij2 = xij * xij + yij * yij + zij * zij;
            
            if (rij2 <= rcut2) {
                
                float r2inv = 1.0f / rij2;
                float r6inv = r2inv * r2inv * r2inv;
        
                float fr = 24.0f * r2inv * r6inv * (2.0f * r6inv - 1.0f);
                
                epoti += 4.0f * r6inv * (r6inv - 1.0f) - ECUT;
                pres_viri +=  fr * rij2;
                
                fxi += fr * xij;
                fyi += fr * yij;
                fzi += fr * zij;
                
                atomicAdd(sh_fx + j, -fr * xij);
                atomicAdd(sh_fy + j, -fr * yij);
                atomicAdd(sh_fz + j, -fr * zij);
            }
        }

        atomicAdd(fx + i, fxi);
        atomicAdd(fy + i, fyi);
        atomicAdd(fz + i, fzi);
    }

    for (int offset = 16; offset > 0; offset /= 2) {
       epoti += __shfl_down_sync(FULL_MASK, epoti, offset);
       pres_viri += __shfl_down_sync(FULL_MASK, pres_viri, offset);
    }
    
    __syncthreads();
    
    if (lid == 0) {
        atomicAdd(&sh_epot, epoti);
        atomicAdd(&sh_pres_vir, pres_viri);
    }
    
    __syncthreads();
    
    if (tid == 0) {
        for (int j = 0; j < TILE && j < N - j_init; j++) {
            atomicAdd(fx + j_init + j, sh_fx[j]);
            atomicAdd(fy + j_init + j, sh_fy[j]);
            atomicAdd(fz + j_init + j, sh_fz[j]);
            
        }

        atomicAdd(epot, sh_epot);
        atomicAdd(pres_vir, sh_pres_vir);
    }
}

extern "C"
void forces(float* rxyz, float* fxyz, float* epot, float* pres,
            const float* temp, const float rho, const float V, const float L)
{

    float pres_vir;
    
    float *device_fxyz, *device_rxyz, *device_epot, *device_pres_vir;

    cudaMalloc((void**)&device_fxyz, 3 * N * sizeof(float));
    cudaMalloc((void**)&device_rxyz, 3 * N * sizeof(float));
    cudaMalloc((void**)&device_epot, sizeof(float));
    cudaMalloc((void**)&device_pres_vir, sizeof(float));
    
    cudaMemcpyToSymbol(device_L, &L, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpy(device_rxyz, rxyz, 3 * N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemset(device_fxyz, 0, 3 * N * sizeof(float));
    cudaMemset(device_epot, 0, sizeof(float));
    cudaMemset(device_pres_vir, 0, sizeof(float));
    
    
    int nthreads = 256;
    dim3 gridDim((N + nthreads - 1) / nthreads, (N + TILE - 1) / TILE);
    dim3 blockDim(nthreads, 1);

    forces_kernel<<<gridDim, blockDim>>>(device_rxyz, device_rxyz + N, device_rxyz + 2*N,
                                         device_fxyz, device_fxyz + N, device_fxyz + 2*N,
                                         device_epot, device_pres_vir);
    
    cudaDeviceSynchronize();
    cudaMemcpy(fxyz, device_fxyz, 3 * N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(epot, device_epot, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pres_vir, device_pres_vir, sizeof(float), cudaMemcpyDeviceToHost);

    pres_vir /= (V * 3.0f);
    *pres = *temp * rho + pres_vir;

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(err));
        abort();
    }
}