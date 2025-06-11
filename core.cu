#include "core.h"
#include "parameters.h"

#include <math.h>
#include <stdlib.h> // rand()

#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifndef CUDA
    #include <immintrin.h>
    #include <omp.h>
#else 
    #include <cuda_runtime.h>
#endif


#define ECUT (float)(4.0 * (pow(RCUT, -12) - pow(RCUT, -6)))


void init_pos(float* rxyz,const float rho)
{
    // inicialización de las posiciones de los átomos en un cristal FCC

    float a = cbrt(4.0f / rho);
    int nucells = ceil(cbrt((double)N / 4.0));
    float *rx = rxyz, *ry = rxyz + N, *rz = rxyz + 2*N;
    int idx = 0;

    for (int i = 0; i < nucells; i++) {
        for (int j = 0; j < nucells; j++) {
            for (int k = 0; k < nucells; k++) {
                rx[idx] = i * a; // x
                ry[idx] = j * a; // y
                rz[idx] = k * a; // z
                    // del mismo átomo
                rx[idx + 1] = (i + 0.5f) * a;
                ry[idx + 1] = (j + 0.5f) * a;
                rz[idx + 1] = k * a;

                rx[idx + 2] = (i + 0.5f) * a;
                ry[idx + 2] = j * a;
                rz[idx + 2] = (k + 0.5f) * a;

                rx[idx + 3] = i * a;
                ry[idx + 3] = (j + 0.5f) * a;
                rz[idx + 3] = (k + 0.5f) * a;

                idx += 4;
            }
        }
    }
}


void init_vel(float* vxyz,float* temp,float* ekin)
{
    // inicialización de velocidades aleatorias

    float sf, sumvx = 0.0f, sumvy = 0.0f, sumvz = 0.0f, sumv2 = 0.0f;
    float *vx = vxyz, *vy = vxyz + N, *vz = vxyz + 2*N;

    for (int i = 0; i < N; i += 1) {
        vx[i] = rand() / (float)RAND_MAX - 0.5f;
        vy[i] = rand() / (float)RAND_MAX - 0.5f;
        vz[i] = rand() / (float)RAND_MAX - 0.5f;

        sumvx += vx[i];
        sumvy += vy[i];
        sumvz += vz[i];
        sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    }

    sumvx /= (float)N;
    sumvy /= (float)N;
    sumvz /= (float)N;
    *temp = sumv2 / (3.0f * N);
    *ekin = 0.5f * sumv2;
    sf = (float)sqrt(T0 / *temp);

    for (int i = 0; i < N; i += 1) { // elimina la velocidad del centro de masa
        // y ajusta la temperatura
        vx[i] = (vx[i] - sumvx) * sf;
        vy[i] = (vy[i] - sumvy) * sf;
        vz[i] = (vz[i] - sumvz) * sf;
    }
}

#ifndef CUDA
static float minimum_image(float cordi, const float cell_length) {
    float t = 0.5f * cell_length;
    return cordi + cell_length * ((cordi <= -t) - (cordi > t));
}

#ifdef SIMD_INTRINSICS
static __m256 minimum_image_v(__m256 cordi, const float cell_length)
{
    // imagen más cercana
/*
    if (cordi <= -0.5f * cell_length) {
        cordi += cell_length;
    } else if (cordi > 0.5f * cell_length) {
        cordi -= cell_length;
    }
    return cordi;
*/
    __m256 t = _mm256_set1_ps(0.5f);
    __m256 op_t = _mm256_set1_ps(-0.5f);
    __m256 L = _mm256_set1_ps(cell_length);
    __m256 mask_lower = _mm256_cmp_ps(cordi, _mm256_mul_ps(op_t, L), _CMP_LE_OQ);
    __m256 mask_higher = _mm256_cmp_ps(cordi, _mm256_mul_ps(t, L), _CMP_GT_OQ);
    cordi = _mm256_add_ps(cordi, _mm256_and_ps(L, mask_lower));
    return _mm256_sub_ps(cordi, _mm256_and_ps(L, mask_higher));
}

static float horizontal_sum(__m256 v) {
    v = _mm256_hadd_ps(v, v);
    v = _mm256_hadd_ps(v, v);
    v = _mm256_add_ps(v, _mm256_permute2f128_ps(v, v, 0x01));

    return _mm256_cvtss_f32(v);
}
#endif

#else
__device__ float minimum_image(float cordi, const float cell_length) {
    if (cordi <= -0.5f * cell_length) {
        cordi += cell_length;
    } else if (cordi > 0.5f * cell_length) {
        cordi -= cell_length;
    };
    return cordi;
}

#endif


#ifndef CUDA
void forces(float* rxyz, float* fxyz, float* epot, float* pres,
            const float* temp, const float rho, const float V, const float L)
{
    // calcula las fuerzas LJ (12-6)

    for (int i = 0; i < 3 * N; i++) {
        fxyz[i] = 0.0f;
    }

    float *rx = rxyz, *ry = rxyz + N, *rz = rxyz + 2 * N;
    float *fx = fxyz, *fy = fxyz + N, *fz = fxyz + 2 * N;

    float rcut2 = RCUT * RCUT;
    float _epot = 0.0f, pres_vir = 0.0f;

#ifndef SIMD_INTRINSICS

    for (int i = 0; i < (N - 1); i += 1) {

        float xi = rx[i];
        float yi = ry[i];
        float zi = rz[i];

        for (int j = i + 1; j < N; j += 1) {

            float xj = rx[j];
            float yj = ry[j];
            float zj = rz[j];

            // distancia mínima entre r_i y r_j
            float xij = xi - xj;
            xij = minimum_image(xij, L);
            float yij = yi - yj;
            yij = minimum_image(yij, L);
            float zij = zi - zj;
            zij = minimum_image(zij, L);

            float rij2 = xij * xij + yij * yij + zij * zij;

            if (rij2 <= rcut2) {
                float r2inv = 1.0f / rij2;
                float r6inv = r2inv * r2inv * r2inv;

                float fr = 24.0f * r2inv * r6inv * (2.0f * r6inv - 1.0f);

                fx[i] += fr * xij;
                fy[i] += fr * yij;
                fz[i] += fr * zij;

                fx[j] -= fr * xij;
                fy[j] -= fr * yij;
                fz[j] -= fr * zij;

                *epot += 4.0f * r6inv * (r6inv - 1.0f) - ECUT;
                pres_vir += fr * rij2;
            }
        }
    }

#else

    __m256 rcut2_v = _mm256_set1_ps(rcut2);
    __m256 ecut = _mm256_set1_ps(ECUT);
    __m256 epot_v = _mm256_setzero_ps();
    __m256 pres_v = _mm256_setzero_ps();

#ifndef TILING

    #pragma omp declare reduction(vadd : __m256 : omp_out = _mm256_add_ps(omp_in, omp_out))\
    initializer(omp_priv=_mm256_setzero_ps())

    #pragma omp parallel for default(shared) \
    reduction(vadd:epot_v) reduction(vadd:pres_v) \
    reduction(+:fx[:N]) reduction(+:fy[:N]) reduction(+:fz[:N]) \
    reduction(+:_epot) reduction(+:pres_vir) \
    schedule (dynamic, 64)
    for (unsigned int i = 0; i < N - 1; i++) {

        __m256 xi = _mm256_set1_ps(rx[i]);
        __m256 yi = _mm256_set1_ps(ry[i]);
        __m256 zi = _mm256_set1_ps(rz[i]);

        __m256 fxi = _mm256_setzero_ps();
        __m256 fyi = _mm256_setzero_ps();
        __m256 fzi = _mm256_setzero_ps();

        unsigned int j;

        for (j = i + 1; j < N - 8; j += 8) {

            __m256 xj = _mm256_loadu_ps(rx + j);
            __m256 yj = _mm256_loadu_ps(ry + j);
            __m256 zj = _mm256_loadu_ps(rz + j);

            // distancia mínima entre r_i y r_j
            __m256 xij = _mm256_sub_ps(xi, xj);
            // rx = rx + L * ((rx <= -t) - (rx > t));
            xij = minimum_image_v(xij, L);

            __m256 yij = _mm256_sub_ps(yi, yj);
            // rx = rx + L * ((rx <= -t) - (rx > t));
            yij = minimum_image_v(yij, L);

            __m256 zij = _mm256_sub_ps(zi, zj);
            // rx = rx + L * ((rx <= -t) - (rx > t));
            zij = minimum_image_v(zij, L);

            // float rij2 = rx * rx + ry * ry + rz * rz;
            __m256 rij2 = _mm256_add_ps(_mm256_mul_ps(xij, xij),
                                        _mm256_add_ps(_mm256_mul_ps(yij, yij),
                                                    _mm256_mul_ps(zij, zij)));

            // if (rij2 < rcut2) ...
            __m256 mask_rcut = _mm256_cmp_ps(rij2, rcut2_v, _CMP_LE_OQ);

	        if(!_mm256_movemask_ps(mask_rcut)) continue;

            // calculo de fr
            __m256 r2inv = _mm256_div_ps(_mm256_set1_ps(1.0f), rij2);
            __m256 r6inv = _mm256_mul_ps(r2inv, _mm256_mul_ps(r2inv, r2inv));
            __m256 fr = _mm256_mul_ps(_mm256_set1_ps(24.0f),
                                    _mm256_mul_ps(r2inv,
                                                    _mm256_mul_ps(r6inv,
                                                                _mm256_sub_ps(_mm256_mul_ps(_mm256_set1_ps(2.0f), r6inv),
                                                                                _mm256_set1_ps(1.0f)))));


            // actualizamos fxyz

            __m256 dfx = _mm256_mul_ps(fr, xij);
            __m256 dfy = _mm256_mul_ps(fr, yij);
            __m256 dfz = _mm256_mul_ps(fr, zij);

            // i
            
            fxi = _mm256_add_ps(fxi, _mm256_and_ps(dfx, mask_rcut));

            fyi = _mm256_add_ps(fyi, _mm256_and_ps(dfy, mask_rcut));

            fzi = _mm256_add_ps(fzi, _mm256_and_ps(dfz, mask_rcut));

            // j
            
            _mm256_storeu_ps(fx + j, _mm256_sub_ps(_mm256_loadu_ps(fx + j),
                                                _mm256_and_ps(dfx, mask_rcut)));

            _mm256_storeu_ps(fy + j, _mm256_sub_ps(_mm256_loadu_ps(fy + j),
                                                _mm256_and_ps(dfy, mask_rcut)));
            _mm256_storeu_ps(fz + j, _mm256_sub_ps(_mm256_loadu_ps(fz + j),
                                                _mm256_and_ps(dfz, mask_rcut)));


	    // actualizamos epot
            __m256 ep = _mm256_sub_ps(_mm256_mul_ps(_mm256_set1_ps(4.0f),
            _mm256_mul_ps(r6inv,
            _mm256_sub_ps(r6inv, _mm256_set1_ps(1.0f)))),
            ecut);
            epot_v = _mm256_add_ps(epot_v, _mm256_and_ps(ep, mask_rcut));
            
            // actualizamos pres
            __m256 pv = _mm256_mul_ps(fr, rij2);
            pres_v = _mm256_add_ps(pres_v, _mm256_and_ps(pv, mask_rcut));

        }

        fx[i] += horizontal_sum(fxi);
        fy[i] += horizontal_sum(fyi);
        fz[i] += horizontal_sum(fzi);
            
        // calculamos el resto (< 7 iteraciones)
        for (; j < N; j++) {

            // distancia mínima entre r_i y r_j
            float xij = rx[i] - rx[j];
            xij = minimum_image(xij, L);

            float yij = ry[i] - ry[j];
            yij = minimum_image(yij, L);

            float zij = rz[i] - rz[j];
            zij = minimum_image(zij, L);

            float rij2 = xij * xij + yij * yij + zij * zij;
            if (rij2 < rcut2) {
                float r2inv = 1.0f / rij2;
                float r6inv = r2inv * r2inv * r2inv;

                float fr = 24.0f * r2inv * r6inv * (2.0f * r6inv - 1.0f);

                fx[i] += fr * xij;
                fy[i] += fr * yij;
                fz[i] += fr * zij;

                fx[j] -= fr * xij;
                fy[j] -= fr * yij;
                fz[j] -= fr * zij;

                _epot += 4.0f * r6inv * (r6inv - 1.0f) - ECUT;
                pres_vir += fr * rij2;
            }
        }
        
    }
#else

    for (unsigned ii = 0; ii < N; ii += BLOCK) {

        unsigned int off = 1;
        for (unsigned jj = ii; jj < N; jj += BLOCK) {

            for (unsigned int i = ii; i < ii + BLOCK - off; i++) {

                __m256 xi = _mm256_set1_ps(rx[i]);
                __m256 yi = _mm256_set1_ps(ry[i]);
                __m256 zi = _mm256_set1_ps(rz[i]);

                __m256 fxi = _mm256_setzero_ps();
                __m256 fyi = _mm256_setzero_ps();
                __m256 fzi = _mm256_setzero_ps();

                unsigned int j;

                for (j = off*(i + 1 - jj) + jj; j + 8 < jj + BLOCK; j += 8) {

                    __m256 xj = _mm256_loadu_ps(rx + j);
                    __m256 yj = _mm256_loadu_ps(ry + j);
                    __m256 zj = _mm256_loadu_ps(rz + j);

                    // distancia mínima entre r_i y r_j
                    __m256 xij = _mm256_sub_ps(xi, xj);
                    // rx = rx + L * ((rx <= -t) - (rx > t));
                    xij = minimum_image_v(xij, L);

                    __m256 yij = _mm256_sub_ps(yi, yj);
                    // rx = rx + L * ((rx <= -t) - (rx > t));
                    yij = minimum_image_v(yij, L);

                    __m256 zij = _mm256_sub_ps(zi, zj);
                    // rx = rx + L * ((rx <= -t) - (rx > t));
                    zij = minimum_image_v(zij, L);

                    // float rij2 = rx * rx + ry * ry + rz * rz;
                    __m256 rij2 = _mm256_add_ps(_mm256_mul_ps(xij, xij),
                                                _mm256_add_ps(_mm256_mul_ps(yij, yij),
                                                            _mm256_mul_ps(zij, zij)));

                    // if (rij2 < rcut2) ...
                    __m256 mask_rcut = _mm256_cmp_ps(rij2, rcut2_v, _CMP_LE_OQ);

                    if(_mm256_testz_ps(mask_rcut, mask_rcut)) continue; // sexo

                    // calculo de fr
                    __m256 r2inv = _mm256_div_ps(_mm256_set1_ps(1.0f), rij2);
                    __m256 r6inv = _mm256_mul_ps(r2inv, _mm256_mul_ps(r2inv, r2inv));
                    __m256 fr = _mm256_mul_ps(_mm256_set1_ps(24.0f),
                                            _mm256_mul_ps(r2inv,
                                                            _mm256_mul_ps(r6inv,
                                                                        _mm256_sub_ps(_mm256_mul_ps(_mm256_set1_ps(2.0f), r6inv),
                                                                                        _mm256_set1_ps(1.0f)))));


                    // actualizamos fxyz

                    __m256 dfx = _mm256_mul_ps(fr, xij);
                    __m256 dfy = _mm256_mul_ps(fr, yij);
                    __m256 dfz = _mm256_mul_ps(fr, zij);

                    // i
                    
                    fxi = _mm256_add_ps(fxi, _mm256_and_ps(dfx, mask_rcut));

                    fyi = _mm256_add_ps(fyi, _mm256_and_ps(dfy, mask_rcut));

                    fzi = _mm256_add_ps(fzi, _mm256_and_ps(dfz, mask_rcut));


                    // j
                    _mm256_storeu_ps(fx + j, _mm256_sub_ps(_mm256_loadu_ps(fx + j),
                                                        _mm256_and_ps(dfx, mask_rcut)));

                    _mm256_storeu_ps(fy + j, _mm256_sub_ps(_mm256_loadu_ps(fy + j),
                                                        _mm256_and_ps(dfy, mask_rcut)));
                    _mm256_storeu_ps(fz + j, _mm256_sub_ps(_mm256_loadu_ps(fz + j),
                                                        _mm256_and_ps(dfz, mask_rcut)));

                    // actualizamos epot
                    __m256 ep = _mm256_sub_ps(_mm256_mul_ps(_mm256_set1_ps(4.0f),
                    _mm256_mul_ps(r6inv,
                    _mm256_sub_ps(r6inv, _mm256_set1_ps(1.0f)))),
                    ecut);
                    epot_v = _mm256_add_ps(epot_v, _mm256_and_ps(ep, mask_rcut));
                    
                    // actualizamos pres
                    __m256 pv = _mm256_mul_ps(fr, rij2);
                    pres_v = _mm256_add_ps(pres_v, _mm256_and_ps(pv, mask_rcut));

                }

                fx[i] += horizontal_sum(fxi);
                fy[i] += horizontal_sum(fyi);
                fz[i] += horizontal_sum(fzi);
                
                // calculamos el resto (< 7 iteraciones)
                for (; j < jj + BLOCK; j++) {

                    // distancia mínima entre r_i y r_j
                    float xij = rx[i] - rx[j];
                    xij = minimum_image(xij, L);

                    float yij = ry[i] - ry[j];
                    yij = minimum_image(yij, L);

                    float zij = rz[i] - rz[j];
                    zij = minimum_image(zij, L);

                    float rij2 = xij * xij + yij * yij + zij * zij;
                    if (rij2 < rcut2) {
                        float r2inv = 1.0f / rij2;
                        float r6inv = r2inv * r2inv * r2inv;

                        float fr = 24.0f * r2inv * r6inv * (2.0f * r6inv - 1.0f);

                        fx[i] += fr * xij;
                        fy[i] += fr * yij;
                        fz[i] += fr * zij;

                        fx[j] -= fr * xij;
                        fy[j] -= fr * yij;
                        fz[j] -= fr * zij;

                        *epot += 4.0f * r6inv * (r6inv - 1.0f) - ECUT;
                        pres_vir += fr * rij2;
                    }
                }
                
            }
            off = 0;
        }
    }

#endif

    // acumular epot y pres
    _epot += horizontal_sum(epot_v);
    pres_vir += horizontal_sum(pres_v);

#endif

    pres_vir /= (V * 3.0f);
    *pres = *temp * rho + pres_vir;
    *epot = _epot;

}

#else

__constant__ float device_rx[N];
__constant__ float device_ry[N];
__constant__ float device_rz[N];

__constant__ float device_L;
__constant__ float device_rcut2;

__global__ void forces_kernel(float* fx, float* fy, float* fz,
                              float* epot, float* pres_vir)
{
    int gtid = blockIdx.x * blockDim.x + threadIdx.x;	// global id

    if (gtid == 0) {
        *epot = 0.0f;
        *pres_vir = 0.0f;
    }

    for (int i = 0; i < gtid; i++) {

        float xi = device_rx[i];
        float yi = device_ry[i];
        float zi = device_rz[i];

        float xj = device_rx[gtid];
        float yj = device_ry[gtid];
        float zj = device_rz[gtid];

        // distancia mínima entre r_i y r_j
        float xij = xi - xj;
        xij = minimum_image(xij, device_L);
        float yij = yi - yj;
        yij = minimum_image(yij, device_L);
        float zij = zi - zj;
        zij = minimum_image(zij, device_L);

        float rij2 = xij * xij + yij * yij + zij * zij;

        if (rij2 <= device_rcut2) {
            float r2inv = 1.0f / rij2;
            float r6inv = r2inv * r2inv * r2inv;

            float fr = 24.0f * r2inv * r6inv * (2.0f * r6inv - 1.0f);

            atomicAdd(fx + i, fr * xij);
            atomicAdd(fy + i, fr * yij);
            atomicAdd(fz + i, fr * zij);

            fx[gtid] -= fr * xij;
            fy[gtid] -= fr * yij;
            fz[gtid] -= fr * zij;

            atomicAdd(epot, 4.0f * r6inv * (r6inv - 1.0f) - ECUT);
            atomicAdd(pres_vir, fr * rij2);
        }
    }
}

void forces(float* rxyz, float* fxyz, float* epot, float* pres,
            const float* temp, const float rho, const float V, const float L)
{
    for (int i = 0; i < 3 * N; i++) {
        fxyz[i] = 0.0f;
    }

    
    float rcut2 = RCUT * RCUT;
    float _epot, pres_vir;
    
    // GPU
    
    float *device_fxyz, *device_epot, *device_pres_vir, *device_L, *device_rcut2;
    cudaMalloc((void**)&device_fxyz, 3 * N * sizeof(float));
    cudaMalloc((void**)&device_epot, sizeof(float));
    cudaMalloc((void**)&device_pres_vir, sizeof(float));
    cudaMalloc((void**)&device_L, sizeof(float));
    cudaMalloc((void**)&device_rcut2, sizeof(float));
    
    float *rx = rxyz, *ry = rxyz + N, *rz = rxyz + 2 * N;
    float *device_fx = device_fxyz, *device_fy = device_fxyz + N, *device_fz = device_fxyz + 2 * N;
    
    cudaMemcpyToSymbol(device_rx, rx, 3 * N * sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(device_ry, ry, 3 * N * sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(device_rz, rz, 3 * N * sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(device_L, &L, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(device_rcut2, &rcut2, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpy(device_fxyz, fxyz, 3 * N * sizeof(float), cudaMemcpyHostToDevice);
    
    

    int nthreads = 256;
    int nblocks = (N + nthreads - 1) / nthreads; 

    forces_kernel<<<nblocks, nthreads>>>(device_fx, device_fy, device_fz,
                                         device_epot, device_pres_vir);
    
    
    cudaDeviceSynchronize();
    cudaMemcpy(fxyz, device_fxyz, 3 * N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&_epot, device_epot, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&pres_vir, device_pres_vir, sizeof(float), cudaMemcpyDeviceToHost);

    pres_vir /= (V * 3.0f);
    *pres = *temp * rho + pres_vir;
    *epot = _epot;

}

#endif



static float pbc(float cordi, const float cell_length)
{
    // condiciones periodicas de contorno coordenadas entre [0,L)
    /*
        if (cordi <= 0) {
            cordi += cell_length;
        } else if (cordi > cell_length) {
            cordi -= cell_length;
        }*/
    return cordi + cell_length * ((cordi <= 0) - (cordi > cell_length));
}


void velocity_verlet(float* rxyz, float* vxyz, float* fxyz, float* epot,
                        float* ekin, float* pres, float* temp, const float rho,
                        const float V, const float L)
{
    float *rx = rxyz, *ry = rxyz + N, *rz = rxyz + 2 * N;
    float *vx = vxyz, *vy = vxyz + N, *vz = vxyz + 2 * N;
    float *fx = fxyz, *fy = fxyz + N, *fz = fxyz + 2 * N;

    for (int i = 0; i < N; i += 1) { // actualizo posiciones
        rx[i] += vx[i] * DT + 0.5f * fx[i] * DT * DT;
        ry[i] += vy[i] * DT + 0.5f * fy[i] * DT * DT;
        rz[i] += vz[i] * DT + 0.5f * fz[i] * DT * DT;

        rx[i] = pbc(rx[i], L);
        ry[i] = pbc(ry[i], L);
        rz[i] = pbc(rz[i], L);

        vx[i] += 0.5f * fx[i] * DT;
        vy[i] += 0.5f * fy[i] * DT;
        vz[i] += 0.5f * fz[i] * DT;
    }

    forces(rxyz, fxyz, epot, pres, temp, rho, V, L); // actualizo fuerzas

    float sumv2 = 0.0f;
    for (int i = 0; i < N; i += 1) { // actualizo velocidades
        vx[i] += 0.5f * fx[i] * DT;
        vy[i] += 0.5f * fy[i] * DT;
        vz[i] += 0.5f * fz[i] * DT;

        sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    }

    *ekin = 0.5f * sumv2;
    *temp = sumv2 / (3.0f * N);
}
