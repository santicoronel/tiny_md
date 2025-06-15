#ifndef FORCES_H
#define FORCES_H

#ifdef __cplusplus
extern "C" {
#endif

void forces(float* rxyz, float* fxyz, float* epot, float* pres,
            const float* temp, const float rho, const float V, const float L);

#ifdef __cplusplus
}
#endif
#endif