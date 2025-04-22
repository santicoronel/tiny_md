void forces(float* rxyz, float* fxyz, float* epot, float* pres,
            const float* temp, const float rho, const float V, const float L)
{
    // Calculate Lennard-Jones forces (12-6)

    for (int i = 0; i < 3 * N; i++) {
        fxyz[i] = 0.0f;
    }

    float *rx = rxyz, *ry = rxyz + N, *rz = rxyz + 2 * N;
    float *fx = fxyz, *fy = fxyz + N, *fz = fxyz + 2 * N;

    float rcut2 = RCUT * RCUT;
    *epot = 0.0f;
    float pres_vir = 0.0f;

    for (int i = 0; i < (N - 1); i++) {

        float xi = rx[i];
        float yi = ry[i];
        float zi = rz[i];

        for (int j = i + 1; j < N; j++) {

            float xj = rx[j];
            float yj = ry[j];
            float zj = rz[j];

            // Minimum image distance between r_i and r_j
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

    pres_vir /= (V * 3.0f);
    *pres = *temp * rho + pres_vir;
}