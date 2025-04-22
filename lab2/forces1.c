void forces(float* rxyz, float* fxyz, float* epot, float* pres,
            const float* temp, const float rho, const float V, const float L) {
    // Calculate Lennard-Jones forces (12-6)

    for (int i = 0; i < 3 * N; i++) {
        fxyz[i] = 0.0f;
    }

    float *rx = rxyz, *ry = rxyz + N, *rz = rxyz + 2 * N;
    float *fx = fxyz, *fy = fxyz + N, *fz = fxyz + 2 * N;

    float rcut2 = RCUT * RCUT;
    *epot = 0.0f;
    float pres_vir = 0.0f;

    __m256 rcut2_v = _mm256_set1_ps(rcut2);
    __m256 ecut = _mm256_set1_ps(ECUT);
    __m256 epot_v = _mm256_setzero_ps();
    __m256 pres_v = _mm256_setzero_ps();

    for (unsigned int i = 0; i < N - 1; i++) {
        __m256 xi = _mm256_set1_ps(rx[i]);
        __m256 yi = _mm256_set1_ps(ry[i]);
        __m256 zi = _mm256_set1_ps(rz[i]);

        __m256 fxi = _mm256_setzero_ps();
        __m256 fyi = _mm256_setzero_ps();
        __m256 fzi = _mm256_setzero_ps();

        unsigned int j;

        for (j = i + 1; j + 8 < N; j += 8) {
            __m256 xj = _mm256_loadu_ps(rx + j);
            __m256 yj = _mm256_loadu_ps(ry + j);
            __m256 zj = _mm256_loadu_ps(rz + j);

            __m256 xij = _mm256_sub_ps(xi, xj);
            xij = minimum_image_v(xij, L);

            __m256 yij = _mm256_sub_ps(yi, yj);
            yij = minimum_image_v(yij, L);

            __m256 zij = _mm256_sub_ps(zi, zj);
            zij = minimum_image_v(zij, L);

            __m256 rij2 = _mm256_add_ps(_mm256_mul_ps(xij, xij),
                                        _mm256_add_ps(_mm256_mul_ps(yij, yij),
                                                      _mm256_mul_ps(zij, zij)));

            __m256 mask_rcut = _mm256_cmp_ps(rij2, rcut2_v, _CMP_LE_OQ);

            if (!_mm256_movemask_ps(mask_rcut)) continue;

            __m256 r2inv = _mm256_div_ps(_mm256_set1_ps(1.0f), rij2);
            __m256 r6inv = _mm256_mul_ps(r2inv, _mm256_mul_ps(r2inv, r2inv));
            __m256 fr = _mm256_mul_ps(_mm256_set1_ps(24.0f),
                                      _mm256_mul_ps(r2inv,
                                                    _mm256_mul_ps(r6inv,
                                                                  _mm256_sub_ps(_mm256_mul_ps(_mm256_set1_ps(2.0f), r6inv),
                                                                                _mm256_set1_ps(1.0f)))));

            __m256 dfx = _mm256_mul_ps(fr, xij);
            __m256 dfy = _mm256_mul_ps(fr, yij);
            __m256 dfz = _mm256_mul_ps(fr, zij);

            fxi = _mm256_add_ps(fxi, _mm256_and_ps(dfx, mask_rcut));
            fyi = _mm256_add_ps(fyi, _mm256_and_ps(dfy, mask_rcut));
            fzi = _mm256_add_ps(fzi, _mm256_and_ps(dfz, mask_rcut));

            _mm256_storeu_ps(fx + j, _mm256_sub_ps(_mm256_loadu_ps(fx + j),
                                                   _mm256_and_ps(dfx, mask_rcut)));
            _mm256_storeu_ps(fy + j, _mm256_sub_ps(_mm256_loadu_ps(fy + j),
                                                   _mm256_and_ps(dfy, mask_rcut)));
            _mm256_storeu_ps(fz + j, _mm256_sub_ps(_mm256_loadu_ps(fz + j),
                                                   _mm256_and_ps(dfz, mask_rcut)));

            __m256 ep = _mm256_sub_ps(_mm256_mul_ps(_mm256_set1_ps(4.0f),
                                                    _mm256_mul_ps(r6inv,
                                                                  _mm256_sub_ps(r6inv, _mm256_set1_ps(1.0f)))),
                                      ecut);
            epot_v = _mm256_add_ps(epot_v, _mm256_and_ps(ep, mask_rcut));

            __m256 pv = _mm256_mul_ps(fr, rij2);
            pres_v = _mm256_add_ps(pres_v, _mm256_and_ps(pv, mask_rcut));
        }

        fx[i] += horizontal_sum(fxi);
        fy[i] += horizontal_sum(fyi);
        fz[i] += horizontal_sum(fzi);

        for (; j < N; j++) {
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

    *epot += horizontal_sum(epot_v);
    pres_vir += horizontal_sum(pres_v);


    pres_vir /= (V * 3.0f);
    *pres = *temp * rho + pres_vir;
}