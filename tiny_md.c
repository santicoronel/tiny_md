#define _XOPEN_SOURCE 500  // M_PI
#include "core.h"
#include "parameters.h"
#include "wtime.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int main()
{
    //FILE *file_xyz, *file_thermo;
    FILE *file_out;
    //file_xyz = fopen("trajectory.xyz", "w");
    //file_thermo = fopen("thermo.log", "w");
    file_out = fopen("log.csv", "a");
    float Ekin, Epot, Temp, Pres; // variables macroscopicas
    float Rho, cell_V, cell_L, tail, Etail, Ptail;
    float *rxyz, *vxyz, *fxyz; // variables microscopicas

    rxyz = (float*)malloc(3 * N * sizeof(float));
    vxyz = (float*)malloc(3 * N * sizeof(float));
    fxyz = (float*)malloc(3 * N * sizeof(float));

    printf("# Número de partículas:      %d\n", N);
    printf("# Temperatura de referencia: %.2f\n", T0);
    printf("# Pasos de equilibración:    %d\n", TEQ);
    printf("# Pasos de medición:         %d\n", TRUN - TEQ);
    printf("# (mediciones cada %d pasos)\n", TMES);
    printf("# densidad, volumen, energía potencial media, presión media\n");
    //fprintf(file_thermo, "# t Temp Pres Epot Etot\n");

    srand(SEED);
    float t = 0.0f, sf;
    float Rhob;
    Rho = RHOI;
    init_pos(rxyz, Rho);
    float start = wtime();
    for (int m = 0; m < 9; m++) {
        Rhob = Rho;
        Rho = RHOI - 0.1f * (float)m;
        cell_V = (float)N / Rho;
        cell_L = cbrt(cell_V);
        tail = 16.0f * M_PI * Rho * ((2.0f / 3.0f) * (float)pow(RCUT, -9) - (float)pow(RCUT, -3)) / 3.0f;
        Etail = tail * (float)N;
        Ptail = tail * Rho;

        int i = 0;
        sf = cbrt(Rhob / Rho);
        for (int k = 0; k < 3 * N; k++) { // reescaleo posiciones a nueva densidad
            rxyz[k] *= sf;
        }
        init_vel(vxyz, &Temp, &Ekin);
        forces(rxyz, fxyz, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);

        for (i = 1; i < TEQ; i++) { // loop de equilibracion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < 3 * N; k++) { // reescaleo de velocidades
                vxyz[k] *= sf;
            }
        }

        int mes = 0;
        float epotm = 0.0f, presm = 0.0f;
        for (i = TEQ; i < TRUN; i++) { // loop de medicion

            velocity_verlet(rxyz, vxyz, fxyz, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < 3 * N; k++) { // reescaleo de velocidades
                vxyz[k] *= sf;
            }

            if (i % TMES == 0) {
                Epot += Etail;
                Pres += Ptail;

                epotm += Epot;
                presm += Pres;
                mes++;

                //fprintf(file_thermo, "%f %f %f %f %f\n", t, Temp, Pres, Epot, Epot + Ekin);
                //fprintf(file_xyz, "%d\n\n", N);
                //for (int k = 0; k < 3 * N; k += 3) {
                //    fprintf(file_xyz, "Ar %e %e %e\n", rxyz[k + 0], rxyz[k + 1], rxyz[k + 2]);
                //}
            }

            t += DT;
        }
        printf("%f\t%f\t%f\t%f\n", Rho, cell_V, epotm / (float)mes, presm / (float)mes);
    }

    float elapsed = wtime() - start;
    printf("# Tiempo total de simulación = %f segundos\n", elapsed);
    printf("# Tiempo simulado = %f [fs]\n", t * 1.6);
    printf("# ns/day = %f\n", (1.6e-6 * t) / elapsed * 86400);
    //                       ^1.6 fs -> ns       ^sec -> day
    // Métrica
    float metrica = N * N / elapsed;
    printf ("# N^2/t = %f\n", metrica);

    fprintf(file_out, "%d,%f,%f\n", N, elapsed, metrica);

    // Cierre de archivos
    //fclose(file_thermo);
    //fclose(file_xyz);
    fclose(file_out);

    // Liberacion de memoria
    free(rxyz);
    free(fxyz);
    free(vxyz);
    return 0;
}
