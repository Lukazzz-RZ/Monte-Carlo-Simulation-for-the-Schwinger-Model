#include "Lattice.h"

///////     VISUALIZACION DEL RETICULO     ///////////////////

/**
 * Se muestra una region del reticulo en consola, incluyendo monomeros,
 * dimeros, loops y numeros de ocupacion de plaqueta.
 */
void show_lattice(void)
{
    const int initial_t = 0;
    const int initial_x = 0;
    const int final_t = 3 * Nt / 4;
    const int final_x = Nx;

    system("cls");

    for (int t = initial_t; t < final_t; ++t) {
        for (int x = initial_x; x < final_x; ++x) {
            if (monomer[t][x]) {
                printf("\033[0;32mx\033[0m");
            } else {
                printf("o");
            }

            if (link_h[t][x] == 1) {
                printf("\033[0;33m>>\033[0m");
            } else if (link_h[t][x] == -1) {
                printf("\033[0;33m<<\033[0m");
            } else if (link_h[t][x] == 2) {
                printf("==");
            } else if (x < final_x - 1) {
                printf("  ");
            }
        }

        printf("\n");

        for (int x = initial_x; x < final_x; ++x) {
            if (link_v[t][x] == -1) {
                printf("\033[0;33m^\033[0m");
            } else if (link_v[t][x] == 1) {
                printf("\033[0;33mv\033[0m");
            } else if (link_v[t][x] == 2) {
                printf(":");
            } else {
                printf(" ");
            }

            if (plaquette[t][x] != 0) {
                printf("\033[96m%2d\033[0m", plaquette[t][x]);
            } else {
                printf("%2d", plaquette[t][x]);
            }
        }

        printf("\n");
    }

    printf("\n");
}

/**
 * Se muestran los pesos de Bessel calculados y almacenados.
 */
void show_probs(void)
{
    printf("Pesos de Bessel:\n");

    for (int order = 0; order < NORDERS; ++order) {
        if (bessel_used[order]) {
            printf(
                "I_%d(beta) = %f\n",
                order,
                bessel_used_orders[order]
            );
        }
    }

    printf("\n");
}

///////     ESCRITURA DE RESULTADOS     ///////////////////

/**
 * Se escriben los observables de una configuracion en un archivo de texto.
 *
 * @param overwrite Indica si se sobrescribe el archivo existente.
 * @param obs Vector de observables.
 * @param n_obs Numero de observables almacenados.
 * @param spec Identificador empleado en el nombre del archivo.
 */
void write_results(
    bool overwrite,
    double *obs,
    int n_obs,
    const char *spec
)
{
    char filename[256];

    _mkdir("Results");

    snprintf(
        filename,
        sizeof(filename),
        "Results/Observables_%s.txt",
        spec
    );

    FILE *file = fopen(filename, overwrite ? "w" : "a");

    if (file == NULL) {
        printf(
            "Error: no se pudo abrir el archivo %s\n",
            filename
        );
        return;
    }

    if (overwrite) {
        fprintf(file, "# ========================================\n");
        fprintf(file, "# Archivo de medidas Monte Carlo\n");
        fprintf(file, "# ========================================\n");
        fprintf(file, "# beta       = %.12e\n", beta);
        fprintf(file, "# m          = %.12e\n", m);
        fprintf(file, "# theta      = %.12e\n", theta);
        fprintf(file, "# theta/pi   = %.12e\n", theta / PI);
        fprintf(file, "# Nt         = %d\n", Nt);
        fprintf(file, "# Nx         = %d\n", Nx);
        fprintf(file, "# CONTORNO_X = %d\n", CONTORNO_X);
        fprintf(file, "# V          = %d\n", Nt * (Nx - CONTORNO_X));
        fprintf(file, "# n_obs      = %d\n", n_obs);
        fprintf(file, "#\n");
        fprintf(file, "# Columnas:\n");
        fprintf(file, "# 0  rho_d\n");
        fprintf(file, "# 1  rho_d^2\n");
        fprintf(file, "# 2  q\n");
        fprintf(file, "# 3  q^2\n");
        fprintf(file, "# 4  U_p\n");
        fprintf(file, "# 5  U_p^2\n");
        fprintf(file, "# 6  rho_m\n");
        fprintf(file, "# 7  rho_m^2\n");
        fprintf(file, "# 8  chiral\n");
        fprintf(file, "# 9  chiral^2\n");
        fprintf(file, "# 10 sign\n");
        fprintf(file, "# ========================================\n");
    }

    for (int i = 0; i < n_obs; ++i) {
        fprintf(file, "%.12e", obs[i]);

        if (i < n_obs - 1) {
            fprintf(file, " ");
        }
    }

    fprintf(file, "\n");
    fclose(file);
}

/**
 * Se guarda la plaqueta media y la densidad de dimeros durante la
 * termalizacion.
 *
 * @param step Indice del paso Monte Carlo.
 * @param U_p Valor medio de la plaqueta.
 * @param rho_d Densidad media de dimeros.
 * @param spec Identificador empleado en el nombre del archivo.
 */
void write_thermal(
    int step,
    double U_p,
    double rho_d,
    const char *spec
)
{
    char filename[256];

    _mkdir("Results");

    snprintf(
        filename,
        sizeof(filename),
        "Results/Thermalization_%s.txt",
        spec
    );

    FILE *file = fopen(filename, step == 0 ? "w" : "a");

    if (file == NULL) {
        return;
    }

    fprintf(file, "%d %f %f\n", step, U_p, rho_d);
    fclose(file);
}

///////     REPRESENTACION DE RESULTADOS     ///////////////////

/**
 * Se representa la evolucion de la plaqueta media y la densidad de dimeros
 * durante la termalizacion.
 *
 * @param spec Identificador empleado en el nombre del archivo.
 */
void plot_thermal(const char *spec)
{
    char filename[256];

    snprintf(
        filename,
        sizeof(filename),
        "Results/Thermalization_%s.txt",
        spec
    );

    FILE *gnuplot = popen("gnuplot -persistent", "w");

    if (gnuplot == NULL) {
        perror("gnuplot");
        return;
    }

    fprintf(gnuplot, "set title 'Termalizado'\n");
    fprintf(gnuplot, "set xlabel 'Paso Monte Carlo'\n");
    fprintf(gnuplot, "set ylabel 'Valor'\n");
    fprintf(gnuplot, "set grid\n");
    fprintf(gnuplot, "set key left top\n");

    fprintf(
        gnuplot,
        "plot '%s' using 1:2 with lines lw 2 title 'U_p', "
        "'%s' using 1:3 with lines lw 2 title 'densidad de dimeros'\n",
        filename,
        filename
    );

    fflush(gnuplot);
    pclose(gnuplot);
}

/**
 * Se representan los datos almacenados en el archivo de observables mediante
 * gnuplot.
 *
 * @param spec Identificador empleado en el nombre del archivo.
 */
void plot_obs(const char *spec)
{
    char filename[256];

    snprintf(
        filename,
        sizeof(filename),
        "Results/Observables_%s.txt",
        spec
    );

    FILE *gnuplot = popen("gnuplot -persistent", "w");

    if (gnuplot == NULL) {
        perror("gnuplot");
        return;
    }

    fprintf(gnuplot, "set title ''\n");
    fprintf(gnuplot, "set xlabel 'beta'\n");
    fprintf(gnuplot, "set ylabel 'U_p'\n");
    fprintf(gnuplot, "set grid\n");
    fprintf(gnuplot, "set key left top\n");

    fprintf(
        gnuplot,
        "plot '%s' using 1:2 with points lw 2 title 'U_p'\n",
        filename
    );

    fflush(gnuplot);
    pclose(gnuplot);
}
