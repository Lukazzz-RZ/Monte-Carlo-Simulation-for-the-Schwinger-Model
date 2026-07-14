#include "Lattice.h"

#include <errno.h>

#ifndef _WIN32
#include <sys/stat.h>
#include <sys/types.h>
#endif

///////     PARAMETROS GLOBALES     ///////////////////

double beta = 10.0;
double m = 0.0;
double theta = 0.0;
double sign = 1.0;

///////     PARAMETROS DE SIMULACION     ///////////////////

static const long N_THERMALIZATION = 100000L;
static const long N_MEASUREMENTS = 50176L;
static const int STEPS_BETWEEN_MEASUREMENTS = 10;
static const int JACKKNIFE_BLOCK_SIZE = 256;
static const long PROGRESS_EVERY = 4096L;

static const double THETA_OVER_PI_VALUES[] = {
    0.00,
    0.25,
    0.50,
    0.75,
    1.00,
    1.25,
    1.50,
    1.75,
    2.00
};

static const double MASS_OVER_G_VALUES[] = {
    1.00
};

static const char *OBSERVABLE_NAMES[N_OBS] = {
    "rho_d",
    "rho_d_sqr",
    "q",
    "q_sqr",
    "U_p",
    "U_p_sqr",
    "rho_m",
    "rho_m_sqr",
    "chiral",
    "chiral_sqr"
};

///////     FUNCIONES AUXILIARES     ///////////////////

/**
 * Se obtiene el numero de elementos de un vector estatico.
 */
#define ARRAY_SIZE(array) (sizeof(array) / sizeof((array)[0]))

/**
 * Se crea el directorio de resultados si no existe.
 *
 * @return true si el directorio esta disponible y false en caso contrario.
 */
static bool create_results_directory(void)
{
#ifdef _WIN32
    const int status = _mkdir("Results");
#else
    const int status = mkdir("Results", 0777);
#endif

    if (status == 0 || errno == EEXIST) {
        return true;
    }

    perror("Results");
    return false;
}

/**
 * Se sustituyen caracteres no deseados en una etiqueta de archivo.
 *
 * @param text Cadena que se modifica.
 */
static void sanitize_label(char *text)
{
    for (size_t i = 0; text[i] != '\0'; ++i) {
        if (text[i] == '.') {
            text[i] = 'p';
        } else if (text[i] == '-') {
            text[i] = 'm';
        }
    }
}

/**
 * Se construye una etiqueta unica a partir de los parametros de simulacion.
 *
 * @param spec Cadena de salida.
 * @param spec_size Tamano disponible en la cadena.
 * @param mass_over_g Masa fermionica en unidades del acoplo.
 * @param theta_over_pi Angulo topologico en unidades de pi.
 */
static void build_spec(
    char *spec,
    size_t spec_size,
    double mass_over_g,
    double theta_over_pi
)
{
    snprintf(
        spec,
        spec_size,
        "%dx%d_beta_%.3f_mg_%.3f_theta_%.3fpi",
        Nt,
        Nx,
        beta,
        mass_over_g,
        theta_over_pi
    );

    sanitize_label(spec);
}

/**
 * Se comprueba que los parametros permiten calcular los pesos duales.
 *
 * @return true si los parametros son validos y false en caso contrario.
 */
static bool validate_parameters(void)
{
    const double lambda = beta / 2.0 - theta / (4.0 * PI);
    const double lambda_conj = beta / 2.0 + theta / (4.0 * PI);

    if (beta <= 0.0) {
        fprintf(stderr, "Error: beta debe ser positiva.\n");
        return false;
    }

    if (m < 0.0) {
        fprintf(stderr, "Error: la masa no puede ser negativa.\n");
        return false;
    }

    if (lambda <= 0.0 || lambda_conj <= 0.0) {
        fprintf(
            stderr,
            "Error: los factores lambda deben ser positivos.\n"
        );
        return false;
    }

    if (Nx % 2 != 0) {
        fprintf(
            stderr,
            "Error: dimmer_start requiere un numero espacial par.\n"
        );
        return false;
    }

    return true;
}

/**
 * Se abre el archivo de produccion asociado a una serie.
 *
 * @param spec Identificador de la serie.
 * @param filename Nombre completo del archivo generado.
 * @param filename_size Tamano disponible para el nombre.
 * @return Puntero al archivo o NULL si no puede abrirse.
 */
static FILE *open_series_file(
    const char *spec,
    char *filename,
    size_t filename_size
)
{
    snprintf(
        filename,
        filename_size,
        "Results/Series_%s.dat",
        spec
    );

    FILE *file = fopen(filename, "w");

    if (file == NULL) {
        fprintf(
            stderr,
            "Error: no se pudo abrir el archivo %s\n",
            filename
        );
    }

    return file;
}

/**
 * Se escribe la cabecera del archivo de produccion.
 *
 * @param file Archivo de salida.
 * @param mass_over_g Masa fermionica en unidades del acoplo.
 * @param theta_over_pi Angulo topologico en unidades de pi.
 */
static void write_series_header(
    FILE *file,
    double mass_over_g,
    double theta_over_pi
)
{
    fprintf(file, "# Simulacion Monte Carlo del modelo de Schwinger\n");
    fprintf(file, "# Nt = %d\n", Nt);
    fprintf(file, "# Nx = %d\n", Nx);
    fprintf(file, "# CONTORNO_X = %d\n", CONTORNO_X);
    fprintf(file, "# beta = %.15e\n", beta);
    fprintf(file, "# m = %.15e\n", m);
    fprintf(file, "# m_over_g = %.15e\n", mass_over_g);
    fprintf(file, "# theta = %.15e\n", theta);
    fprintf(file, "# theta_over_pi = %.15e\n", theta_over_pi);
    fprintf(file, "# thermalization_steps = %ld\n", N_THERMALIZATION);
    fprintf(file, "# measurements = %ld\n", N_MEASUREMENTS);
    fprintf(
        file,
        "# steps_between_measurements = %d\n",
        STEPS_BETWEEN_MEASUREMENTS
    );
    fprintf(
        file,
        "# jackknife_block_size = %d\n",
        JACKKNIFE_BLOCK_SIZE
    );

    fprintf(file, "# columns: mc_step");

    for (int observable = 0; observable < N_OBS; ++observable) {
        fprintf(file, " %s", OBSERVABLE_NAMES[observable]);
    }

    fprintf(file, " sign");

    for (int observable = 0; observable < N_OBS; ++observable) {
        fprintf(
            file,
            " sign_%s",
            OBSERVABLE_NAMES[observable]
        );
    }

    fprintf(file, "\n");
}

/**
 * Se escribe una medida y los productos necesarios para la reponderacion.
 *
 * @param file Archivo de salida.
 * @param mc_step Paso Monte Carlo de produccion.
 * @param obs Vector de observables.
 */
static void write_measurement(
    FILE *file,
    long mc_step,
    const double *obs
)
{
    fprintf(file, "%ld", mc_step);

    for (int observable = 0; observable < N_OBS; ++observable) {
        fprintf(file, " %.15e", obs[observable]);
    }

    fprintf(file, " %.15e", sign);

    for (int observable = 0; observable < N_OBS; ++observable) {
        fprintf(file, " %.15e", sign * obs[observable]);
    }

    fprintf(file, "\n");
}

/**
 * Se muestra el estado de una serie durante la fase de produccion.
 *
 * @param measurement Numero de medida realizada.
 * @param obs Vector de observables.
 * @param mean_sign Signo medio acumulado.
 */
static void print_progress(
    long measurement,
    const double *obs,
    double mean_sign
)
{
    printf(
        "Medida %ld / %ld | U_p = %.9f | rho_d = %.9f | "
        "q = %.9e | signo medio = %.6f\n",
        measurement,
        N_MEASUREMENTS,
        obs[OBS_U_P],
        obs[OBS_RHO_D],
        obs[OBS_Q],
        mean_sign
    );
}

///////     SIMULACION DE UN PUNTO     ///////////////////

/**
 * Se termaliza y se mide un punto del plano de parametros.
 *
 * @param mass_over_g Masa fermionica en unidades del acoplo.
 * @param theta_over_pi Angulo topologico en unidades de pi.
 * @return true si la simulacion termina correctamente y false si falla.
 */
static bool run_parameter_point(
    double mass_over_g,
    double theta_over_pi
)
{
    char spec[256];
    char filename[320];

    theta = theta_over_pi * PI;
    m = mass_over_g / sqrt(beta);
    sign = 1.0;

    if (!validate_parameters()) {
        return false;
    }

    build_spec(
        spec,
        sizeof(spec),
        mass_over_g,
        theta_over_pi
    );

    dimmer_start();

    if (!check_lattice()) {
        fprintf(
            stderr,
            "Error: configuracion inicial no valida para %s\n",
            spec
        );
        return false;
    }

    printf("\n========================================\n");
    printf("Serie: %s\n", spec);
    printf(
        "beta = %.6f | m = %.9f | m/g = %.6f | theta/pi = %.3f\n",
        beta,
        m,
        mass_over_g,
        theta_over_pi
    );
    printf(
        "Termalizacion: %ld pasos Monte Carlo\n",
        N_THERMALIZATION
    );

    thermalizate((int)N_THERMALIZATION, spec);

    if (!check_lattice()) {
        fprintf(
            stderr,
            "Error: configuracion no valida tras termalizar %s\n",
            spec
        );
        return false;
    }

    FILE *series = open_series_file(
        spec,
        filename,
        sizeof(filename)
    );

    if (series == NULL) {
        return false;
    }

    write_series_header(
        series,
        mass_over_g,
        theta_over_pi
    );

    double obs[N_OBS];
    double sign_sum = 0.0;
    long mc_step = 0L;

    for (
        long measurement = 1L;
        measurement <= N_MEASUREMENTS;
        ++measurement
    ) {
        for (
            int step = 0;
            step < STEPS_BETWEEN_MEASUREMENTS;
            ++step
        ) {
            monte_carlo_step();
            ++mc_step;
        }

        measure_all(obs);
        sign_sum += sign;

        write_measurement(series, mc_step, obs);

        if (
            measurement % PROGRESS_EVERY == 0 ||
            measurement == N_MEASUREMENTS
        ) {
            const double mean_sign =
                sign_sum / (double)measurement;

            print_progress(
                measurement,
                obs,
                mean_sign
            );

            fflush(series);
        }
    }

    fclose(series);

    if (!check_lattice()) {
        fprintf(
            stderr,
            "Error: configuracion final no valida para %s\n",
            spec
        );
        return false;
    }

    printf("Archivo generado: %s\n", filename);
    printf(
        "Signo medio bruto: %.12e\n",
        sign_sum / (double)N_MEASUREMENTS
    );
    printf("Serie terminada correctamente.\n");

    return true;
}

///////     PROGRAMA PRINCIPAL     ///////////////////

/**
 * Se ejecuta el barrido definido sobre masa y angulo topologico.
 *
 * @return EXIT_SUCCESS si todas las series terminan correctamente.
 */
int main(void)
{
    srand((unsigned int)time(NULL));

    if (!create_results_directory()) {
        return EXIT_FAILURE;
    }

    if (N_MEASUREMENTS % JACKKNIFE_BLOCK_SIZE != 0) {
        fprintf(
            stderr,
            "Error: el numero de medidas debe ser multiplo de %d.\n",
            JACKKNIFE_BLOCK_SIZE
        );
        return EXIT_FAILURE;
    }

    init_moves(dx_p, dx_n, dt_p, dt_n);

    WORM_DEBUG = false;
    WORMCORE_DEBUG = false;

    for (
        size_t theta_index = 0;
        theta_index < ARRAY_SIZE(THETA_OVER_PI_VALUES);
        ++theta_index
    ) {
        for (
            size_t mass_index = 0;
            mass_index < ARRAY_SIZE(MASS_OVER_G_VALUES);
            ++mass_index
        ) {
            if (!run_parameter_point(
                    MASS_OVER_G_VALUES[mass_index],
                    THETA_OVER_PI_VALUES[theta_index])) {
                return EXIT_FAILURE;
            }
        }
    }

    printf("\nTodas las series han terminado correctamente.\n");

    return EXIT_SUCCESS;
}
