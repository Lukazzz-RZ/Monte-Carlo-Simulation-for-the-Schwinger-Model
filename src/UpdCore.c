#include "Lattice.h"

///////     CONTROL DE DEPURACION     ///////////////////

bool WORMCORE_DEBUG = false;

///////     PESOS Y METROPOLIS     ///////////////////

/**
 * Se calcula el peso asociado a un numero de ocupacion de plaqueta.
 *
 * @param p Numero de ocupacion de plaqueta.
 * @return Peso de la plaqueta.
 */
double comp_weight(int p)
{
    const int order = abs(p);
    const double bessel_argument = 2.0 * eta * eta_conj;

    if (!bessel_used[order]) {
        bessel_used_orders[order] = besselmod(p, bessel_argument);
        bessel_used[order] = true;
    }

    return pow(eta / eta_conj, p) * bessel_used_orders[order];
}

/**
 * Se aplica el criterio de Metropolis al cambio entre dos ocupaciones de
 * plaqueta.
 *
 * @param p_ini Ocupacion inicial.
 * @param p_fin Ocupacion propuesta.
 * @return true si se acepta la propuesta y false en caso contrario.
 */
bool strong_coup_metropolis(int p_ini, int p_fin)
{
    const double initial_weight = comp_weight(p_ini);
    const double final_weight = comp_weight(p_fin);
    const double ratio = final_weight / initial_weight;
    const double acceptance = (ratio < 1.0) ? ratio : 1.0;

    return rand_double() < acceptance;
}

/**
 * Se aplica el criterio de Metropolis a dos pesos positivos.
 *
 * @param p_ini Peso inicial.
 * @param p_fin Peso propuesto.
 * @return true si se acepta la propuesta y false en caso contrario.
 */
bool metropolis(double p_ini, double p_fin)
{
    if (p_ini < 0.0 || p_fin <= 0.0) {
        return false;
    }

    if (p_ini <= 0.0 || p_fin >= p_ini) {
        return true;
    }

    return rand_double() < (p_fin / p_ini);
}

///////     BARRIDO DE PLAQUETAS     ///////////////////

/**
 * Se realiza un barrido secuencial completo de actualizaciones de plaqueta.
 */
void plaquette_sweeep(void)
{
    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < Nx - CONTORNO_X; ++x) {
            selected_plaq_update(t, x);
        }
    }
}

///////     PASO MONTE CARLO     ///////////////////

/**
 * Se realiza un paso Monte Carlo formado por un barrido de plaquetas,
 * maximo entre uno y V dividido entre diez intentos de gusano, y un segundo
 * barrido de plaquetas.
 */
void monte_carlo_step(void)
{
    const int volume = Nt * Nx;
    const int worm_attempts =
        (volume / 10 > 0) ? volume / 10 : 1;

    plaquette_sweeep();

    for (int attempt = 0; attempt < worm_attempts; ++attempt) {
        run_worm();
    }

    plaquette_sweeep();
}

///////     TERMALIZACION     ///////////////////

/**
 * Se termaliza la configuracion mediante pasos Monte Carlo completos y se
 * registra la plaqueta media y la densidad de dimeros cada diez pasos.
 *
 * @param n_steps Numero de pasos Monte Carlo de termalizacion.
 * @param spec Identificador empleado para guardar los datos.
 */
void thermalizate(int n_steps, const char *spec)
{
    printf("Termalizando...\n");

    for (int step = 0; step < n_steps; ++step) {
        monte_carlo_step();

        if (step % 10 == 0) {
            const double plaquette_mean = measure_U_p();
            const double dimer_density = measure_rho_d();

            write_thermal(
                step,
                plaquette_mean,
                dimer_density,
                spec
            );
        }
    }

    printf("Termalizacion terminada.\n\n");
}
