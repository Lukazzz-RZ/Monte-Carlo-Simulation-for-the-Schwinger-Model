#include "Lattice.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

///////     DIMENSIONES FISICAS     ///////////////////

#define X_PHYS (Nx - CONTORNO_X)
#define VOLUME (Nt * X_PHYS)

///////     FUNCIONES AUXILIARES     ///////////////////

/**
 * Se calcula la funcion de Bessel modificada usando un orden no negativo.
 *
 * @param n Orden entero de la funcion.
 * @param x Argumento real de la funcion.
 * @return Valor aproximado de I_n(x).
 */
double besselmod_safe(int n, double x)
{
    return besselmod(abs(n), x);
}

///////     DENSIDAD DE DIMEROS     ///////////////////

/**
 * Se calcula la densidad media de dimeros horizontales y verticales.
 *
 * @return Densidad media de dimeros.
 */
double measure_rho_d(void)
{
    double count = 0.0;

    if (X_PHYS <= 0) {
        fprintf(
            stderr,
            "Error: volumen espacial invalido en measure_rho_d\n"
        );
        return 0.0;
    }

    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < X_PHYS; ++x) {
            if (link_h[t][x] == 2) {
                count += 1.0;
            }

            if (link_v[t][x] == 2) {
                count += 1.0;
            }
        }
    }

    return count / (2.0 * Nt * X_PHYS);
}

/**
 * Se calcula el cuadrado de la densidad media de dimeros.
 *
 * @return Cuadrado de la densidad media de dimeros.
 */
double measure_rho_d_sqr(void)
{
    const double rho_d = measure_rho_d();

    return rho_d * rho_d;
}

///////     DENSIDAD TOPOLOGICA     ///////////////////

/**
 * Se calcula la contribucion local exacta a la densidad topologica.
 *
 * @param p Numero de ocupacion de plaqueta.
 * @return Contribucion local a la densidad topologica.
 */
double local_q_exact(int p)
{
    const int order = abs(p);
    const double argument = 2.0 * eta * eta_conj;

    const double bessel = besselmod_safe(order, argument);
    const double bessel_minus = besselmod_safe(order - 1, argument);
    const double bessel_plus = besselmod_safe(order + 1, argument);

    if (bessel == 0.0) {
        fprintf(stderr, "Error: Bessel nulo en local_q_exact\n");
        return 0.0;
    }

    const double term_bessel =
        (1.0 / (8.0 * PI)) *
        (eta_conj / eta - eta / eta_conj) *
        (bessel_minus + bessel_plus) / bessel;

    const double term_occupation =
        ((double)p / (8.0 * PI)) *
        (
            1.0 / (eta * eta) +
            1.0 / (eta_conj * eta_conj)
        );

    return term_bessel + term_occupation;
}

/**
 * Se calcula la densidad topologica media del reticulo.
 *
 * @return Densidad topologica media.
 */
double measure_q(void)
{
    double q_sum = 0.0;

    if (VOLUME <= 0) {
        fprintf(stderr, "Error: volumen invalido en measure_q\n");
        return 0.0;
    }

    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < X_PHYS; ++x) {
            q_sum += local_q_exact(plaquette[t][x]);
        }
    }

    return q_sum / (double)VOLUME;
}

/**
 * Se calcula el cuadrado de la densidad topologica media.
 *
 * @return Cuadrado de la densidad topologica media.
 */
double measure_q_sqr(void)
{
    const double q = measure_q();

    return q * q;
}

///////     PLAQUETA MEDIA     ///////////////////

/**
 * Se calcula la contribucion local exacta a la plaqueta media.
 *
 * @param p Numero de ocupacion de plaqueta.
 * @return Contribucion local a la plaqueta media.
 */
double local_U_p_exact(int p)
{
    const int order = abs(p);
    const double argument = 2.0 * eta * eta_conj;

    const double bessel = besselmod_safe(order, argument);
    const double bessel_minus = besselmod_safe(order - 1, argument);
    const double bessel_plus = besselmod_safe(order + 1, argument);

    if (bessel == 0.0) {
        fprintf(stderr, "Error: Bessel nulo en local_U_p_exact\n");
        return 0.0;
    }

    const double term_bessel =
        (beta / (4.0 * eta * eta_conj)) *
        (bessel_minus + bessel_plus) / bessel;

    const double term_occupation =
        0.25 * (double)p *
        (
            1.0 / (eta * eta) -
            1.0 / (eta_conj * eta_conj)
        );

    return term_bessel + term_occupation;
}

/**
 * Se calcula el valor medio de la plaqueta en el reticulo.
 *
 * @return Valor medio de la plaqueta.
 */
double measure_U_p(void)
{
    double plaquette_sum = 0.0;

    if (VOLUME <= 0) {
        fprintf(stderr, "Error: volumen invalido en measure_U_p\n");
        return 0.0;
    }

    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < X_PHYS; ++x) {
            plaquette_sum += local_U_p_exact(plaquette[t][x]);
        }
    }

    return plaquette_sum / (double)VOLUME;
}

/**
 * Se calcula el cuadrado del valor medio de la plaqueta.
 *
 * @return Cuadrado del valor medio de la plaqueta.
 */
double measure_U_p_sqr(void)
{
    const double plaquette_mean = measure_U_p();

    return plaquette_mean * plaquette_mean;
}

///////     DENSIDAD DE MONOMEROS     ///////////////////

/**
 * Se calcula la densidad media de monomeros.
 *
 * @return Densidad media de monomeros.
 */
double measure_monomer_density(void)
{
    double count = 0.0;

    if (VOLUME <= 0) {
        fprintf(
            stderr,
            "Error: volumen invalido en measure_monomer_density\n"
        );
        return 0.0;
    }

    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < X_PHYS; ++x) {
            if (monomer[t][x] == 1) {
                count += 1.0;
            }
        }
    }

    return count / (double)VOLUME;
}

/**
 * Se calcula el cuadrado de la densidad media de monomeros.
 *
 * @return Cuadrado de la densidad media de monomeros.
 */
double measure_monomer_density_sqr(void)
{
    const double rho_m = measure_monomer_density();

    return rho_m * rho_m;
}

///////     MEDIDA DE OBSERVABLES     ///////////////////

/**
 * Se calculan los observables principales y sus cuadrados para la
 * configuracion actual.
 *
 * @param obs Vector de salida indexado mediante las constantes OBS.
 */
void measure_all(double *obs)
{
    const double rho_d = measure_rho_d();
    const double q = measure_q();
    const double plaquette_mean = measure_U_p();
    const double rho_m = measure_monomer_density();

    double chiral = 0.0;

    if (m == 0.0) {
        fprintf(
            stderr,
            "Error: masa nula al calcular condensado quiral\n"
        );
    } else {
        chiral = rho_m / m;
    }

    obs[OBS_RHO_D] = rho_d;
    obs[OBS_RHO_D_SQR] = rho_d * rho_d;

    obs[OBS_Q] = q;
    obs[OBS_Q_SQR] = q * q;

    obs[OBS_U_P] = plaquette_mean;
    obs[OBS_U_P_SQR] = plaquette_mean * plaquette_mean;

    obs[OBS_RHO_M] = rho_m;
    obs[OBS_RHO_M_SQR] = rho_m * rho_m;

    obs[OBS_CHIRAL] = chiral;
    obs[OBS_CHIRAL_SQR] = chiral * chiral;
}
