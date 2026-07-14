#include "Lattice.h"

///////     VARIABLES DEL RETICULO     ///////////////////

int link_h[Nt][Nx];
int link_v[Nt][Nx];
int plaquette[Nt][Nx];
int monomer[Nt][Nx];
bool fixed[Nt][Nx];
bool block[Nt][Nx];

///////     VARIABLES DE PRECALCULO     ///////////////////

double bessel_used_orders[NORDERS];
bool bessel_used[NORDERS];
double eta;
double eta_conj;
int dx_p[Nx];
int dx_n[Nx];
int dt_p[Nt];
int dt_n[Nt];

///////     PARAMETROS NUMERICOS     ///////////////////

enum { BESSEL_MAX_TERMS = 200 };
static const double BESSEL_TOLERANCE = 1e-64;

///////     PRECALCULOS     ///////////////////

/**
 * Se inicializan los desplazamientos hacia los sitios vecinos del reticulo.
 *
 * @param dx_plus Desplazamientos espaciales en sentido positivo.
 * @param dx_minus Desplazamientos espaciales en sentido negativo.
 * @param dt_plus Desplazamientos temporales en sentido positivo.
 * @param dt_minus Desplazamientos temporales en sentido negativo.
 */
void init_moves(int *dx_plus, int *dx_minus, int *dt_plus, int *dt_minus)
{
    for (int i = 0; i < Nx - 1; ++i) {
        dx_plus[i] = 1;
        dx_minus[i + 1] = -1;
    }

    dx_plus[Nx - 1] = -(Nx - 1) * (1 - CONTORNO_X);
    dx_minus[0] = (Nx - 1) * (1 - CONTORNO_X);

    for (int i = 0; i < Nt - 1; ++i) {
        dt_plus[i] = 1;
        dt_minus[i + 1] = -1;
    }

    dt_plus[Nt - 1] = -(Nt - 1);
    dt_minus[0] = Nt - 1;
}

/**
 * Se calcula la funcion de Bessel modificada de primera especie mediante su
 * expansion en serie de potencias.
 *
 * @param n Orden entero de la funcion.
 * @param x Argumento real de la funcion.
 * @return Aproximacion de I_n(x).
 */
double besselmod(int n, double x)
{
    const int order = (n < 0) ? -n : n;
    const double half_x = x / 2.0;
    double term = pow(half_x, order) / tgamma(order + 1);
    double sum = term;

    for (int k = 1; k < BESSEL_MAX_TERMS; ++k) {
        term *= (half_x * half_x) / (k * (order + k));
        sum += term;

        if (fabs(term) < BESSEL_TOLERANCE) {
            break;
        }
    }

    return sum;
}

///////     CONSTRUCCION DE CONFIGS     ///////////////////

/**
 * Se reinician las variables del reticulo y las tablas de precalculo, y se
 * actualizan los factores dependientes de beta y theta.
 */
void clear_lattice(void)
{
    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < Nx; ++x) {
            link_h[t][x] = 0;
            link_v[t][x] = 0;
            plaquette[t][x] = 0;
            monomer[t][x] = 0;
            fixed[t][x] = false;
            block[t][x] = false;
        }
    }

    for (int order = 0; order < NORDERS; ++order) {
        bessel_used_orders[order] = 0.0;
        bessel_used[order] = false;
    }

    eta = sqrt(beta / 2.0 - theta / (4.0 * PI));
    eta_conj = sqrt(beta / 2.0 + theta / (4.0 * PI));
}

/**
 * Se inserta un loop elemental orientado sobre una plaqueta y se actualizan
 * los sitios ocupados por la configuracion.
 *
 * @param t0 Coordenada temporal de la plaqueta.
 * @param x0 Coordenada espacial de la plaqueta.
 * @param dir Orientacion del loop, con valor 1 o -1.
 */
void set_loop_1x1(int t0, int x0, int dir)
{
    const int next_t = t0 + dt_p[t0];
    const int next_x = x0 + dx_p[x0];

    link_h[t0][x0] = dir;
    link_v[t0][next_x] = dir;
    link_h[next_t][x0] = -dir;
    link_v[t0][x0] = -dir;
    plaquette[t0][x0] += dir;

    block[t0][x0] = true;
    block[t0][next_x] = true;
    block[next_t][x0] = true;
    block[next_t][next_x] = true;

    monomer[t0][x0] = 0;
    monomer[t0][next_x] = 0;
    monomer[next_t][x0] = 0;
    monomer[next_t][next_x] = 0;
}

/**
 * Se distribuyen dimeros verticales a lo largo de una fila espacial.
 *
 * @param x0 Coordenada espacial inicial.
 * @param t0 Coordenada temporal de la fila.
 * @param freq Separacion espacial entre dimeros.
 */
void set_vstrip_spatial(int x0, int t0, int freq)
{
    for (int i = 0; i < Nx / freq; ++i) {
        const int x = (x0 + freq * i) % Nx;
        link_v[t0][x] = 2;
    }
}

/**
 * Se distribuyen dimeros horizontales a lo largo de una fila espacial.
 *
 * @param x0 Coordenada espacial inicial.
 * @param t0 Coordenada temporal de la fila.
 * @param freq Separacion espacial entre dimeros.
 */
void set_hstrip(int x0, int t0, int freq)
{
    for (int i = 0; i < Nx / freq; ++i) {
        const int x = (x0 + freq * i) % Nx;
        link_h[t0][x] = 2;
    }
}

/**
 * Se construye una configuracion inicial formada por dimeros horizontales.
 */
void dimmer_start(void)
{
    clear_lattice();

    for (int t = 0; t < Nt; ++t) {
        set_hstrip(0, t, 2);
    }
}

/**
 * Se construye una configuracion inicial con un monomero en cada sitio.
 */
void monomer_start(void)
{
    clear_lattice();

    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < Nx; ++x) {
            monomer[t][x] = 1;
        }
    }
}
