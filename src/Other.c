#include "Lattice.h"

///////     GENERACION ALEATORIA     ///////////////////

/**
 * Se genera un numero real pseudoaleatorio en el intervalo cerrado entre cero
 * y uno.
 *
 * @return Numero real pseudoaleatorio.
 */
double rand_double(void)
{
    unsigned int random_value;

    rand_s(&random_value);

    return (double)random_value / (double)UINT_MAX;
}

/**
 * Se genera un numero entero pseudoaleatorio entre cero y N menos uno.
 *
 * @param N Limite superior no incluido.
 * @return Numero entero pseudoaleatorio.
 */
int rand_int(int N)
{
    double random_value = rand_double();

    if (random_value == 1.0) {
        random_value = 0.0;
    }

    return (int)(random_value * N);
}

///////     COMPROBACION DEL RETICULO     ///////////////////

/**
 * Se comprueba la condicion de saturacion fermionica en un sitio.
 *
 * @param t Coordenada temporal del sitio.
 * @param x Coordenada espacial del sitio.
 * @return true si el sitio cumple la condicion y false en caso contrario.
 */
bool check_site(int t, int x)
{
    int connections = 0;

    if (dt_n[t] != 0) {
        const int previous_t = t + dt_n[t];

        connections += abs(link_v[previous_t][x]);
    }

    if (dx_n[x] != 0) {
        const int previous_x = x + dx_n[x];

        connections += abs(link_h[t][previous_x]);
    }

    connections += abs(link_h[t][x]);
    connections += abs(link_v[t][x]);

    if (monomer[t][x] == 1) {
        return connections == 0;
    }

    return connections == 2;
}

/**
 * Se comprueba la compatibilidad local entre los numeros de plaqueta y los
 * enlaces orientados.
 *
 * @param t Coordenada temporal de la plaqueta.
 * @param x Coordenada espacial de la plaqueta.
 * @return true si se cumple la condicion de flujo y false en caso contrario.
 */
bool check_plaq(int t, int x)
{
    if (dt_n[t] == 0 || dx_n[x] == 0) {
        return true;
    }

    const int current_p = plaquette[t][x];
    const int previous_t_p = plaquette[t + dt_n[t]][x];
    const int previous_x_p = plaquette[t][x + dx_n[x]];

    const int horizontal_link =
        (link_h[t][x] != 2) ? link_h[t][x] : 0;

    const int vertical_link =
        (link_v[t][x] != 2) ? link_v[t][x] : 0;

    return current_p - previous_t_p == horizontal_link &&
           previous_x_p - current_p == vertical_link;
}

/**
 * Se comprueban las condiciones locales de saturacion y flujo en todo el
 * reticulo.
 *
 * @return true si la configuracion es valida y false en caso contrario.
 */
bool check_lattice(void)
{
    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < Nx - CONTORNO_X; ++x) {
            if (!check_site(t, x) || !check_plaq(t, x)) {
                show_lattice();

                printf(
                    "Fallo de saturado en t: %d, x: %d\n",
                    t,
                    x
                );

                Sleep(10000);

                return false;
            }
        }
    }

    return true;
}

///////     RELACION ENTRE LOOPS     ///////////////////

/**
 * Se comprueba si dos enlaces pertenecen al mismo loop orientado.
 *
 * @param t1 Coordenada temporal del primer enlace.
 * @param x1 Coordenada espacial del primer enlace.
 * @param h1 Indica si el primer enlace es horizontal.
 * @param t2 Coordenada temporal del segundo enlace.
 * @param x2 Coordenada espacial del segundo enlace.
 * @param h2 Indica si el segundo enlace es horizontal.
 * @return true si ambos enlaces pertenecen al mismo loop y false en caso
 * contrario.
 */
bool same_loop(
    int t1,
    int x1,
    bool h1,
    int t2,
    int x2,
    bool h2
)
{
    int current_t = t1;
    int current_x = x1;

    if (link_h[current_t][current_x] == 1) {
        current_x = (current_x + 1) % Nx;
    } else if (link_v[(current_t - 1 + Nt) % Nt][current_x] == -1) {
        current_t = (current_t - 1 + Nt) % Nt;
    } else if (link_h[current_t][(current_x - 1 + Nx) % Nx] == -1) {
        current_x = (current_x - 1 + Nx) % Nx;
    } else if (link_v[current_t][current_x] == 1) {
        current_t = (current_t + 1) % Nt;
    }

    for (int step = 0; step < 2 * Nt * Nx; ++step) {
        if (h2) {
            if (link_h[current_t][current_x] == link_h[t2][x2]) {
                return true;
            }
        } else if (link_v[current_t][current_x] == link_v[t2][x2]) {
            return true;
        }

        if (h1) {
            if (link_h[current_t][current_x] == link_h[t1][x1]) {
                return false;
            }
        } else if (link_v[current_t][current_x] == link_v[t1][x1]) {
            return false;
        }

        if (link_h[current_t][current_x] == 1) {
            current_x = (current_x + 1) % Nx;
        } else if (
            link_v[(current_t - 1 + Nt) % Nt][current_x] == -1
        ) {
            current_t = (current_t - 1 + Nt) % Nt;
        } else if (
            link_h[current_t][(current_x - 1 + Nx) % Nx] == -1
        ) {
            current_x = (current_x - 1 + Nx) % Nx;
        } else if (link_v[current_t][current_x] == 1) {
            current_t = (current_t + 1) % Nt;
        }
    }

    return false;
}

/**
 * Se comprueba si dos plaquetas estan relacionadas por un mismo loop o por
 * una linea de enlaces orientados que atraviesa el reticulo.
 *
 * @param t1 Coordenada temporal de la primera plaqueta.
 * @param x1 Coordenada espacial de la primera plaqueta.
 * @param t2 Coordenada temporal de la segunda plaqueta.
 * @param x2 Coordenada espacial de la segunda plaqueta.
 * @param h Indica si la comprobacion se realiza en direccion horizontal.
 * @return true si las plaquetas estan relacionadas y false en caso contrario.
 */
bool plaq_related(
    int t1,
    int x1,
    int t2,
    int x2,
    bool h
)
{
    if (same_loop(t1, x1, h, t2, x2, h)) {
        return true;
    }

    if (h) {
        for (int offset = 0; offset < Nx; ++offset) {
            const int first_link = link_h[t1][(x1 + offset) % Nx];
            const int second_link = link_h[t2][(x2 + offset) % Nx];

            if (first_link == 2 || first_link == 0) {
                return false;
            }

            if (second_link == 2 || second_link == 0) {
                return false;
            }
        }
    } else {
        for (int offset = 0; offset < Nt; ++offset) {
            const int first_link = link_v[(t1 + offset) % Nt][x1];
            const int second_link = link_v[(t2 + offset) % Nt][x2];

            if (first_link == 2 || first_link == 0) {
                return false;
            }

            if (second_link == 2 || second_link == 0) {
                return false;
            }
        }
    }

    return true;
}
