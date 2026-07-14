#include "Lattice.h"

///////     CONTROL DE DEPURACION     ///////////////////

bool WORM_DEBUG = false;

///////     CONSTANTES DEL GUSANO     ///////////////////

enum { WORM_STOP = 91 };

///////     PROBABILIDADES     ///////////////////

/**
 * Se restringe una probabilidad al intervalo cerrado entre cero y uno.
 *
 * @param p Probabilidad de entrada.
 * @return Probabilidad restringida al intervalo permitido.
 */
double clamp_prob(double p)
{
    if (p < 0.0) {
        return 0.0;
    }

    if (p > 1.0) {
        return 1.0;
    }

    return p;
}

/**
 * Se selecciona aleatoriamente un sitio no bloqueado para iniciar el gusano.
 *
 * @param t0 Coordenada temporal del sitio seleccionado.
 * @param x0 Coordenada espacial del sitio seleccionado.
 * @return true si existe un sitio disponible y false en caso contrario.
 */
bool start_worm(int *t0, int *x0)
{
    int free_sites = 0;

    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < Nx; ++x) {
            if (!block[t][x]) {
                ++free_sites;
            }
        }
    }

    if (free_sites == 0) {
        return false;
    }

    const int target = rand_int(free_sites);
    int counter = 0;

    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < Nx; ++x) {
            if (block[t][x]) {
                continue;
            }

            if (counter == target) {
                *t0 = t;
                *x0 = x;
                return true;
            }

            ++counter;
        }
    }

    return false;
}

/**
 * Se determina si el gusano debe operar en el sector abierto o cerrado.
 *
 * @param t0 Coordenada temporal del sitio inicial.
 * @param x0 Coordenada espacial del sitio inicial.
 * @param p_add Probabilidad de anadir un monomero.
 * @param p_remove Probabilidad de eliminar un monomero.
 * @return true para ejecutar un gusano abierto y false para uno cerrado.
 */
bool choose_worm_mode(int t0, int x0, double p_add, double p_remove)
{
    const double probability =
        (monomer[t0][x0] == 1) ? p_remove : p_add;

    return rand_double() < probability;
}

///////     EXPLORACION DEL RETICULO     ///////////////////

/**
 * Se identifican las direcciones disponibles desde la posicion actual.
 *
 * @param t Coordenada temporal actual.
 * @param x Coordenada espacial actual.
 * @param next_dimmer Estado requerido para el siguiente enlace.
 * @param disponibles Direcciones disponibles encontradas.
 * @return Numero de direcciones disponibles.
 */
int explorar_vecinos(int t, int x, bool next_dimmer, int disponibles[4])
{
    int n_available = 0;

    if (dx_p[x] != 0 && !block[t][x + dx_p[x]]) {
        const int link = link_h[t][x];

        if (next_dimmer ? link == 0 : link == 2) {
            disponibles[n_available++] = 0;
        }
    }

    if (dx_n[x] != 0 && !block[t][x + dx_n[x]]) {
        const int link = link_h[t][x + dx_n[x]];

        if (next_dimmer ? link == 0 : link == 2) {
            disponibles[n_available++] = 1;
        }
    }

    if (dt_p[t] != 0 && !block[t + dt_p[t]][x]) {
        const int link = link_v[t][x];

        if (next_dimmer ? link == 0 : link == 2) {
            disponibles[n_available++] = 2;
        }
    }

    if (dt_n[t] != 0 && !block[t + dt_n[t]][x]) {
        const int link = link_v[t + dt_n[t]][x];

        if (next_dimmer ? link == 0 : link == 2) {
            disponibles[n_available++] = 3;
        }
    }

    if (WORM_DEBUG) {
        printf("\nVecinos disponibles: %d\n", n_available);
    }

    return n_available;
}

/**
 * Se comprueba si un gusano abierto puede detenerse en el sitio actual.
 *
 * @param worm_con_monomero Indica si se ejecuta un gusano abierto.
 * @param next_dimmer Estado requerido para el siguiente enlace.
 * @param t Coordenada temporal actual.
 * @param x Coordenada espacial actual.
 * @return true si el cierre es admisible y false en caso contrario.
 */
bool can_stop_mono(
    bool worm_con_monomero,
    bool next_dimmer,
    int t,
    int x
)
{
    if (!worm_con_monomero) {
        return false;
    }

    return (!next_dimmer && monomer[t][x] == 1) ||
           (next_dimmer && monomer[t][x] == 0);
}

/**
 * Se intenta cerrar un gusano abierto modificando el monomero del sitio.
 *
 * @param t Coordenada temporal actual.
 * @param x Coordenada espacial actual.
 * @param p_add Probabilidad de anadir un monomero.
 * @param p_remove Probabilidad de eliminar un monomero.
 * @return true si se acepta el cierre y false en caso contrario.
 */
bool try_to_stop(int t, int x, double p_add, double p_remove)
{
    const double random_value = rand_double();

    if (monomer[t][x] == 0) {
        if (random_value < p_add) {
            monomer[t][x] = 1;
            return true;
        }
    } else if (random_value < p_remove) {
        monomer[t][x] = 0;
        return true;
    }

    return false;
}

///////     MOVIMIENTO DEL GUSANO     ///////////////////

/**
 * Se desplaza un extremo del gusano y se alterna la ocupacion del enlace
 * recorrido.
 *
 * @param dir Direccion del desplazamiento.
 * @param t Coordenada temporal del extremo.
 * @param x Coordenada espacial del extremo.
 */
void move_worm(int dir, int *t, int *x)
{
    if (dir == 0) {
        link_h[*t][*x] = 2 - link_h[*t][*x];
        *x += dx_p[*x];
    } else if (dir == 1) {
        const int previous_x = *x + dx_n[*x];

        link_h[*t][previous_x] = 2 - link_h[*t][previous_x];
        *x = previous_x;
    } else if (dir == 2) {
        link_v[*t][*x] = 2 - link_v[*t][*x];
        *t += dt_p[*t];
    } else if (dir == 3) {
        const int previous_t = *t + dt_n[*t];

        link_v[previous_t][*x] = 2 - link_v[previous_t][*x];
        *t = previous_t;
    }
}

/**
 * Se actualiza el tipo de enlace requerido tras alcanzar un nuevo sitio.
 *
 * @param next_dimmer Estado requerido para el siguiente enlace.
 * @param t Coordenada temporal actual.
 * @param x Coordenada espacial actual.
 */
void upd_current_path(bool *next_dimmer, int t, int x)
{
    if (monomer[t][x] == 1) {
        *next_dimmer = false;
    } else {
        *next_dimmer = !*next_dimmer;
    }
}

/**
 * Se selecciona una direccion disponible o la opcion de detener el gusano.
 *
 * @param n_dir Numero de direcciones disponibles.
 * @param disps Direcciones disponibles.
 * @param can_stop Indica si puede seleccionarse la opcion de parada.
 * @return Direccion seleccionada, WORM_STOP o -1 si no hay movimientos.
 */
int worm_decision(int n_dir, int disps[4], bool can_stop)
{
    const int moves = can_stop ? n_dir + 1 : n_dir;

    if (moves == 0) {
        return -1;
    }

    const int decision = rand_int(moves);

    if (can_stop && decision == n_dir) {
        return WORM_STOP;
    }

    return disps[decision];
}

/**
 * Se selecciona el extremo del gusano cerrado que debe desplazarse.
 *
 * @param n_head Movimientos disponibles para la cabeza.
 * @param n_tail Movimientos disponibles para la cola.
 * @return Cero para la cabeza, uno para la cola o -1 si no pueden moverse.
 */
int choose_end_to_move(int n_head, int n_tail)
{
    if (n_head > 0 && n_tail > 0) {
        return rand_int(2);
    }

    if (n_head > 0) {
        return 0;
    }

    if (n_tail > 0) {
        return 1;
    }

    return -1;
}

///////     GUSANO ABIERTO     ///////////////////

/**
 * Se ejecuta un gusano abierto que conecta dos modificaciones de monomeros
 * mediante una trayectoria alternada de enlaces.
 *
 * @param t0 Coordenada temporal inicial.
 * @param x0 Coordenada espacial inicial.
 * @param p_add Probabilidad de anadir un monomero.
 * @param p_remove Probabilidad de eliminar un monomero.
 * @return true si el gusano se cierra y false si queda bloqueado.
 */
bool run_open_worm(int t0, int x0, double p_add, double p_remove)
{
    int t = t0;
    int x = x0;
    int disps[4] = {0, 0, 0, 0};

    monomer[t0][x0] = 1 - monomer[t0][x0];

    bool next_dimmer = monomer[t0][x0] == 0;

    while (true) {
        const int n_dir =
            explorar_vecinos(t, x, next_dimmer, disps);

        const bool can_stop =
            can_stop_mono(true, next_dimmer, t, x);

        if (n_dir == 0 && !can_stop) {
            return false;
        }

        const int action =
            worm_decision(n_dir, disps, can_stop);

        if (action == WORM_STOP) {
            if (try_to_stop(t, x, p_add, p_remove)) {
                return true;
            }
        } else if (action >= 0 && action <= 3) {
            move_worm(action, &t, &x);
            upd_current_path(&next_dimmer, t, x);
        } else {
            return false;
        }
    }
}

///////     GUSANO CERRADO     ///////////////////

/**
 * Se ejecuta un gusano cerrado desplazando de forma alternada la cabeza y la
 * cola hasta que ambos extremos vuelven a coincidir.
 *
 * @param t0 Coordenada temporal inicial.
 * @param x0 Coordenada espacial inicial.
 * @return true si el gusano se cierra y false si ambos extremos se bloquean.
 */
bool run_closed_worm_head_tail(int t0, int x0)
{
    int head_t = t0;
    int head_x = x0;
    int tail_t = t0;
    int tail_x = x0;

    bool head_next_dimmer = false;
    bool tail_next_dimmer = false;

    int disps_head[4] = {0, 0, 0, 0};
    int disps_tail[4] = {0, 0, 0, 0};
    int step = 0;

    while (true) {
        const int n_head = explorar_vecinos(
            head_t,
            head_x,
            head_next_dimmer,
            disps_head
        );

        const int n_tail = explorar_vecinos(
            tail_t,
            tail_x,
            tail_next_dimmer,
            disps_tail
        );

        ++step;

        const int which_end =
            choose_end_to_move(n_head, n_tail);

        if (which_end < 0) {
            return false;
        }

        if (which_end == 0) {
            const int action =
                worm_decision(n_head, disps_head, false);

            if (action < 0 || action > 3) {
                return false;
            }

            move_worm(action, &head_t, &head_x);
            upd_current_path(
                &head_next_dimmer,
                head_t,
                head_x
            );
        } else {
            const int action =
                worm_decision(n_tail, disps_tail, false);

            if (action < 0 || action > 3) {
                return false;
            }

            move_worm(action, &tail_t, &tail_x);
            upd_current_path(
                &tail_next_dimmer,
                tail_t,
                tail_x
            );
        }

        if (head_t == tail_t &&
            head_x == tail_x &&
            step > 0) {
            if (WORM_DEBUG) {
                printf("Steps totales: %d\n", step);
            }

            return true;
        }
    }
}

///////     EJECUCION DEL GUSANO     ///////////////////

/**
 * Se selecciona el sitio inicial y se ejecuta un gusano abierto o cerrado
 * segun la ocupacion local y la masa fermionica.
 *
 * @return true si la actualizacion termina correctamente y false en caso
 * contrario.
 */
bool run_worm(void)
{
    const double p_add =
        (m != 0.0) ? clamp_prob(2.0 * m) : 0.0;

    const double p_remove =
        (m != 0.0) ? clamp_prob(1.0 / (2.0 * m)) : 0.0;

    int t0;
    int x0;

    if (!start_worm(&t0, &x0)) {
        return false;
    }

    const bool worm_con_monomero =
        choose_worm_mode(t0, x0, p_add, p_remove);

    if (worm_con_monomero) {
        return run_open_worm(t0, x0, p_add, p_remove);
    }

    return run_closed_worm_head_tail(t0, x0);
}
