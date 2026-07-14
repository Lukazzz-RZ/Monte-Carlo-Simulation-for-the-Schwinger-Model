#include "Lattice.h"

///////     ACTUALIZACIONES LOCALES DE PLAQUETA     ///////////////////

/**
 * Se intercambia una configuracion local entre monomeros, dimeros paralelos
 * y un loop elemental.
 *
 * @param t0 Coordenada temporal de la plaqueta.
 * @param x0 Coordenada espacial de la plaqueta.
 * @param links Enlaces del contorno de la plaqueta.
 * @return true si se realiza la actualizacion y false en caso contrario.
 */
bool paralel_dimmers_1x1(int t0, int x0, int links[4])
{
    const int next_t = t0 + dt_p[t0];
    const int next_x = x0 + dx_p[x0];
    int local_mono[4] = {
        monomer[t0][x0],
        monomer[t0][next_x],
        monomer[next_t][next_x],
        monomer[next_t][x0]
    };

    int node_ini = -1;
    int dextro = 0;
    int chirality = 0;
    double w_ini = -1.0;
    double w_fin = -1.0;

    if (local_mono[0] + local_mono[1] +
        local_mono[2] + local_mono[3] == 4) {
        w_ini = 16.0 * m * m * m * m *
                comp_weight(plaquette[t0][x0]);
        node_ini = 0;
    } else if ((links[0] == 2 && links[2] == 2) ||
               (links[1] == 2 && links[3] == 2)) {
        w_ini = 2.0 * comp_weight(plaquette[t0][x0]);
        node_ini = 1;
    } else if (links[0] * links[1] * links[2] * links[3] == 1) {
        dextro = (links[0] == 1) ? 1 : -1;
        w_ini = 2.0 * comp_weight(plaquette[t0][x0]);
        node_ini = 2;
    }

    if (w_ini < 0.0 || node_ini < 0) {
        return false;
    }

    const int node_fin = (node_ini + 1 + rand_int(2)) % 3;

    if (node_fin == 0) {
        w_fin = 16.0 * m * m * m * m *
                comp_weight(plaquette[t0][x0] - dextro);
    } else if (node_fin == 1) {
        chirality = (rand_int(2) == 0) ? 1 : -1;
        w_fin = 2.0 * comp_weight(plaquette[t0][x0] - dextro);
    } else {
        chirality = (rand_int(2) == 0) ? 1 : -1;
        w_fin = 2.0 * comp_weight(plaquette[t0][x0] + chirality);
    }

    if (!metropolis(w_ini, w_fin)) {
        return false;
    }

    if (node_fin == 0) {
        block[t0][x0] = false;
        block[t0][next_x] = false;
        block[next_t][x0] = false;
        block[next_t][next_x] = false;

        monomer[t0][x0] = 1;
        monomer[t0][next_x] = 1;
        monomer[next_t][next_x] = 1;
        monomer[next_t][x0] = 1;

        plaquette[t0][x0] -= dextro;

        link_h[t0][x0] = 0;
        link_h[next_t][x0] = 0;
        link_v[t0][x0] = 0;
        link_v[t0][next_x] = 0;
    } else if (node_fin == 1) {
        if (chirality == 1) {
            link_h[t0][x0] = 2;
            link_h[next_t][x0] = 2;
            link_v[t0][x0] = 0;
            link_v[t0][next_x] = 0;
        } else {
            link_v[t0][x0] = 2;
            link_v[t0][next_x] = 2;
            link_h[t0][x0] = 0;
            link_h[next_t][x0] = 0;
        }

        plaquette[t0][x0] -= dextro;

        monomer[t0][x0] = 0;
        monomer[t0][next_x] = 0;
        monomer[next_t][next_x] = 0;
        monomer[next_t][x0] = 0;

        block[t0][x0] = false;
        block[t0][next_x] = false;
        block[next_t][x0] = false;
        block[next_t][next_x] = false;
    } else {
        set_loop_1x1(t0, x0, chirality);
    }

    return true;
}

/**
 * Se transforma una pareja de enlaces antiparalelos conservando la condicion
 * local de flujo.
 *
 * @param t0 Coordenada temporal de la plaqueta.
 * @param x0 Coordenada espacial de la plaqueta.
 * @param links Enlaces del contorno de la plaqueta.
 * @return true si se realiza la actualizacion y false en caso contrario.
 */
bool antiparalel_links(int t0, int x0, int links[4])
{
    const int next_t = t0 + dt_p[t0];
    const int prev_t = t0 + dt_n[t0];
    const int next_x = x0 + dx_p[x0];
    const int prev_x = x0 + dx_n[x0];

    const bool horizontal =
        links[0] * links[0] == 1 &&
        links[2] == -links[0] &&
        links[1] == 0 &&
        links[3] == 0;

    const bool vertical =
        links[1] * links[1] == 1 &&
        links[3] == -links[1] &&
        links[0] == 0 &&
        links[2] == 0;

    if (!horizontal && !vertical) {
        return false;
    }

    if (links[0] != 0) {
        int reference_link;

        if (link_h[t0][prev_x] != 0 && dx_n[x0] != 0) {
            reference_link = link_h[t0][prev_x];
        } else {
            reference_link = link_v[prev_t][x0];
        }

        if (!strong_coup_metropolis(
                plaquette[t0][x0],
                plaquette[t0][x0] - reference_link)) {
            return false;
        }

        link_h[t0][x0] = 0;
        link_h[next_t][x0] = 0;
        link_v[t0][x0] = reference_link;
        link_v[t0][next_x] = -reference_link;
        plaquette[t0][x0] -= reference_link;
    } else {
        int reference_link;

        if (link_v[prev_t][x0] != 0 && dt_n[t0] != 0) {
            reference_link = link_v[prev_t][x0];
        } else {
            reference_link = link_h[t0][prev_x];
        }

        if (!strong_coup_metropolis(
                plaquette[t0][x0],
                plaquette[t0][x0] + reference_link)) {
            return false;
        }

        link_h[t0][x0] = reference_link;
        link_h[next_t][x0] = -reference_link;
        link_v[t0][x0] = 0;
        link_v[t0][next_x] = 0;
        plaquette[t0][x0] += reference_link;
    }

    return true;
}

/**
 * Se intercambia una configuracion local entre monomeros con enlace, dimero
 * con enlace y una configuracion de tipo C.
 *
 * @param t0 Coordenada temporal de la plaqueta.
 * @param x0 Coordenada espacial de la plaqueta.
 * @param links Enlaces del contorno de la plaqueta.
 * @return true si se realiza la actualizacion y false en caso contrario.
 */
bool c_config(int t0, int x0, int links[4])
{
    const int next_t = t0 + dt_p[t0];
    const int next_x = x0 + dx_p[x0];

    for (int i = 0; i < 4; ++i) {
        int local_mono[4] = {
            monomer[t0][x0],
            monomer[t0][next_x],
            monomer[next_t][next_x],
            monomer[next_t][x0]
        };

        int node_ini = -1;
        int dextro = 0;
        double w_ini = -1.0;
        double w_fin = -1.0;

        if (local_mono[i] + local_mono[(i + 1) % 4] == 2 &&
            abs(links[(i + 2) % 4]) == 1) {
            w_ini = 4.0 * m * m *
                    comp_weight(plaquette[t0][x0]);
            node_ini = 0;
        } else if (links[i] == 2 &&
                   abs(links[(i + 2) % 4]) == 1) {
            w_ini = comp_weight(plaquette[t0][x0]);
            node_ini = 1;
        } else if (links[(i + 2) % 4] == 0 &&
                   abs(links[(i + 1) % 4]) == 1 &&
                   abs(links[i]) == 1 &&
                   abs(links[(i + 3) % 4]) == 1) {
            w_ini = comp_weight(plaquette[t0][x0]);
            node_ini = 2;
            dextro = links[i];

            if (i >= 2) {
                dextro = -dextro;
            }
        }

        if (node_ini < 0 || w_ini < 0.0) {
            continue;
        }

        const int node_fin = (node_ini + 1 + rand_int(2)) % 3;

        if (node_fin == 0) {
            w_fin = 4.0 * m * m *
                    comp_weight(plaquette[t0][x0] - dextro);
        } else if (node_fin == 1) {
            w_fin = comp_weight(plaquette[t0][x0] - dextro);
        } else {
            dextro = links[(i + 2) % 4];

            if (i >= 2) {
                dextro = -dextro;
            }

            w_fin = comp_weight(plaquette[t0][x0] + dextro);
        }

        if (!metropolis(w_ini, w_fin)) {
            return false;
        }

        if (node_fin == 2) {
            set_loop_1x1(t0, x0, dextro);

            if (i == 0) {
                link_h[next_t][x0] = 0;
            } else if (i == 1) {
                link_v[t0][x0] = 0;
            } else if (i == 2) {
                link_h[t0][x0] = 0;
            } else {
                link_v[t0][next_x] = 0;
            }

            if (node_ini == 0) {
                local_mono[i] = 0;
                local_mono[(i + 1) % 4] = 0;

                monomer[t0][x0] = local_mono[0];
                monomer[t0][next_x] = local_mono[1];
                monomer[next_t][next_x] = local_mono[2];
                monomer[next_t][x0] = local_mono[3];
            }

            return true;
        }

        plaquette[t0][x0] -= dextro;

        if (i == 2) {
            block[next_t][next_x] = false;
            block[next_t][x0] = false;

            if (node_ini == 2) {
                link_h[t0][x0] = -dextro;
            }

            if (node_fin == 1) {
                link_h[next_t][x0] = 2;
                monomer[next_t][next_x] = 0;
                monomer[next_t][x0] = 0;
            } else {
                monomer[next_t][next_x] = 1;
                monomer[next_t][x0] = 1;
                link_h[next_t][x0] = 0;
            }

            link_v[t0][x0] = 0;
            link_v[t0][next_x] = 0;
        } else if (i == 3) {
            block[t0][x0] = false;
            block[next_t][x0] = false;

            if (node_ini == 2) {
                link_v[t0][next_x] = -dextro;
            }

            if (node_fin == 1) {
                link_v[t0][x0] = 2;
                monomer[t0][x0] = 0;
                monomer[next_t][x0] = 0;
            } else {
                monomer[t0][x0] = 1;
                monomer[next_t][x0] = 1;
                link_v[t0][x0] = 0;
            }

            link_h[next_t][x0] = 0;
            link_h[t0][x0] = 0;
        } else if (i == 0) {
            block[t0][x0] = false;
            block[t0][next_x] = false;

            if (node_ini == 2) {
                link_h[next_t][x0] = dextro;
            }

            if (node_fin == 1) {
                link_h[t0][x0] = 2;
                monomer[t0][next_x] = 0;
                monomer[t0][x0] = 0;
            } else {
                monomer[t0][x0] = 1;
                monomer[t0][next_x] = 1;
                link_h[t0][x0] = 0;
            }

            link_v[t0][x0] = 0;
            link_v[t0][next_x] = 0;
        } else {
            block[t0][next_x] = false;
            block[next_t][next_x] = false;

            if (node_ini == 2) {
                link_v[t0][x0] = dextro;
            }

            if (node_fin == 1) {
                link_v[t0][next_x] = 2;
                monomer[next_t][next_x] = 0;
                monomer[t0][next_x] = 0;
            } else {
                monomer[next_t][next_x] = 1;
                monomer[t0][next_x] = 1;
                link_v[t0][next_x] = 0;
            }

            link_h[next_t][x0] = 0;
            link_h[t0][x0] = 0;
        }

        return true;
    }

    return false;
}

/**
 * Se desplaza un monomero entre esquinas opuestas de una plaqueta y se
 * reorganizan los enlaces asociados.
 *
 * @param t0 Coordenada temporal de la plaqueta.
 * @param x0 Coordenada espacial de la plaqueta.
 * @param links Enlaces del contorno de la plaqueta.
 * @return true si se realiza la actualizacion y false en caso contrario.
 */
bool flip_corner(int t0, int x0, int links[4])
{
    const int next_t = t0 + dt_p[t0];
    const int next_x = x0 + dx_p[x0];
    int local_mono[4] = {
        monomer[t0][x0],
        monomer[t0][next_x],
        monomer[next_t][next_x],
        monomer[next_t][x0]
    };

    for (int i = 0; i < 4; ++i) {
        const bool isolated_monomer =
            local_mono[i] == 1 &&
            !local_mono[(i + 1) % 4] &&
            !local_mono[(i + 2) % 4] &&
            !local_mono[(i + 3) % 4];

        const bool occupied_corner =
            abs(links[(i + 1) % 4] * links[(i + 2) % 4]) == 1;

        if (!isolated_monomer || !occupied_corner) {
            continue;
        }

        const int delta_t = (i >= 2) ? dt_p[t0] : dt_n[t0];
        const int out_t = (i < 2) ? dt_p[t0] : dt_n[t0];

        const int new_p =
            (plaquette[t0][x0] == plaquette[t0 + delta_t][x0])
                ? plaquette[t0 + out_t][x0]
                : plaquette[t0 + delta_t][x0];

        if (!strong_coup_metropolis(plaquette[t0][x0], new_p)) {
            return false;
        }

        local_mono[i] = 0;
        local_mono[(i + 2) % 4] = 1;

        links[i] = links[(i + 2) % 4];
        links[(i + 3) % 4] = links[(i + 1) % 4];
        links[(i + 1) % 4] = 0;
        links[(i + 2) % 4] = 0;

        monomer[t0][x0] = local_mono[0];
        monomer[t0][next_x] = local_mono[1];
        monomer[next_t][next_x] = local_mono[2];
        monomer[next_t][x0] = local_mono[3];

        link_h[t0][x0] = links[0];
        link_v[t0][next_x] = links[1];
        link_h[next_t][x0] = links[2];
        link_v[t0][x0] = links[3];

        plaquette[t0][x0] = new_p;
        sign *= -1;

        return true;
    }

    return false;
}

///////     ACTUALIZACION GLOBAL     ///////////////////

/**
 * Se propone desplazar simultaneamente todos los numeros de plaqueta en una
 * unidad con orientacion aleatoria.
 *
 * @return true si se acepta la actualizacion y false en caso contrario.
 */
bool global_jump(void)
{
    const int sigma = (rand_int(2) == 0) ? 1 : -1;
    double log_probability = 0.0;

    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < Nx; ++x) {
            const int occupation = plaquette[t][x];

            log_probability +=
                log(comp_weight(occupation + sigma)) -
                log(comp_weight(occupation));
        }
    }

    if (log(rand_double()) >= log_probability) {
        return false;
    }

    for (int t = 0; t < Nt; ++t) {
        for (int x = 0; x < Nx; ++x) {
            plaquette[t][x] += sigma;
        }
    }

    return true;
}

///////     SELECCION DE PLAQUETA     ///////////////////

/**
 * Se intenta actualizar una plaqueta seleccionada mediante los movimientos
 * locales disponibles.
 *
 * @param t0 Coordenada temporal de la plaqueta.
 * @param x0 Coordenada espacial de la plaqueta.
 * @return true si se realiza o bloquea la actualizacion y false si no existe
 * una propuesta admisible.
 */
bool selected_plaq_update(int t0, int x0)
{
    const int next_t = t0 + dt_p[t0];
    const int next_x = x0 + dx_p[x0];
    const int prev_x = x0 + dx_n[x0];

    if (fixed[t0][x0]) {
        return false;
    }

    int links[4] = {
        link_h[t0][x0],
        link_v[t0][next_x],
        link_h[next_t][x0],
        link_v[t0][x0]
    };

    if (fixed[t0][prev_x]) {
        return true;
    }

    return paralel_dimmers_1x1(t0, x0, links) ||
           flip_corner(t0, x0, links) ||
           antiparalel_links(t0, x0, links) ||
           c_config(t0, x0, links);
}

/**
 * Se selecciona una plaqueta aleatoria y se intenta aplicar una actualizacion
 * local.
 *
 * @return true si se realiza o bloquea la actualizacion y false si no existe
 * una propuesta admisible.
 */
bool plaq_update(void)
{
    const int t0 = rand() % Nt;
    const int x0 = rand() % (Nx - CONTORNO_X);

    return selected_plaq_update(t0, x0);
}
