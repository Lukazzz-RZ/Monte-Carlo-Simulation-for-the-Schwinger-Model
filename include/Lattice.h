#pragma once

#define _CRT_SECURE_NO_WARNINGS

///////     LIBRERIAS     ///////////////////

#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h>
#include <windows.h>

#define popen _popen
#define pclose _pclose
#endif

///////     PARAMETROS DEL RETICULO     ///////////////////

#define Nx 10
#define Nt 10
#define CONTORNO_X 0
#define NORDERS 100

///////     CONSTANTES NUMERICAS     ///////////////////

#define PI 3.14159265358979323846264338327950288419716939937510
#define N_OBS 10

///////     TIPOS DE DATOS     ///////////////////

typedef struct {
    double value;
    double error;
} medida;

enum {
    OBS_RHO_D = 0,
    OBS_RHO_D_SQR,
    OBS_Q,
    OBS_Q_SQR,
    OBS_U_P,
    OBS_U_P_SQR,
    OBS_RHO_M,
    OBS_RHO_M_SQR,
    OBS_CHIRAL,
    OBS_CHIRAL_SQR
};

///////     PARAMETROS GLOBALES     ///////////////////

extern double beta;
extern double m;
extern double theta;
extern double eta;
extern double eta_conj;
extern double sign;

extern int K;

extern bool WORM_DEBUG;
extern bool WORMCORE_DEBUG;

///////     VARIABLES DEL RETICULO     ///////////////////

extern int link_h[Nt][Nx];
extern int link_v[Nt][Nx];
extern int plaquette[Nt][Nx];
extern int monomer[Nt][Nx];

extern bool block[Nt][Nx];
extern bool fixed[Nt][Nx];

///////     VARIABLES DE PRECALCULO     ///////////////////

extern double bessel_used_orders[NORDERS];
extern bool bessel_used[NORDERS];

extern int dx_p[Nx];
extern int dx_n[Nx];
extern int dt_p[Nt];
extern int dt_n[Nt];

///////     INICIALIZACION Y PRECALCULOS     ///////////////////

/**
 * Se inicializan los desplazamientos entre sitios vecinos.
 */
void init_moves(
    int *dx_plus,
    int *dx_minus,
    int *dt_plus,
    int *dt_minus
);

/**
 * Se reinicia el reticulo y se actualizan los factores de precalculo.
 */
void clear_lattice(void);

double besselmod(int n, double x);
void set_loop_1x1(int t0, int x0, int dir);
void set_vstrip_spatial(int x0, int t0, int freq);
void set_hstrip(int x0, int t0, int freq);
void dimmer_start(void);
void monomer_start(void);

///////     NUCLEO DE ACTUALIZACION     ///////////////////

double comp_weight(int p);
bool strong_coup_metropolis(int p_ini, int p_fin);
bool metropolis(double p_ini, double p_fin);
void plaquette_sweeep(void);

/**
 * Se realiza un paso Monte Carlo completo.
 */
void monte_carlo_step(void);

/**
 * Se termaliza la configuracion mediante pasos Monte Carlo completos.
 */
void thermalizate(int n_steps, const char *spec);

///////     ACTUALIZACIONES DE PLAQUETA     ///////////////////

bool paralel_dimmers_1x1(int t0, int x0, int links[4]);
bool antiparalel_links(int t0, int x0, int links[4]);
bool c_config(int t0, int x0, int links[4]);
bool flip_corner(int t0, int x0, int links[4]);
bool global_jump(void);
bool selected_plaq_update(int t0, int x0);

/**
 * Se selecciona una plaqueta aleatoria y se intenta actualizar.
 */
bool plaq_update(void);

///////     ACTUALIZACIONES DE GUSANO     ///////////////////

double clamp_prob(double p);
bool start_worm(int *t0, int *x0);
bool choose_worm_mode(
    int t0,
    int x0,
    double p_add,
    double p_remove
);
int explorar_vecinos(
    int t,
    int x,
    bool next_dimmer,
    int disponibles[4]
);
bool can_stop_mono(
    bool worm_con_monomero,
    bool next_dimmer,
    int t,
    int x
);
bool try_to_stop(
    int t,
    int x,
    double p_add,
    double p_remove
);
void move_worm(int dir, int *t, int *x);
void upd_current_path(bool *next_dimmer, int t, int x);
int worm_decision(int n_dir, int disps[4], bool can_stop);
int choose_end_to_move(int n_head, int n_tail);
bool run_open_worm(
    int t0,
    int x0,
    double p_add,
    double p_remove
);
bool run_closed_worm_head_tail(int t0, int x0);

/**
 * Se ejecuta una actualizacion de gusano abierto o cerrado.
 */
bool run_worm(void);

///////     MEDIDAS     ///////////////////

double besselmod_safe(int n, double x);

double measure_rho_d(void);
double measure_rho_d_sqr(void);

double local_q_exact(int p);
double measure_q(void);
double measure_q_sqr(void);

double local_U_p_exact(int p);
double measure_U_p(void);
double measure_U_p_sqr(void);

double measure_monomer_density(void);
double measure_monomer_density_sqr(void);

/**
 * Se calculan los observables principales de la configuracion actual.
 */
void measure_all(double *obs);

///////     VISUALIZACION Y RESULTADOS     ///////////////////

void show_lattice(void);
void show_probs(void);

/**
 * Se escriben los observables de una configuracion en un archivo de texto.
 */
void write_results(
    bool overwrite,
    double *obs,
    int n_obs,
    const char *spec
);

void write_thermal(
    int step,
    double U_p,
    double rho_d,
    const char *spec
);

void plot_thermal(const char *spec);
void plot_obs(const char *spec);

///////     FUNCIONES AUXILIARES     ///////////////////

double rand_double(void);
int rand_int(int N);

bool check_site(int t, int x);
bool check_plaq(int t, int x);

/**
 * Se comprueba la saturacion fermionica y la conservacion local del flujo.
 */
bool check_lattice(void);

bool same_loop(
    int t1,
    int x1,
    bool h1,
    int t2,
    int x2,
    bool h2
);

bool plaq_related(
    int t1,
    int x1,
    int t2,
    int x2,
    bool h
);
