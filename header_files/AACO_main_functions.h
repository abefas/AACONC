#ifndef AACO_main_functions
#define AACO_main_functions

#include "structs.h"

void AACONC(SON *G, int n_ants, int n_freq, int n_size, int n_sect, int n_prim, 
                    double T_update, double a_update, double p_min, double p_max, double d, double a, double b);


void antSolution(SON *G, int **K, double *t, vt_solution *R, 
                    int n_size, int n_prim, double a, double b);

int selectVehicleType(int idepot, vt_solution *Ra, int *v_free, SON *G, int **K, double *phMatrix, 
                        int *launch_count, int n_size, int n_prim);

int selectDepot(vt_solution *Ra, int *v_free, SON *G, int **K, double *phMatrix, 
                int n_size, int n_prim);

int selectCluster(int idepot, int ilast, int *v_free, int *K, double *phMatrix,
                  SON *G, int n_size, int n_prim, double a, double b);

int selectCustomer(int idep, int ilast, node *v_cand, double *phMatrix, SON *G, 
                    double a, double b);


void createClusters(int *K, SON *G, int vi, int n_size, int n_prim, int n_sect);

#endif
