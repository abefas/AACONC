#ifndef AACO_misc_functions
#define AACO_misc_functions

#include <stdbool.h>
#include "structs.h"

double get_route_distance(SON *G, node *list);

double calculate_pheromone_sum(int idep, int ilast, node *v_cand, double *phMatrix, int n_nodes);

double calculate_h(SON *G, node* v_cand, int ilast);

double calculate_angle(int x, int y);

void find_polar_coordinates(SON *G, int vi, polar *polarArray);

void assign_vertices_to_sectors(SON *G, int vi, polar *polarArray, sector *sec, int n_sect);

bool not_empty(int *array, int length);

bool search_icluster(int *K_v, int icluster, int *v_free, int n_size);

int find_closest_free_vertex_in_sector(polar *polarArray, node *vertexList);

int find_closest_free_vertex(int *v_free, polar *polarArray, int numOfCustomers);

void find_free_in_icluster(node **head, int *v_free, int *K_v, int icluster, int n_size);

int errCheck(int value, int lower, int upper);

void new_best(vt_solution *new_r, vt_solution *r, SON *G);

int store_edge_count(vt_solution *Ra, SON *G, int *edge_matrix);

void update_pheromones(SON *G, double *phMatrix, vt_solution *R, vt_solution *R_best, double tupdate, double d);

double evaporate_pheromones(SON *G, int *edge_matrix, double *phMatrix, int edge_sum, int n_ants, double p_min, double p_max);

void fprint_results(vt_solution *R, SON *G);

void fprint_data(int iterations, int best_iter, double foundtime, double runtime);

#endif
