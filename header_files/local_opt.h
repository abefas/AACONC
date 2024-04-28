#ifndef local_opt
#define local_opt 
#include "structs.h"

double local_opt_full(vt_solution *R, SON *G);

double k_optimization2(route *r, SON *G, int n_max);

double mutual_k_optimization(vt_solution *R, SON *G, int condition, int n_max);

double mutual_drone(vt_solution *R, SON *G);

double depot_VT_optimization(vt_solution *R, SON *G, int n_max);

void remove_duplicate_nodes(node **list);

int check_route_feasibility(node *route, SON *G);

#endif
