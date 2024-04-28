#ifndef structs
#define structs
#include "listFunctions.h"

#define epsilon 0.000001

typedef struct Customer{
    int id;
    int x;
    int y;
    int demand;
} customer;

typedef struct Combined{
    int id;
    int x;
    int y;
} combined;

typedef struct Depot{
    int id;
    int x;
    int y;
} depot;

typedef struct SetOfNodes{
    int capacity;
    int n_customers;
    int n_depots;
    int n_nodes;
    double **d_matrix;
    customer *a_customers;
    depot *a_depots;
    combined *a_combined;
} SON;

typedef struct Route{
    int v_d;
    int current_load;
    double distance;
    int quantity_served;
    int depot_id;
    node *routelist;
} route;

typedef struct vt_solution{
    double total_distance;
    route *a_depots;
} vt_solution;


//Sectoring structs
typedef struct setOfClusters{
    int *array;
} setOfClusters;

typedef struct polarCoordinates{
    double distance;
    double degrees;
} polar;

typedef struct sector{
    node *vertexList;
} sector;

#endif
