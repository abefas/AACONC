#include <stdio.h>
#include <stdlib.h>
#include "header_files/listFunctions.h"
#include "header_files/structs.h"
#include "header_files/AACO_main_functions.h"
#include "header_files/AACO_misc_functions.h"

void antSolution(SON *G, int **K, double *phMatrix, vt_solution *Ra, int n_size, int n_prim, double a, double b){

    /* Initialize v_free array which contains not yet visited Customer vertices */
    int v_free[G->n_customers];
    for(int i = 0; i < G->n_customers; i++)
        v_free[i] = 1; //free = 1, not free = -1

    int idepot, ilast, icluster, icustomer;

    node *v_candidates = NULL;

    Ra->total_distance = 0.0;
    for(int idep = 0; idep < G->n_depots; idep++){
        deleteList(&Ra->a_depots[idep].routelist);
        Ra->a_depots[idep].distance = 0.0;
        Ra->a_depots[idep].quantity_served = 0;
        Ra->a_depots[idep].v_d = Ra->a_depots[idep].depot_id;
        Ra->a_depots[idep].current_load = G->capacity;
        push(&Ra->a_depots[idep].routelist, Ra->a_depots[idep].depot_id);
    }

    /*********** Main loop ***********/

    while(not_empty(v_free, G->n_customers))
    { 
        idepot = (icluster = (icustomer = -1));

        while(idepot == -1 || icluster == -1 || icustomer == -1){

            idepot = selectDepot(Ra, v_free, G, K, phMatrix, n_size, n_prim);
            if(idepot == -1){
                printf("idepot -1");
                continue;
            }


            ilast = Ra->a_depots[idepot].v_d - 1;
            if(errCheck(ilast, 0, G->n_nodes-1) == 1){ 
                perror("ilast error in antSolution.c\n"); 
                exit(1); 
            }

            icluster = selectCluster(idepot, ilast, v_free, K[ilast], phMatrix, G, n_size, n_prim, a, b);
            if(icluster == -1){
                printf("icluster -1\n");
                continue;
            }

            find_free_in_icluster(&v_candidates, v_free, K[ilast], icluster, n_size); 
            if(!v_candidates){
                perror("Could not find free in cluster k in function antSolution\n");
                exit(1);
            }

            icustomer = selectCustomer(idepot, ilast, v_candidates, phMatrix, G, a, b);
            if(icustomer >= G->n_customers || icustomer < 0){
                printf("idepot %d icluster %d icustomer %d\n", idepot, icluster, icustomer);
                exit(1);
            }

            deleteList(&v_candidates);
        }

        //Vehicle returns to depot to reload
        if(Ra->a_depots[idepot].current_load < G->a_customers[icustomer].demand){

            append(&Ra->a_depots[idepot].routelist, Ra->a_depots[idepot].depot_id);

            Ra->a_depots[idepot].current_load = G->capacity;

            Ra->a_depots[idepot].distance += G->d_matrix[ilast][Ra->a_depots[idepot].depot_id - 1];

            Ra->total_distance += G->d_matrix[ilast][Ra->a_depots[idepot].depot_id - 1];

            ilast = Ra->a_depots[idepot].depot_id - 1;

        }

        /* Insert selected customer to route */
        append(&Ra->a_depots[idepot].routelist, icustomer+1);

        Ra->a_depots[idepot].distance += G->d_matrix[ilast][icustomer];
        Ra->a_depots[idepot].v_d = icustomer + 1;
        Ra->a_depots[idepot].current_load -= G->a_customers[icustomer].demand;
        Ra->a_depots[idepot].quantity_served += G->a_customers[icustomer].demand;

        Ra->total_distance += G->d_matrix[ilast][icustomer];


        v_free[icustomer] = -1;      //Mark as visited
    }

    /* Vehicles return to depot */
    int last, depotID;
    for(int idep = 0; idep < G->n_depots; idep++){
        if(Ra->a_depots[idep].v_d != Ra->a_depots[idep].depot_id){
            append(&Ra->a_depots[idep].routelist, Ra->a_depots[idep].depot_id);
            Ra->a_depots[idep].distance += G->d_matrix[Ra->a_depots[idep].v_d - 1][Ra->a_depots[idep].depot_id - 1];
            Ra->total_distance += G->d_matrix[Ra->a_depots[idep].v_d - 1][Ra->a_depots[idep].depot_id - 1];
        }
    }

    return;
}

