#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include "header_files/AACO_misc_functions.h"
#include "header_files/listFunctions.h"
#include "header_files/structs.h"

extern int instance_id;

double get_route_distance(SON *G, node *list){
    if(!list) return 0;
    double distance = 0.0;
    node *temp = list;
    while(temp && temp->next){
        distance += G->d_matrix[temp->data - 1][temp->next->data - 1];
        temp = temp->next;
    }
    return distance;
}

double calculate_pheromone_sum(int idep, int ilast, node *v_cand, double *phMatrix, int n_nodes){
    if(!v_cand){ perror("v_cand == NULL at calculate_pheromone_sum\n"); exit(1); }
    if(!phMatrix){ perror("phMatrix == NULL at calculate_pheromone_sum\n"); exit(1); }
    double sum = 0.0;
    node *temp = v_cand;

    while(temp){
        if(!phMatrix[idep*n_nodes*n_nodes + ilast*(n_nodes) + (temp->data-1)])
        {
            perror("phMatrix NULL!\n");
            exit(1);
        }

        if(phMatrix[idep*n_nodes*n_nodes + ilast*(n_nodes) + (temp->data-1)] == -1.0)
        {
            perror("calculate_pheromone_sum error!\n");
            exit(1);
        }

        sum += phMatrix[idep*n_nodes*n_nodes + ilast*(n_nodes) + (temp->data-1)];

        temp = temp->next;
    }

    return sum;
}

double calculate_h(SON *G, node *v_cand, int ilast){
    if(!G){ perror("G == NULL at calculate_h\n"); exit(1); }
    if(!v_cand){ perror("v_cand == NULL at calculate_h\n"); exit(1); }
    double distance = 0.0;
    struct Node *temp = v_cand;
    while(temp){
        distance += G->d_matrix[ilast][temp->data - 1];
        temp = temp->next;
    }
    return 1.0/distance;
}

//Calculate arctan (inverse tangent)
double calculate_angle(int x, int y){
    return atan2(y,x) * (180.0/M_PI);   //* (180.0/M_PI) is used to transform radians to degrees
}

void find_polar_coordinates(SON *G, int vi, polar *polarArray){
    if(!polarArray){ perror("polarArray == NULL at find_polar_coordinates\n"); exit(1); }
    //Cluster Vertex is out of bounds
    if(vi < 0 || vi >= G->n_nodes){
        perror("Error in findPolarCoordinates in CreateClusters\n");
        exit(1);
    }
    // Assign Cluster Vertex as origin point instead of 0,0 in the cartesian plane 
    int x0, y0;
    x0 = G->a_combined[vi].x;
    y0 = G->a_combined[vi].y;

    int dx, dy;
    for(int i = 0; i < G->n_customers; i++){
        //If Cluster Vertex
        if(i == vi){
            continue;
        //If not Cluster Vertex then find polar coordinates based on Cluster Vertex position
        }else{
            dx = G->a_combined[i].x - x0;
            dy = G->a_combined[i].y - y0;
            polarArray[i].distance = sqrt((dx*dx) + (dy*dy));
            polarArray[i].degrees = calculate_angle(dx, dy);
            
            if(polarArray[i].degrees < 0.0)     //Transform negative angle
                polarArray[i].degrees += 360.0;
        }
    }
    return;
}

void assign_vertices_to_sectors(SON *G, int vi, polar *polarArray, sector *sec, int n_sect){
    if(!polarArray){ perror("polarArray == NULL at assign_vertices_to_sectors\n"); exit(1); }
    if(!sec){ perror("sec == NULL at assign_vertices_to_sectors\n"); exit(1); }
    //Initialize sectors
    for(int k = 0; k < n_sect; k++)
        sec[k].vertexList = NULL;
    
    //Assign nodes to sectors
    for(int i = 0; i < G->n_customers; i++){
        int sc = 0; //sector index [0, n_sect-1]
        if(i == vi) //Cluster Vertex is not in any sector
            continue;
        
        for(double j = 360.0/n_sect; j <= 360.0; j += 360.0/n_sect){
            if(polarArray[i].degrees <= j){
                append(&sec[sc].vertexList, G->a_customers[i].id);   // Add vertex id to its sector
                break;
            }
            sc++;
        }
    }
    return;
}

bool not_empty(int *array, int length){
    if(!array){ perror("array == NULL at not_empty()\n"); exit(1); }
    for(int i = 0; i < length; i++){
        if(array[i] == 1){
            return true;
        }
    }
    return false;
}

bool search_icluster(int *K_v, int icluster, int *v_free, int n_size){

    if(!K_v){ perror("K_v == NULL at search_cluster in misc_functions\n"); exit(1); }
    for(int i = 0; i < n_size; i++){ 
        if(K_v[icluster*n_size+i] > 0 && v_free[K_v[icluster*n_size+i]-1] == 1)
            return true;
    }
    return false;
}

int find_closest_free_vertex_in_sector(polar *polarArray, node *vertexList){
    if(!polarArray){ perror("polarArray == NULL at find_closest_free_vertex_in_sector\n"); exit(1); }
    if(!vertexList) { perror("vertexList == NULL at find_closest_free_vertex_in_sector()"); exit(1); }   //if sector j is empty return -1
    node *temp = vertexList; 
    double minDistance = HUGE_VAL;
    int minID = -1;
    //Search sector for closest vertex (At this point all vertices are still free anyway)
    while(temp){
        if(polarArray[temp->data - 1].distance < minDistance - epsilon){   //Compare vertices in sector
            minDistance = polarArray[temp->data - 1].distance;
            minID = temp->data;
        }
        temp = temp->next;
    }     
    return minID;    //return ID of closest vertex in sector j
}

int find_closest_free_vertex(int *v_free, polar *polarArray, int n_customers){
    if(!v_free){ perror("v_free == NULL at find_closest_free_vertex\n"); exit(1); }
    if(!polarArray){ perror("polarArray == NULL at find_closest_free_vertex\n"); exit(1); }
    double minDistance = HUGE_VAL;
    int minID = -1;
    for(int i = 0; i < n_customers; i++){
        if(v_free[i] == 1 && polarArray[i].distance < minDistance - epsilon){
            minDistance = polarArray[i].distance;
            minID = i+1;
        }
    }
    if(minID < 1){
        printf("minID < 1\n");
        exit(1);
    }
    return minID;       //return ID of closest unvisited vertex (no matter the sector)
}

void find_free_in_icluster(node **head, int *v_free, int *K_v, int icluster, int n_size){
    if(!v_free){ perror("v_free == NULL at find_free_in_cluster_k\n"); exit(1); }
    if(!K_v){ perror("K_pos == NULL at find_free_in_cluster_k\n"); exit(1); }
    for(int j = 0; j < n_size; j++){
        //If vertex is customer ID and is free then push into v_cand
        if(K_v[(icluster * n_size) + j] > 0){
            if(v_free[K_v[(icluster * n_size) + j] - 1] == 1)
                push(head, K_v[(icluster * n_size) + j]);
        }
    }
   return; 
}

int errCheck(int value, int lower, int upper){
    if(value < lower || value > upper)
        return 1;
    return 0;
}

void new_best(vt_solution *new_r, vt_solution *r, SON *G){
    for(int idep = 0; idep < G->n_depots; idep++){
        new_r->a_depots[idep].distance = r->a_depots[idep].distance;
        new_r->a_depots[idep].quantity_served = r->a_depots[idep].quantity_served;
        deleteList(&new_r->a_depots[idep].routelist);                                         //Delete old 
        new_r->a_depots[idep].routelist = copyList(r->a_depots[idep].routelist);
    }
    new_r->total_distance = r->total_distance;
    return;
}

int store_edge_count(vt_solution *Ra, SON *G, int *edge_matrix){

    if(!edge_matrix){ perror("p == NULL at store_edge_count\n"); exit(1); }
    int edge_sum = 0;

    for(int idep = 0; idep < G->n_depots; idep++){
        node *temp = Ra->a_depots[idep].routelist;
        if(!temp)
        {
            perror("temp == NULL at function store_edge_count...\n"); 
            exit(1);
        }
        while(temp->next){
            if(temp->data < 1 || temp->data > G->n_nodes || temp->next->data < 1 || temp->next->data > G->n_nodes){
                perror("Wrong values at temp->data at store_edge_count\n");
                exit(1);
            }
            edge_matrix[(G->n_nodes)*(temp->data - 1) + (temp->next->data - 1)] += 1;
            edge_sum++;

            temp = temp->next;
        }
    }
    return edge_sum;
}


void update_pheromones(SON *G, double *phMatrix, vt_solution *R, vt_solution *R_best, double T_update, double d){
    
    if(!phMatrix){ perror("phMatrix == NULL at update_pheromone\n"); exit(1); }
    
    double probability_R_best, t, x, R_update; 
    bool flag = false;

    if(R->total_distance < R_best->total_distance - epsilon){           
        t = (R_best->total_distance - R->total_distance)/R->total_distance;
        probability_R_best = exp(-t/T_update);                      
        printf("probability_R_best = %0.3lf\n", probability_R_best);
        if((double)random()/(double)RAND_MAX < probability_R_best - epsilon){
            R_update = R_best->total_distance;                          //Use worse solution R_best
            flag = true;
        }else{
            R_update = R->total_distance;
        }
    }else{                                                          
        R_update = R->total_distance;                                   //In this case R == R_best
    }

    node *temp = NULL;
    for(int idep = 0; idep < G->n_depots; idep++){
        if(flag){
            temp = R_best->a_depots[idep].routelist;
        }else{
            temp = R->a_depots[idep].routelist;
        }
        if(!temp){ perror("temp = NULL at update_pheromone \n"); printf("idep = %d\n", idep); continue; }
        while(temp && temp->next){
            if(phMatrix[idep*G->n_nodes*G->n_nodes + (temp->data-1)*G->n_nodes + (temp->next->data - 1)] < 0.0){
                perror("pheromone update error\n"); 
                printf("temp->data %d temp->next->data %d\n", temp->data, temp->next->data);
                printf("\nphMatrix = %0.2lf\n", 
                       phMatrix[idep*G->n_nodes*G->n_nodes + (temp->data-1)*G->n_nodes + (temp->next->data - 1)]);
                exit(1);
            }

            phMatrix[idep*G->n_nodes*G->n_nodes + (temp->data-1)*G->n_nodes + (temp->next->data - 1)] += 
                    d*(R->total_distance/R_update);

            phMatrix[idep*G->n_nodes*G->n_nodes + (temp->next->data-1)*G->n_nodes + (temp->data - 1)] += 
                    d*(R->total_distance/R_update);

            temp = temp->next;
        }
    }

    return;
}

double evaporate_pheromones(SON *G, int *edge_matrix, double *phMatrix, int edge_sum, int n_ants, double p_min, double p_max){
    
    if(!phMatrix){ perror("phMatrix == NULL at evaporate_pheromones\n"); exit(1); }
    if(!edge_matrix){ perror("edge_matrix == NULL at evaporate_pheromones\n"); exit(1); }

    double H = 0.0, *p;
    if(NULL == (p = calloc(G->n_nodes*G->n_nodes, sizeof *p))){
        printf("Error callocing p[]\n");
        exit(1);
    }
    
    //Upper triangle of matrix used
    for(int i = 1; i < G->n_nodes; i++){
        for(int j = 0; j < i; j++){
            if(i != j){
                p[(G->n_nodes)*j + i] = (double) edge_matrix[(G->n_nodes)*j + i] / edge_sum +
                                        (double) edge_matrix[(G->n_nodes)*i + j] / edge_sum;
            }
        }
    }

    for(int i = 1; i < G->n_nodes; i++){
        for(int j = 0; j < i; j++){
            if(i != j && p[(G->n_nodes*j + i)] > 0){
                H -= p[(G->n_nodes)*j + i] * log2(p[(G->n_nodes)*j + i]);
            }
        }
    }
    free(p);

    double H_min, H_max;
    H_min = -log2((double)n_ants/edge_sum);
    H_max = -log2(1.0/edge_sum);

    double rho = p_min + (p_max - p_min)*((H - H_min)/(H_max - H_min));

    for(int idep = 0; idep < G->n_depots; idep++){
        for(int c1 = 1; c1 < G->n_nodes; c1++){
            for(int c2 = 0; c2 < c1; c2++){
                if(phMatrix[idep*G->n_nodes*G->n_nodes + G->n_nodes*c1 + c2] != -1.0 &&
                    phMatrix[idep*G->n_nodes*G->n_nodes + G->n_nodes*c1 + c2] > 0.0001
                ){
                    phMatrix[idep*G->n_nodes*G->n_nodes + G->n_nodes*c1 + c2] *= (1.0 - rho);
                }

                if( phMatrix[idep*G->n_nodes*G->n_nodes + G->n_nodes*c2 + c1] != -1.0 &&
                    phMatrix[idep*G->n_nodes*G->n_nodes + G->n_nodes*c2 + c1] > 0.0001
                ){
                    phMatrix[idep*G->n_nodes*G->n_nodes + G->n_nodes*c2 + c1] *= (1.0 - rho);
                }
            }
        }
    }

    printf("rho = %0.3lf\n", rho);
    return (H-H_min)/H_min; //Termination condition
}

static int loop_iteration = 1;

void fprint_results(vt_solution *R, SON *G){

    int detect_dup[G->n_customers];
    for(int i = 0; i < G->n_customers; i++)
        detect_dup[i] = 0;

    char file_name[18], fn[18];
    sprintf(file_name, "p%02d_%d.res", instance_id, loop_iteration);

    FILE *fp, *fp_1;
    if(NULL == (fp = fopen(file_name, "w")))
    {
        perror("Couldn't open file fp at fprint_results\n");
        exit(1);
    }
    fprintf(fp, "total distance: %0.2lf\n", R->total_distance);

    node *temp = NULL, *vehicleRoute = NULL;
    push(&vehicleRoute, 0);
    sprintf(fn, "p%02d_%d.csv", instance_id, loop_iteration);
    if(NULL == (fp_1 = fopen(fn, "w"))){
        perror("Error opening fp_1!\n");
        exit(1);
    }

    double distance = 0.0;
    fprintf(fp_1, "StartNode,EndNode\n");
    for(int idep = 0; idep < G->n_depots; idep++){
        int vehicle = 1;
        int q_served = 0;
        if(!R->a_depots[idep].routelist){
            continue;
        } 
        temp = R->a_depots[idep].routelist;
        while(temp && temp->next){
            distance += G->d_matrix[temp->data - 1][temp->next->data - 1];
            if(temp->next->data > G->n_customers){
                append(&vehicleRoute, 0);
                fprintf(fp, "depot %d  route %2d  q_served %3d\t", idep+G->n_customers+1, vehicle, q_served);
                fprintf(fp_1, "%d,%d\n", temp->data, temp->next->data);
                fprintList(vehicleRoute, fp);
                deleteList(&vehicleRoute);
                vehicle++;	
                q_served = 0;
                push(&vehicleRoute, 0);
            }else{
                q_served += G->a_customers[temp->next->data - 1].demand;
                append(&vehicleRoute, temp->next->data);
                detect_dup[temp->next->data - 1]++;
                fprintf(fp_1, "%d,%d\n", temp->data, temp->next->data);
            }
            temp = temp->next;
        }
    }

    fprintf(fp, "total_distance = %0.3lf\n", distance);

    for(int i = 0; i < G->n_customers; i++){
        if(detect_dup[i] != 1){
            fprintf(fp, "detect_dup[%d] = %d\n", i, detect_dup[i]);
        }
    }

    if(fclose(fp_1) != 0){
        perror("Error closing fp_1!\n");
        exit(1);
    }

    if(fclose(fp) != 0){
        perror("Error closing fp!\n");
        exit(1);
    }

    free(vehicleRoute);

    return;
}

void fprint_data(int iterations, int best_iter, double foundtime, double runtime){
    char file_name[18];
    sprintf(file_name, "p%02d_%d.data", instance_id, loop_iteration);
    FILE *fp;
    if(NULL == (fp = fopen(file_name, "w"))){
        perror("Error opening data file!\n");
        exit(1);
    }

    fprintf(fp, "Solution found in: %.2lf seconds in iteration %d\n\nTotal runtime: %.2lf seconds\nTotal iterations: %d",foundtime, best_iter, runtime, iterations);
    
    if(fclose(fp) != 0){
        perror("Error closing data file!\n");
        exit(1);
    }
    
    loop_iteration++;

    return;
}
