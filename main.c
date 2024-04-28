#include "header_files/AACO_main_functions.h"
#include "header_files/structs.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int instance_id;
int main(int argc, char **argv) {

    if(argc != 2){
        printf("Usage: <program> <instance_id>\n");
        exit(EXIT_SUCCESS);
    }else{
        if( sscanf(argv[1], "%02d", &instance_id) != 1) {
            printf("Invalid instance ID\n");
            exit(EXIT_FAILURE);
        }
    }

    /* Read file and store all info */
    FILE *fp;
    if(NULL == (fp = fopen("../../Cordeau_dataset/p01.txt", "r"))){
        perror("Instance file error!\n");
        exit(1);
    }

    SON G;

    /* Get first line and display data info */
    char line[100];
    fgets(line, sizeof line, fp);
    sscanf(line, "%*d %*d %d %d", &G.n_customers, &G.n_depots);
    printf("numOfCustomers = %d, numOfDepots = %d\n", G.n_customers, G.n_depots);

    G.n_nodes = G.n_customers + G.n_depots;
    /* Initialize arrays of SetOfNodes G */
    if (NULL == (G.a_depots = malloc(sizeof *G.a_depots * G.n_depots)) ||
        NULL == (G.a_customers = malloc(sizeof *G.a_customers * G.n_customers)) ||
        NULL == (G.a_combined = malloc(sizeof *G.a_combined * G.n_nodes))) {
        perror("Error mallocing G arrays!\n");
        exit(1);
    }

    for (int i = 0; i < G.n_depots; i++) {
        fgets(line, sizeof line, fp);
        sscanf(line, "%*d %d", &G.capacity);
    }

    for(int i = 0; i < G.n_customers; i++){
        fgets(line, sizeof line, fp);
        sscanf(line, "%d %d %d %*d %d", &G.a_customers[i].id, &G.a_customers[i].x, &G.a_customers[i].y, &G.a_customers[i].demand);

        G.a_combined[i].id = G.a_customers[i].id;
        G.a_combined[i].x = G.a_customers[i].x;
        G.a_combined[i].y = G.a_customers[i].y;
    }

    /* Get depots info */
    for(int i = 0; i < G.n_depots; i++){
        fgets(line, sizeof line, fp);
        sscanf(line, "%d %d %d %*d", &G.a_depots[i].id, &G.a_depots[i].x, &G.a_depots[i].y);

        G.a_combined[G.n_customers + i].id = G.a_depots[i].id;
        G.a_combined[G.n_customers + i].x = G.a_depots[i].x;
        G.a_combined[G.n_customers + i].y = G.a_depots[i].y;
    }

    if (fclose(fp) != 0) {
        perror("fclose error at main\n");
        exit(1);
    }
    if(NULL == (G.d_matrix = malloc(sizeof *G.d_matrix * G.n_nodes))){
        printf("Error mallocing d_matrix!\n");
        exit(1);
    }
    for(int i = 0; i < G.n_nodes; i++){
        if(NULL == (G.d_matrix[i] = malloc(sizeof *G.d_matrix[i] * G.n_nodes))){
            printf("Error mallocing d_matrix[]!\n");
            exit(1);
        }
    }

    for(int i = 0; i < G.n_nodes; i++){
        for(int j = 0; j < G.n_nodes; j++){
            if(i == j || (i >= G.n_customers && j >= G.n_customers)){
                G.d_matrix[i][j] = 0.0;
            }else{
                G.d_matrix[i][j] = sqrt(pow(G.a_combined[i].x - G.a_combined[j].x, 2) + pow(G.a_combined[i].y - G.a_combined[j].y, 2));
            }

        }
    }

    // n_prim * n_size must not exceed total number of customers
    int n_ants = 192;
    int n_freq = 10;
    int n_size = 24;
    int n_sect = 16;
    int n_prim = 4;
    //int n_prim = (int)ceil((double)G.n_customers/n_size);
    if((int)ceil((double)G.n_customers/n_size) < n_prim)
        n_prim = (int)ceil((double)G.n_customers/n_size);
    double T_update = 0.1; // The higher it is - the higher probability to update
    // pheromones using the worse solution
    double a_update = 0.0001;
    double p_min = 0.001;
    double p_max = 0.1;
    double d = 3.0;
    double a = 1.0;
    double b = 1.0;


    srand(time(NULL));

    for (int i = 0; i < 10; i++){
        AACONC(&G, n_ants, n_freq, n_size, n_sect, n_prim, T_update,
               a_update, p_min, p_max, d, a, b);
    }

    /* End of Program */
    for(int i = 0; i < G.n_nodes; i++){
        free(G.d_matrix[i]);
    }
    free(G.d_matrix);


    free(G.a_customers);
    free(G.a_depots);
    free(G.a_combined);

    return 0;
}
