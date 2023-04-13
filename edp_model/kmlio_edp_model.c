# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

double kmeans_one_iteration_nb_instructions = 10 ;
double kmeans_one_iteration_avg_CPI = 10 ; 
C_eff = 0.8 ;
CST = 0.5 ;
int IT_MAX = 150 ;


double kmeans_one_iteration_edp_estimation( int M, int dim, int k, double frequency ){
    double D , E ;
    D = kmeans_one_iteration_nb_instructions * kmeans_one_iteration_avg_CPI / frequency ;
    E = C_eff * pow(CST,2) * pow(frequency,3) * D ;
    return E*D ;
}

int kmeans_nb_iteration_estimation( int N , int M, int dim, int k, int sepVal ){
    return IT_MAX ;
}