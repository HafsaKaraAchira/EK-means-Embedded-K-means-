# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>


double getmatrix_1conversion_nb_instructions = 14496 ;
double getmatrix_1conversion_avg_CPI = 0.44 ;

double affect_x_nb_instructions = 1959 ;
double affect_x_avg_CPI = 0.326 ;

double find_best_center_nb_instructions = 13098.13647 ;
double find_best_center_avg_CPI = 0.334 ;

double km_update_nb_instructions = 10 ;
double km_update_avg_CPI = 10 ;

double km_dist_mat_nb_instructions = 10 ;
double km_dist_mat_avg_CPI = 10 ;

double km_reassign_nb_instructions = 10 ;
double km_reassign_avg_CPI = 10 ;

double km_setup_nb_instructions = 10 ;
double km_setup_avg_CPI = 10 ;

double form_index_nb_instructions = 10 ;
double form_index_avg_CPI = 10 ;

double C_eff = 0.8 ;
double alpha = 0.5 ;
int IT_MAX = 150 ;
double T_1read = 15.7E-9 ;

//decide best buffer size
int MAX_DOUBLE_LEN = 17 ;

int MAX_FREQ = 2900000 ;
int MIN_FREQ = 400000 ;

double T_off_mark( int N,int dim){
    double D;
    // the average size of chunk line in file
	int avg_line_size = (MAX_DOUBLE_LEN+1) * dim ;

    D = N * avg_line_size * T_1read ;
    return D ;
}

double T_off_getmatrixbuf(int M,int dim,int k,int *L_opt,int *it_num){
    double D;
    // calculate the rest size allowed for the chunk km
	int km_size = (sizeof(double)*k + 3*sizeof(int)) * M ;
	
	// the average size of chunk line in file
	int avg_line_size = (MAX_DOUBLE_LEN+1) * dim ;
	// the maximum number of lines to store in buffer
	int L_max = km_size/(avg_line_size+1+sizeof(char *)) ;
	
    // ajust the number of lines in buffer 
	*it_num = ceil((double)M/L_max) ;
	
    // balance the number of read lines in each iteration
	*L_opt = ceil((double)M / (*it_num)) ;

    D = (*L_opt) * (*it_num) * avg_line_size * T_1read ;
    
    return D ;
}

double T_on_getmatrixbuf(int M,int dim,int k,int L_opt,int it_num,double frequency){
    double D;
    D = it_num * L_opt * getmatrix_1conversion_nb_instructions * getmatrix_1conversion_avg_CPI / frequency ;
    return D ;
}

double T_getmatrixbuf(int M,int dim,int k,double frequency){
    double D;
    int L_opt,it_num ;
    double T_on , T_off ;
    T_off = T_off_getmatrixbuf(M,dim,k,&L_opt,&it_num) ;
    T_on = T_on_getmatrixbuf(M,dim,k,L_opt,it_num,frequency) ; 

    D = T_off + T_on ;

    printf("getmatrix estimated time :\n\t \n\tL_opt = %d\n\tit_num = %d\n\toff = %f s\n\ton = %f s\n\ttotal = %f s\n",L_opt,it_num,T_off,T_on,D) ;

    return D ;
}

double T_on_init_plus_plus(int M,int dim,int k, double frequency ){
    double T_on_affect_x = ( affect_x_nb_instructions * affect_x_avg_CPI ) / frequency ;
    double T_on_find_best_center = ( find_best_center_nb_instructions * find_best_center_avg_CPI ) / frequency ;
    double D;
    D = k * ((M * T_on_affect_x) + T_on_find_best_center) ;
    return D ;
}

double T_on_km_update( int M, int dim, int k, double frequency ){
    double D ;
    D = M * dim * km_update_nb_instructions * km_update_avg_CPI / frequency ;
    return D ;
}

double T_on_km_dist_mat( int M, int dim, int k, double frequency ){
    double D ;
    D = M * k * dim * km_dist_mat_nb_instructions * km_dist_mat_avg_CPI / frequency ;
    return D ;
}

double T_on_km_reassign( int M, int dim, int k, double frequency ){
    double D ;
    D = M * k * km_reassign_nb_instructions * km_reassign_avg_CPI / frequency ;
    return D ;
}

double T_on_km_setup( int M, int dim, int k, double frequency ){
    double D ;
    D = T_on_km_dist_mat(M,dim,k,frequency) + T_on_km_reassign(M,dim,k,frequency) ;
    return D ;
}

double T_on_km_1iteration( int M, int dim, int k, double frequency ){
    double D ;
    D = T_on_km_update(M,dim,k,frequency) + T_on_km_dist_mat(M,dim,k,frequency) + T_on_km_reassign(M,dim,k,frequency) ;
    return D ;
}

int km_nb_iteration_estimation(int M, int dim, int k, int sepVal ){
    return IT_MAX ;
}

double T_on_km_iterations( int M, int dim, int k, double frequency,double sepVal){
    double D;
    D = T_on_km_setup(M,dim,k,frequency) + ( km_nb_iteration_estimation(M,dim,k,sepVal) * T_on_km_1iteration(M,dim,k,frequency)) ;
    return D ;
}

double T_on_form_index(int N , int M,int dim,int k,double frequency){
    double D;
    D = k * ((int)N/M) * form_index_nb_instructions * form_index_avg_CPI / frequency ;
    return D ;
}

double T_form_final_chunk(int N,int M,int dim,int k,double frequency){
    double D;
    int L_opt,it_num ;
    double T_form_index , T_dataset ;
    T_dataset = ((int)N/M) * T_getmatrixbuf(M,dim,k,frequency) ;
    T_form_index = T_on_form_index(N,M,dim,k,frequency) ; 
    D = T_form_index + T_dataset ;
    return D ;
}

int main (int argc, char **argv){
    size_t N,M,dim,k;
    double sepVal ;

    sepVal = atof(argv[1]) ;
	N = atoi(argv[2]); 
	M = atoi(argv[3]);
	dim = atoi(argv[4]); 
    k = atoi(argv[5]);

    long frequency = 2300 * pow(10,6) ;
    
    printf("mark estimated time  = %f s\n",T_off_mark(N,dim)) ;
    T_getmatrixbuf(M,dim,k,frequency);
    printf("init++ estimated time  = %f s\n",T_on_init_plus_plus(M,dim,k,frequency)) ;
    // printf("getmatrix estimated time :\n\t off = %f s\n\t on = %f s\n\t total = %f s",T_getmatrixbuf(M,dim,k,frequency)) ;
    // T_on_km_iterations(M,dim,k,frequency,sepVal) ;
    // T_on_form_index(N,M,dim,k,frequency) ;

    // double T_total = T_off_mark(N,dim) + 
    //                 ((int)(N/M)+1) * 
    //                     ( 
    //                         T_getmatrixbuf(M,dim,k,frequency) + 
    //                         T_on_init_plus_plus(M,dim,k,frequency) + 
    //                         T_on_km_iterations(M,dim,k,frequency,sepVal) 
    //                     ) +
    //                         T_on_form_index(N,M,dim,k,frequency) ;

	
    return 0 ;
}


///////////////////////

// E = C_eff * pow(CST,2) * pow(frequency,3) * D ;
// return E*D ;