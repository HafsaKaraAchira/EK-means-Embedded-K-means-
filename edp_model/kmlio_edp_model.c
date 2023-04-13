# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

// (clear && gcc -g edp_model/kmlio_edp_model.c -o edp_model/kmlio_edp_model -lm -D _GNU_SOURCE) && (edp_model/kmlio_edp_model 0.2 13421800 6710900 10 10)

// get mat on NI * CPI = 1491.543

// get mat cpi 40  * 0.67 * M
// get mat conv  1138.6 * 0.442 * M *dim
//1248.932 * 0.4418

// 1492.5

double getmatrix_1copy_IC = 40 ;
double getmatrix_1copy_CPI = 0.67 ;

double getmatrix_1conversion_IC = 1248.932 ;//1138.6;
double getmatrix_1conversion_CPI = 0.4418 ;//0.442 ;

double affect_x_nb_instructions = 1959 ;
double affect_x_avg_CPI = 0.326 ;

double find_best_center_nb_instructions = 13098.13647 ;
double find_best_center_avg_CPI = 0.334 ;

double km_update_nb_instructions = 76.7545 ;    // 36.656 ;  // 184.68 ;
double km_update_avg_CPI =  0.3257 ;    //0.3195 ;//0.35 ;

double km_dist_mat_nb_instructions = 7.9955 ; //18.4 ;
double km_dist_mat_avg_CPI = 0.3209 ;

double km_reassign_nb_instructions = 5.5299 ;   //184.68 ;
double km_reassign_avg_CPI = 0.5119 ; //0.45...   // 0.35 ;

double km_copy_nb_instructions = 28.9054 ; //28.9200 ; 
double km_copy_avg_CPI = 0.4640; //0.5119 ; 

double km_IC = 1678.9800 ;  //1666.5542 ;
double km_CPI = 0.3402 ;    //0.3268 ;

double form_index_nb_instructions = 150 ;
double form_index_avg_CPI = 0.326 ;

// //////////////////////////////////////////////////////////////

double C_eff = 0.8 ;
double alpha = 0.5 ;

int IT_MAX = 150 ;
double T_1read = 14 ; //12.23 ; //octet read latency in nano seconds

//decide best buffer size
int MAX_DOUBLE_LEN = 17 ;

int MAX_FREQ = 2900000 ;
int MIN_FREQ = 800000 ;
int FREQ_STEP = 100000 ;

int * available_frequencies=NULL ;

// /////////////////////////////////////////////////////////////

double T_off_mark( int N,int dim){
    double D;
    // the average size of chunk line in file
    unsigned int avg_line_size = (MAX_DOUBLE_LEN+1) * dim ;

    D = N * avg_line_size * T_1read  ;
    D *= pow(10,-9) ;
    printf("\n-----------------------------------------------------------\n") ;
    printf("*** mark estimated time  = %f s\n",D) ;
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

    D = (*L_opt) * (*it_num) * avg_line_size * T_1read * pow(10,-9) ;
    
    return D ;
}

double T_on_getmatrixbuf(int M,int dim,int k,int L_opt,int it_num,long frequency){
    double D;
    D = it_num * L_opt * ( dim * getmatrix_1conversion_IC * getmatrix_1conversion_CPI + getmatrix_1copy_IC * getmatrix_1copy_CPI ) / frequency ;
    return D ;
}

double T_getmatrixbuf(int M,int dim,int k,long frequency){
    double D;
    int L_opt,it_num ;
    double T_on , T_off ;
    T_off = T_off_getmatrixbuf(M,dim,k,&L_opt,&it_num) ;
    T_on = T_on_getmatrixbuf(M,dim,k,L_opt,it_num,frequency) ; 

    D = T_off + T_on ;
    printf("\n-----------------------------------------------------------\n") ;
    printf("*** getmatrix estimated time :\n\tL_opt = %d ,\tit_num = %d\n\toff = %f s\n\ton = %f s\n\ttotal = %f s\n",L_opt,it_num,T_off,T_on,D) ;

    return D ;
}

double T_on_init_plus_plus(int M,int dim,int k, long frequency ){
    double T_on_affect_x = ( affect_x_nb_instructions * affect_x_avg_CPI ) / frequency ;
    double T_on_find_best_center = ( find_best_center_nb_instructions * find_best_center_avg_CPI ) / frequency ;
    double D;
    D = k * ((M * T_on_affect_x) + T_on_find_best_center) ;
    printf("\n-----------------------------------------------------------\n") ;
    printf("*** init++ estimated time  = %f s\n",D) ;
    return D ;
}

double T_on_km_update( int M, int dim, int k, long frequency ){
    double D ;
    D = M * dim * km_update_nb_instructions * km_update_avg_CPI / frequency ;
    return D ;
}

double T_on_km_dist_mat( int M, int dim, int k, long frequency ){
    double D ;
    D = M * k * dim * km_dist_mat_nb_instructions * km_dist_mat_avg_CPI / frequency ;
    return D ;
}

double T_on_km_reassign( int M, int dim, int k, long frequency ){
    double D ;
    D = M * k * km_reassign_nb_instructions * km_reassign_avg_CPI / frequency ;
    return D ;
}

double T_on_km_copy( int M, int dim, int k, long frequency ){
    double D ;
    D = M * km_copy_nb_instructions * km_copy_avg_CPI / frequency ;
    return D ;
}

double T_on_km_setup( int M, int dim, int k, long frequency ){
    double D ;
    D = T_on_km_dist_mat(M,dim,k,frequency) + T_on_km_reassign(M,dim,k,frequency) ;
    return D ;
}

double T_on_km_1iteration( int M, int dim, int k, long frequency ){
    double D ;
    
    D = T_on_km_update(M,dim,k,frequency) + T_on_km_dist_mat(M,dim,k,frequency) + T_on_km_reassign(M,dim,k,frequency) + T_on_km_copy(M,dim,k,frequency) ; 
    printf("\n-----------------------------------------------------------\n") ;
    printf("*** kmeans 1 iteration estimated time :\n") ;
    printf("\tcenters update = %f s\n",T_on_km_update(M,dim,k,frequency)) ;
    printf("\tdistance matrix calculus = %f s\n",T_on_km_dist_mat(M,dim,k,frequency)) ;
    printf("\tpoints reassign time = %f s\n",T_on_km_reassign(M,dim,k,frequency)) ;
    printf("\n\tpoints copy time = %f s\n",T_on_km_copy(M,dim,k,frequency)) ;
    printf("\ttotal time in details = %f s\n",D) ;
    
    D = M * km_IC *km_CPI /frequency ;
    
    printf("\ttotal time in average = %f s\n",D) ;
    return D ;
}

int km_nb_iteration_estimation(int M, int dim, int k, int sepVal ){
    return IT_MAX ;
}

double T_on_km_iterations( int M, int dim, int k, long frequency,double sepVal){
    double D;
    D = T_on_km_setup(M,dim,k,frequency) + ( km_nb_iteration_estimation(M,dim,k,sepVal) * T_on_km_1iteration(M,dim,k,frequency)) ;
    return D ;
}

double T_on_form_index(int N , int M,int dim,int k,long frequency){
    double D ;
    D = k * ((int)N/M) * N * form_index_nb_instructions * form_index_avg_CPI / frequency ;
    return D ;
}

double T_form_final_chunk(int N,int M,int dim,int k,long frequency){
    double D;
    int L_opt,it_num ;
    double T_form_index , T_dataset ;

    T_dataset = ((int)N/M) * T_getmatrixbuf(M,dim,k,frequency) ;
    T_form_index = T_on_form_index(N,M,dim,k,frequency) ; 
    D = T_form_index + T_dataset ;
    
    printf("\n-----------------------------------------------------------\n") ;
    printf("*** form final chunk estimated time :\n") ;
    printf("\tdataset read = %f s\n",T_dataset) ;
    printf("\tgroups index forming = %f s\n",T_form_index) ;
    printf("\ttotal = %f s\n\n",D) ;
    
    return D ;
}

double T_cp_max(int sepVal,int N,int M,int dim,int k,long frequency){
    double D ;
    D = T_getmatrixbuf(M,dim,k,frequency) + 
        T_on_init_plus_plus(M,dim,k,frequency) + 
        T_on_km_iterations(M,dim,k,frequency,sepVal)  ;

    printf("\n-----------------------------------------------------------\n") ;
    printf("*** Partiel chunk calculus estimated time  = %f s\n",D) ;

    return D ;
}

double T_cf_max(int sepVal,int N,int M,int dim,int k,long frequency){
    double D ;
    D = T_form_final_chunk(N,M,dim,k,frequency) + 
        T_on_init_plus_plus(M,dim,k,frequency) + 
        T_on_km_iterations(M,dim,k,frequency,sepVal)  ;

    printf("\n-----------------------------------------------------------\n") ;
    printf("*** Final chunk calculus estimated time  = %f s\n\n",D) ;

    return D ;
}

double estimate_Kmeans_time(int sepVal,int N,int M,int dim,int k,long frequency){
    double T_total = T_off_mark(N,dim) + 
                    ( 
                        T_getmatrixbuf(M,dim,k,frequency) + 
                        T_on_init_plus_plus(M,dim,k,frequency) + 
                        T_on_km_iterations(M,dim,k,frequency,sepVal) 
                    ) ;

    return T_total ;
}

double estimate_KMLIO_time(int sepVal,int N,int M,int dim,int k,long frequency){
    double T_total = T_off_mark(N,dim) + 
                    ((int)(N/M)) *  T_cp_max(sepVal,N,M,dim,k,frequency) +
                    T_cf_max(sepVal,N,M,dim,k,frequency) ;


    return T_total ;
}

int * get_available_frequencies(){
    int * available_frequencies = malloc(sizeof(int)*50) ;
}

int main (int argc, char **argv){
    size_t N,M,dim,k;
    double D_max ;
    double sepVal ;

    sepVal = atof(argv[1]) ;
	N = atoi(argv[2]); 
	M = atoi(argv[3]);
	dim = atoi(argv[4]); 
    k = atoi(argv[5]);
    // D_max = atoi(argv[6]);  // maximal K-MLIO execution delay constraint in seconds

    long frequency = 2250 * pow(10,6) ;

    printf("Execution estimation for N=%lu , M=%lu , dim = %lu , k=%lu\n",N,M,dim,k) ;
    // T_getmatrixbuf(M,dim,k,frequency);
    // T_on_km_1iteration(M,dim,k,frequency) ;
    // T_form_final_chunk(N,M,dim,k,frequency) ;
    
    printf("\n==========================================================================\n") ;

    printf("KMLIO total execution estimated time = %f s\n\n",estimate_KMLIO_time(sepVal,N,M,dim,k,frequency)) ;






    // ///////////////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////////////



    // frequency = MAX_FREQ ;

    // int abd_min  = -1 ;

    // int M_max = M ;
    // double T_est ;

    // while (abd_min == -1 )  //search continusly for a minimal skipped chunk number 
    // {
    //     // while the skip rate is lower than the total number of chunks  = partiels + final
    //     for(int abd = 0 ; abd <((int)(N/M)+1) ; abd++){
    //         // search if there is a maximal frequency 
    //         // that realize an execution within the maximal delay constraint
    //         for(long freq=MAX_FREQ; freq<=MIN_FREQ ; freq-=FREQ_STEP){
    //             // if skip rate is equal to partiel chunks number - 1
    //             // then we estimate for classic kmeans otherwise we stimate time for K-MLIO
    //             if(abd == ((int)(N/M)-1) ){
    //                 T_est = estimate_Kmeans_time(sepVal,N,M,dim,k,freq) ;
    //             }else{
    //                 T_est = estimate_KMLIO_time(sepVal,N,M,dim,k,freq) ;
    //             }
    //             // if estimated time is less than maximal delay constraint
    //             // then we fix minimal skipped chunk number
    //             if(T_est < D_max){
    //                 abd_min = abd ;
    //                 break ;
    //             }
    //         }
    //         // if we can't even execute K-mlio nor Kmeans 
    //         // on a signle chunk within D-max
    //         // the we choose a new smaller chunk size M
    //         if (abd == ((int)(N/M)))
    //         {
    //             for ( int m = ((int)(N/M)) ; m<N ; m++){
    //                 if (N%m == 0 )
    //                 {
    //                     M = N/m ;
    //                     break ;
    //                 }
                    
    //             }
    //         }
    //     }
    // }





    // ///////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////

    
    return 0 ;
}


///////////////////////

// E = C_eff * pow(CST,2) * pow(frequency,3) * D ;
// return E*D ;