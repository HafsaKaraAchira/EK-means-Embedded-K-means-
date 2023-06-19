// kmeans.c
// Ethan Brodsky
// modified by Camélia SLIMANI
// modified by Hafsa KARA ACHIRA

/*

CMD :

PARAMETERS : $K $N $DIM $MEM_SIZE $NB_REPEAT

scripts/prog_skip_chunk/prog_script_beta 10 3355450 10 75 5

*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <fcntl.h>
#include <malloc.h>
#include <sys/mman.h>

// #include <conio.h>
#include <sys/stat.h>

#include <limits.h>
#include <sys/param.h>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <sched.h>

#define sqr(x) ((x) * (x))

#define MAX_CLUSTERS 100000

#define MAX_ITERATIONS 150

#define BIG_double (INFINITY)

struct kmlio_time_estimation
{
	double mark_time;
	double getmat_time;
	double init_time;
	double form_chunk_time;
	double km_1_iteration_time;
	double var_copy_time;

	double C_1_it_elem;
	double C_init_elem;
	double C_get_matrix_on_elem;
	double C_var_copy;
	double T_1_read;

	int nb_group ;
};
typedef struct kmlio_time_estimation kmlio_time_estimation;

struct skip_chunk_solution
{
	double skip_chunk;
	double freq;
	double estimated_curr_chunk_delay;
};
typedef struct skip_chunk_solution skip_chunk_solution;

struct chunk_stats
{
	double chunk_estimated_delay;
	double chunk_real_delay;
	double chunk_rem_checkpoint;
	int km_nb_iterations;
	int skp_stat;
	long freq;
};
typedef struct chunk_stats chunk_stats;

int chunk_ind = 0;

int skip_chunk = 0;
int D_max;
double beta;
int nb_it_sample = 1;
int corrected = 0;

struct timeval km_it_start, kmlio_start, kmlio_end;

double *kmeans_iterations_durations;

chunk_stats *kmlio_chunks_stats;

long MAX_FREQ = 2900000000, MIN_FREQ = 798000000, BASE_FREQ = 2294000000, FREQ_STEP = 100000000;

long set_speed;

const int nb_available_frequencies = 16;

long available_frequencies[] = {798000000, 897000000, 997000000, 1097000000, 1197000000, 1297000000, 1396000000, 1496000000, 1596000000, 1696000000, 1795000000, 1895000000, 1995000000, 2095000000, 2194000000, 2294000000};

kmlio_time_estimation real_time;

void set_frequency(long freq)
{
	set_speed = freq;
	char *cmd;
	asprintf(&cmd, "echo %ld | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_setspeed", (set_speed / 1000));
	system(cmd);

	cmd = NULL;
}

void set_governor(char *gov)
{
	char *cmd;
	asprintf(&cmd, "echo %s | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor", gov);
	system(cmd);
}

double calc_delai_time(struct timeval tv_start, struct timeval tv_end)
{
	return ((double)(tv_end.tv_usec - tv_start.tv_usec) / 1000000 + (double)(tv_end.tv_sec - tv_start.tv_sec));
}

double calc_history_past_time(struct timeval tv_start)
{
	struct timeval tv_now;
	gettimeofday(&tv_now, NULL);
	return calc_delai_time(tv_start, tv_now);
}

double calc_kmlio_past_time()
{
	struct timeval tv_now;
	gettimeofday(&tv_now, NULL);
	return calc_delai_time(kmlio_start, tv_now);
}

double calc_remaining_Time()
{
	return (D_max - calc_kmlio_past_time());
}

double estimate_T_get_matrix_time(size_t M, long freq)
{
	return (real_time.T_1_read * M + real_time.C_get_matrix_on_elem / freq);
}

double estimate_T_init_time(long freq)
{
	return (real_time.C_init_elem / freq);
}

double estimate_T_1_iteration_time(long freq)
{
	return (real_time.C_1_it_elem / freq);
}

double estimate_T_var_copy_time(long freq)
{
	return (real_time.C_var_copy / freq);
}

double estimate_wcet_partial_chunk(size_t M, long freq)
{
	return ( estimate_T_get_matrix_time(M, freq) + estimate_T_init_time(freq) + ((MAX_ITERATIONS + 1) * estimate_T_1_iteration_time(freq)) + estimate_T_var_copy_time(freq) );
}

double estimate_wcet_final_chunk(size_t N, size_t M, int skip_chunk, long freq)
{
	return (estimate_T_get_matrix_time(M, freq) * ((N / M) - skip_chunk) + estimate_T_init_time(freq) + ((MAX_ITERATIONS + 1) * estimate_T_1_iteration_time(freq)));
}

double calc_cycles_per_elem(size_t N, size_t M)
{
	real_time.T_1_read = real_time.mark_time / N ;
	real_time.C_get_matrix_on_elem = (real_time.getmat_time - real_time.T_1_read * M) * set_speed ;
	real_time.C_init_elem = real_time.init_time * set_speed ;
	real_time.C_1_it_elem = real_time.km_1_iteration_time * set_speed ;
	real_time.C_var_copy = real_time.var_copy_time * set_speed ;
}


// a local function to return the optimal frequency before processing each chunk
long get_optimal_freq(size_t N, size_t M, double *estimated_curr_chunk_delay, double *rem_time)
{
	// we suppose the optimal frequency is the maximum one
	long freq = available_frequencies[nb_available_frequencies - 1];
	// we supoose that we can't find a better optimization within thetime constraint
	int (*found) = 0;

	int exec_chunk_nb = (N / M) - (chunk_ind == 0) - chunk_ind , exec_nb_iterations = 0, final_chunk_exist = 1, rem_chunks = 0 , curr_chunk_read_time = 0 ;
	double req_time, T_read;

	// initialize chunks opartions number of cycles
	double c_km = real_time.C_1_it_elem , c_init = real_time.C_init_elem , c_getmat = real_time.C_get_matrix_on_elem , c_var_cpy = real_time.C_var_copy ;

	if (exec_chunk_nb == skip_chunk)	// if we reached final chunk
	{
		exec_nb_iterations += MAX_ITERATIONS + 1;
		c_getmat *= ((N / M) - skip_chunk); // fix

		// save chunk's estimated computation cycles and I/O time
		(*estimated_curr_chunk_delay) = (exec_nb_iterations * c_km + c_init + c_getmat);
		curr_chunk_read_time = ((N / M) - skip_chunk) * real_time.T_1_read * M;
	}
	else // execute partiel chunks other than the first + final chunk
	{
		exec_nb_iterations += MAX_ITERATIONS + 1;

		// save chunk's estimated computation cycles and I/O time
		(*estimated_curr_chunk_delay) = exec_nb_iterations * c_km + c_init + c_getmat;
		curr_chunk_read_time = real_time.T_1_read * M;
		
		// add other chunks values
		rem_chunks = exec_chunk_nb - skip_chunk + final_chunk_exist;
		exec_nb_iterations *= rem_chunks;
		c_init *= rem_chunks;
		c_var_cpy *= rem_chunks;
		c_getmat *= (exec_chunk_nb - skip_chunk + final_chunk_exist * ((N / M) - skip_chunk));
	}

	// calculate the dataset I/O time
	T_read = (exec_chunk_nb - skip_chunk + final_chunk_exist * ((N / M) - skip_chunk)) * real_time.T_1_read * M;
	// calculate the rrquired cycles to performe remaining chunks computations
	double req_cycles = c_getmat + c_init + c_var_cpy + c_km * exec_nb_iterations;
	// calculate remaining time
	(*rem_time) = calc_remaining_Time();
	// deduce the required frequincy constraint
	double req_freq = req_cycles / ((*rem_time) - T_read);
	// verify if the required frequency is possible
	(*found) = (req_freq > 0 && freq >= req_freq);

	if (*found) // IF A REQUIRED FREQ IS POSSIBLE , WE SELECT THE MINIMUM EXISTANT FREQ IN THE AVAILABLE FREQUENCIES VALUES
	{
		int min = 0;
		int max = nb_available_frequencies - 1;
		int optimal = (min + max) / 2;

		while (min <= max) // perform a simple dichotomic serach for the optimal frequency
		{
			if (available_frequencies[optimal] < req_freq)
				min = optimal + 1;
			else
				max = optimal - 1;

			optimal = (min + max) / 2;
		}
		freq = available_frequencies[min];
	}

	// save the estimated chunk time with choosen frequency
	(*estimated_curr_chunk_delay) = (*estimated_curr_chunk_delay) / freq + curr_chunk_read_time;

	//
	return freq;
}

void decide_skip_chunk(size_t N, size_t M)
{
	long freq = available_frequencies[nb_available_frequencies - 1];
	int skp = 0, found = 0 ;	// sol_ind = 0;
	double estimated_curr_chunk_delay;
	skip_chunk_solution skip_chk_sols[(N / M) - chunk_ind];
	double rem_time = calc_remaining_Time() ;

	while( ( rem_time = calc_remaining_Time() - estimate_wcet_final_chunk(N,M,skp,freq) ) <= 0 )
		skp ++ ;

	skip_chunk = MIN( N/M - 1 , rem_time / estimate_wcet_partial_chunk(M,freq) ) ;

	// while (skp < (N / M) - chunk_ind)
	// {
	// 	f = get_optimal_freq(N, M, skp, &found, &estimated_curr_chunk_delay, &rem_time);
	// 	if (found == 1)
	// 	{
	// 		skip_chunk = skp;
	// 		(*freq) = f;
	// 		kmlio_chunks_stats[chunk_ind].chunk_estimated_delay = estimated_curr_chunk_delay;
	// 		kmlio_chunks_stats[chunk_ind].chunk_rem_checkpoint = rem_time;
	// 		kmlio_chunks_stats[chunk_ind].freq = f;
	// 		kmlio_chunks_stats[chunk_ind].skp_stat = skp;

	// 		return;
	// 	}
	// 	else
	// 		skp++;
	// }

	// skip_chunk = skp;
	// (*freq) = available_frequencies[nb_available_frequencies - 1];
	// kmlio_chunks_stats[chunk_ind].chunk_estimated_delay = estimated_curr_chunk_delay;
	// kmlio_chunks_stats[chunk_ind].chunk_rem_checkpoint = rem_time;
}



void kmlio_diag(size_t k,size_t dim, size_t N, size_t taille,double D_max,double * centroid){
	char*  dirname ;
	asprintf(&dirname,"SKIP_CHUNKS/reports/%dN_%dM_%dD_%dK_%dL_%.2fbeta",N,taille,dim,k,(int)D_max,beta) ;
	int check = mkdir(dirname,0777);
	sleep(5);

	FILE * fl = fopen("SKIP_CHUNKS/reports/log_skip_chunk.csv","at") ;

	// log kmlio dataset properties : N, DIM , K , M , DMAX , TOTAL kmlio time
	fprintf(fl,"%d,%d,%d,%d,%.3f,%.3f,%d,%lf,",N,taille,dim,k,D_max,beta,skip_chunk,calc_delai_time(kmlio_start,kmlio_end));

	// LOG PROBLEM CONSTANT :

	fprintf(fl,"%d,%lf,%lf,%lf,%lf,%lf,{",real_time.nb_group,real_time.T_1_read,real_time.getmat_time,real_time.init_time,real_time.km_1_iteration_time,real_time.var_copy_time);

	// log each skip chunk checkpoint decision : CURRENT CHUNK , C_cp, C_cf, REMAINING TIME, REQUIERED TIME , OPTIMAL FREQ , SKIP_CHUNK , LOG EACH CHUNK REAL TIME
	char*  log_file_name ;
	asprintf(&log_file_name,"%s/log_chunks_iterations.csv",dirname) ;
	FILE * fp = fopen(log_file_name,"wt") ;

	double T_all_iteration ;
	for(int m = 0 ; m<(N/taille)+1;m++){
		T_all_iteration = 0 ;
		for(int it=0;it<kmlio_chunks_stats[m].km_nb_iterations;it++){
			fprintf(fp,"%d,%d,%lf\n",m,it,kmeans_iterations_durations[m*MAX_ITERATIONS+it]);
			T_all_iteration += kmeans_iterations_durations[m*MAX_ITERATIONS+it] ;
		}
		double chunk_delay_error = ( kmlio_chunks_stats[m].chunk_estimated_delay - kmlio_chunks_stats[m].chunk_real_delay ) / kmlio_chunks_stats[m].chunk_estimated_delay ;
		// fprintf(fl,"\"chunk_%d\":(\"nb_it\":%d;\"freq\":%ld;\"skp\":%d;T_all_it\":%f;\"remaining_time\":%f;\"estimated_delay\":%f;\"real_delay\":%f;\"delay_error\":%.4f)%s",m,kmlio_chunks_stats[m].km_nb_iterations,kmlio_chunks_stats[m].freq/1000000,T_all_iteration,kmlio_chunks_stats[m].chunk_rem_checkpoint,kmlio_chunks_stats[m].chunk_estimated_delay,kmlio_chunks_stats[m].chunk_real_delay,chunk_delay_error,(m==(N/taille)?"":",")) ;
		fprintf(fl,"\"chunk_%d\":(%d;%ld;%d;%f;%f;%f;%f;%.4f)%s",
																m,
																kmlio_chunks_stats[m].km_nb_iterations,
																kmlio_chunks_stats[m].freq/1000000,
																kmlio_chunks_stats[m].skp_stat,
																T_all_iteration,
																kmlio_chunks_stats[m].chunk_rem_checkpoint,
																kmlio_chunks_stats[m].chunk_estimated_delay,
																kmlio_chunks_stats[m].chunk_real_delay,
																chunk_delay_error,
																(m==(N/taille)?"":",")
															) ;
	}
	fprintf(fl,"}");

	// LOG_CENTERS
	char*  result_file_name ;
	asprintf(&result_file_name,"%s/result_centers.csv",dirname) ;
	r8mat_write(result_file_name,dim,k,centroid) ;
}


// //////////////////////////////////////////////////
// //////////////////////////////////////////////////
// //////////////////////////////////////////////////

struct groupe
{
	int *means;
	int nb;
	int nb_members;
	int *members;
};
typedef struct groupe groupe;

struct index_chunk_elem
{
	int element;
	struct index_chunk_elem *next;
};
typedef struct index_chunk_elem index_chunk_elem;

void r8mat_write(char *output_filename, int m, int n, double table[]);

void fail(char *str)
{
	// printf("%s", str);
	exit(-1);
}

double calc_distance(int dim, double *p1, double *p2)
{
	double distance_sq_sum = 0;
	int ii;
	for (ii = 0; ii < dim; ii++)
		distance_sq_sum += sqr(p1[ii] - p2[ii]);

	return distance_sq_sum;
}

void calc_all_distances(int dim, int n, int k, double *X, double *centroid, double *distance_output)
{
	int ii, jj;
	for (ii = 0; ii < n; ii++)	   // for each point
		for (jj = 0; jj < k; jj++) // for each cluster
		{
			// calculate distance between point and cluster centroid
			distance_output[ii * k + jj] = calc_distance(dim, &X[ii * dim], &centroid[jj * dim]);
		}
}

double calc_total_distance(int dim, int n, int k, double *X, double *centroids, int *cluster_assignment_index)
// NOTE: a point with cluster assignment -1 is ignored
{
	double tot_D = 0;
	int ii;
	// for every point
	for (ii = 0; ii < n; ii++)
	{
		// which cluster is it in?
		int active_cluster = cluster_assignment_index[ii];

		// sum distance
		if (active_cluster != -1)
			tot_D += calc_distance(dim, &X[ii * dim], &centroids[active_cluster * dim]);
	}

	return tot_D;
}

void choose_all_clusters_from_distances(int dim, int n, int k, double *distance_array, int *cluster_assignment_index)
{
	int ii, jj;
	// for each point
	for (ii = 0; ii < n; ii++)
	{
		int best_index = -1;
		double closest_distance = BIG_double;

		// for each cluster
		for (jj = 0; jj < k; jj++)
		{
			// distance between point and cluster centroid

			double cur_distance = distance_array[ii * k + jj];
			if (cur_distance < closest_distance)
			{
				best_index = jj;
				closest_distance = cur_distance;
			}
		}

		// record in array
		cluster_assignment_index[ii] = best_index;
	}
}

void calc_cluster_centroids(int dim, int n, int k, double *X, int *cluster_assignment_index, double *new_cluster_centroid)
{
	int cluster_member_count[MAX_CLUSTERS];
	int ii, jj;
	// initialize cluster centroid coordinate sums to zero
	for (ii = 0; ii < k; ii++)
	{
		cluster_member_count[ii] = 0;

		for (jj = 0; jj < dim; jj++)
			new_cluster_centroid[ii * dim + jj] = 0;
	}

	// sum all points
	// for every point
	for (ii = 0; ii < n; ii++)
	{
		// which cluster is it in?
		int active_cluster = cluster_assignment_index[ii];

		// update count of members in that cluster
		cluster_member_count[active_cluster]++;

		// sum point coordinates for finding centroid
		for (jj = 0; jj < dim; jj++)
			new_cluster_centroid[active_cluster * dim + jj] += X[ii * dim + jj];
	}

	// now divide each coordinate sum by number of members to find mean/centroid
	// for each cluster
	for (ii = 0; ii < k; ii++)
	{
		/*   if (cluster_member_count[ii] == 0)
			// printf("WARNING: Empty cluster %d! \n", ii);
			 */
		// for each dimension
		for (jj = 0; jj < dim; jj++)
			new_cluster_centroid[ii * dim + jj] /= cluster_member_count[ii]; /// XXXX will divide by zero here for any empty clusters!
	}
}

void get_cluster_member_count(int n, int k, int *cluster_assignment_index, int *cluster_member_count)
{
	int ii;
	// initialize cluster member counts
	for (ii = 0; ii < k; ii++)
		cluster_member_count[ii] = 0;
	// count members of each cluster
	for (ii = 0; ii < n; ii++)
	{
		cluster_member_count[cluster_assignment_index[ii]]++;
	}
}

void update_delta_score_table(int dim, int n, int k, double *X, int *cluster_assignment_cur, double *cluster_centroid, int *cluster_member_count, double *point_move_score_table, int cc)
{
	int ii, kk;
	// for every point (both in and not in the cluster)
	for (ii = 0; ii < n; ii++)
	{
		double dist_sum = 0;
		for (kk = 0; kk < dim; kk++)
		{
			double axis_dist = X[ii * dim + kk] - cluster_centroid[cc * dim + kk];
			dist_sum += sqr(axis_dist);
		}

		double mult = ((double)cluster_member_count[cc] / (cluster_member_count[cc] + ((cluster_assignment_cur[ii] == cc) ? -1 : +1)));

		point_move_score_table[ii * dim + cc] = dist_sum * mult;
	}
}

void perform_move(int dim, int n, int k, double *X, int *cluster_assignment, double *cluster_centroid, int *cluster_member_count, int move_point, int move_target_cluster)
{
	int cluster_old = cluster_assignment[move_point];
	int cluster_new = move_target_cluster;
	int ii;
	// update cluster assignment array
	cluster_assignment[move_point] = cluster_new;

	// update cluster count array
	cluster_member_count[cluster_old]--;
	cluster_member_count[cluster_new]++;

	if (cluster_member_count[cluster_old] <= 1)
		// printf("WARNING: Can't handle single-member clusters! \n");

		// update centroid array
		for (ii = 0; ii < dim; ii++)
		{
			cluster_centroid[cluster_old * dim + ii] -= (X[move_point * dim + ii] - cluster_centroid[cluster_old * dim + ii]) / cluster_member_count[cluster_old];
			cluster_centroid[cluster_new * dim + ii] += (X[move_point * dim + ii] - cluster_centroid[cluster_new * dim + ii]) / cluster_member_count[cluster_new];
		}
}

void copy_assignment_array(int n, int *src, int *tgt)
{
	int ii;
	for (ii = 0; ii < n; ii++)
		tgt[ii] = src[ii];
}

int assignment_change_count(int n, int a[], int b[])
{
	int change_count = 0;
	int ii;
	for (ii = 0; ii < n; ii++)
		if (a[ii] != b[ii])
			change_count++;

	return change_count;
}

void kmeans(
	int dim, // dimension of data

	double *X, // pointer to data
	int N,
	int n, // number of elements

	int k,						   // number of clusters
	double *cluster_centroid,	   // initial cluster centroids
	int *cluster_assignment_final, // output
	int *batch_iteration,
	int *change_count)
{
	double T_1_iteration, T_all_iteration;
	double *dist = (double *)malloc(sizeof(double) * n * k);
	int *cluster_assignment_prev = NULL;
	cluster_assignment_prev = (int *)malloc(sizeof(int) * n);
	int *cluster_assignment_cur = NULL;
	cluster_assignment_cur = (int *)malloc(sizeof(int) * n);
	*change_count = 0;

	// initial setup
	calc_all_distances(dim, n, k, X, cluster_centroid, dist);
	choose_all_clusters_from_distances(dim, n, k, dist, cluster_assignment_cur);
	copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);

	// BATCH UPDATE
	double prev_totD = BIG_double;
	(*batch_iteration) = 0;

	struct timeval km_start, km_end;

	gettimeofday(&km_start, NULL);

	while ((*batch_iteration) < MAX_ITERATIONS)
	{

		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////

		// update cluster centroids
		calc_cluster_centroids(dim, n, k, X, cluster_assignment_cur, cluster_centroid);

		// see if we've failed to improve
		double totD = calc_total_distance(dim, n, k, X, cluster_centroid, cluster_assignment_cur);
		if (totD > prev_totD)
		// failed to improve - currently solution worse than previous
		{
			// restore old assignments
			copy_assignment_array(n, cluster_assignment_prev, cluster_assignment_cur);
			// recalc centroids
			calc_cluster_centroids(dim, n, k, X, cluster_assignment_cur, cluster_centroid);
			// printf("  negative progress made on this step - iteration completed (%.2f) \n", totD - prev_totD);
			// done with this phase
			break;
		}

		// save previous step
		copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);
		// move all points to nearest cluster
		calc_all_distances(dim, n, k, X, cluster_centroid, dist);
		// reassign points
		choose_all_clusters_from_distances(dim, n, k, dist, cluster_assignment_cur);
		// calculate change count
		int change_count = assignment_change_count(n, cluster_assignment_cur, cluster_assignment_prev);

		// done with this phase if nothing has changed
		if (change_count == 0){	break;}

		prev_totD = totD;
		(*batch_iteration)++;

	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	gettimeofday(&km_end, NULL);
	// kmeans_iterations_durations[(chunk_ind)*MAX_ITERATIONS + (*batch_iteration) - 1] = calc_delai_time(km_start, km_end);
	if (chunk_ind == 0)
	{
		// get the average iteration delay
		real_time.km_1_iteration_time += calc_delai_time(km_start, km_end) / (*batch_iteration) ;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	kmlio_chunks_stats[chunk_ind].km_nb_iterations = (*batch_iteration);

	// write to output array
	copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_final);

	free(dist);
	dist = NULL;

	free(cluster_assignment_cur);
	cluster_assignment_cur = NULL;

	free(cluster_assignment_prev);
	cluster_assignment_prev = NULL;
}

int *mark(char *source, size_t dim, size_t N, int chunk_size)
{
	int *marks = (int *)malloc(sizeof(int) * (N / chunk_size) + 1);
	char line[100000];
	FILE *src = fopen(source, "r+");
	char delim[3] = "\t";
	char *token;
	int j = 0, mark_index = 1;
	marks[0] = ftell(src);
	// // printf ("marks[0] = %d\n", marks[0]);
	while (!feof(src))
	{
		fgets(line, 100000, src);
		// // printf("j = %d , ",j) ;
		j = (j + 1) % chunk_size;
		// if(strlen(line)>max_line_size) max_line_size = strlen(line) ;
		if (j == 0)
		{
			marks[mark_index] = ftell(src);
			// // printf ("marks[%d] = %ld , %d\n", mark_index, ftell(src),j);
			mark_index++;
		}
	}
	// // printf("max line size = %ld\n",max_line_size) ;

	fclose(src);
	return marks;
}

double *getmatrix(char *source, size_t dim, size_t N, int sub, int offset, int *marks)
{
	FILE *src = fopen(source, "r");
	char line[100000];
	double *X = NULL;
	size_t i = 0, j = 0;
	int stop = 0;

	X = (double *)malloc(sizeof(double) * sub * dim);
	fseek(src, marks[offset / sub], SEEK_CUR);
	char delim[3] = "\t";
	char *token;

	j = 0;
	while (!feof(src) && !stop)
	{
		fgets(line, 100000, src);
		token = strtok(line, delim);
		while (token != NULL)
		{
			X[j] = atof(token);
			j++;
			token = strtok(NULL, delim);
		}

		if (sub != 0)
			if ((j / (dim)) == sub)
				stop = 1;
	}
	fclose(src);

	return X;
}

double *getmatrix_buf(char *source, size_t dim, int k, size_t N, int sub, int offset, int *marks)
{

	FILE *src = fopen(source, "r");

	double *X = NULL;
	size_t i = 0, j = 0;
	int stop = 0;

	// decide best buffer size
	int MAX_DOUBLE_LEN = 24;
	int MAX_LINE_SIZE = (MAX_DOUBLE_LEN + 1) * dim;
	char line[MAX_LINE_SIZE + 1];

	// calculate the rest size allowed for the chunk kmeans
	int kmeans_size = (sizeof(double) * k + 3 * sizeof(int)) * sub;

	// the size of chunk data in file
	int chunk_text_size = marks[(offset / sub) + 1] - marks[offset / sub];
	// the average size of chunk line in file
	int avg_line_size = ceil(chunk_text_size / (double)sub);

	// the maximum number of lines to store in buffer
	int L_max = kmeans_size / (avg_line_size + 1 + sizeof(char *));
	// ajust the number of lines in buffer
	int it_num = ceil((double)sub / L_max);
	// balance the number of read lines in each iteration
	int L_opt = ceil((double)sub / it_num);

	// allocate matrix
	X = (double *)malloc(sizeof(double) * sub * dim);
	// // printf("stop after X malloc\n") ;

	// sleep(20) ;
	// begin reading the file
	fseek(src, marks[offset / sub], SEEK_CUR);
	char delim[3] = "\t";
	// var to handle strtok tokens
	char *token;
	// store the buffered lines read form files at once
	char **buf = (char **)calloc(L_opt, sizeof(char *));
	// // printf("stop after buf malloc\n") ;

	// load the chunk
	j = 0;
	int L = 0;
	int to_read_lines = sub;
	int s;
	size_t max_len;

	while (!stop)
	{
		// allocate buffer size
		if (L_opt > to_read_lines)
			L_opt = to_read_lines;
		// fill the buffer as much as possible
		L = 0;
		size_t buf_len = 0;
		while (!feof(src) && L < L_opt)
		{
			fgets(line, sizeof(line), src);
			// buf[L] = strndup(line,strlen(line)) ;
			buf[L] = calloc(strlen(line), sizeof(char));
			strncpy(buf[L], line, strlen(line));
			buf[L][strlen(line) - 1] = '\0';
			// buf_len += strlen(line) ;
			L++;
		}

		// // printf("stop after read\n") ;
		// sleep(15) ;

		// // printf("Lmax = %d \tLopt = %d \tL = %d \tread lines = %ld \tmax line size = %ld \tstrings real size %f MB + buffer pointers size %f MB\n",L_max,L_opt, L,(j/(dim)),max_line_size+1+sizeof(char *),buf_len/pow(1024,2),(L_opt*sizeof(char*))/pow(1024,2)) ;

		L = L_opt;
		to_read_lines -= L;

		// convert data strings to matrix values
		for (int l = 0; l < L; l++)
		{
			// for(int d = 0 ; d<dim ; d++){
			// 	X[j] = 0 ;
			//  	j++;
			// }
			token = strtok(buf[l], delim);
			while (token != NULL)
			{
				X[j] = atof(token);
				j++;
				token = strtok(NULL, delim);
			}
			free(buf[l]);
			buf[l] = NULL;
		}

		// s = malloc_trim(0);
		// // printf("stop after fill\n") ;
		// sleep(15) ;

		// if we have read all the chunk lines
		if (sub != 0)
		{
			if ((j / (dim)) == sub)
			{
				stop = 1;
			}
		}
	}

	if (buf != NULL)
	{
		free(buf);
		buf = NULL;
	}

	fclose(src);

	s = malloc_trim(0);

	// // printf("stop after free\n") ;
	// sleep(10) ;

	return X;
}

int token_center(int *centers, int nb, int center)
{
	int found = 0, i = 0;
	while (i < nb && !found)
	{
		if (centers[i] == center)
			found = 1;
		i++;
	}
	return found;
}

int find_best_center(double *distance_cur_center, int *centers, size_t N, size_t dim, size_t k, double sum)
{
	int best, i;
	double max = 0;
	// iterate over points
	for (i = 0; i < N; i++)
	{
		// if the distance point-center/sumof distance is bigger than max and the current center of the point is not taken
		if (distance_cur_center[i] / sum > max && !token_center(centers, k, i))
		{
			// reassign the best found center and the max distance
			best = i;
			max = distance_cur_center[i] / sum;
		}
	}
	return best;
}

double affect_x(double *X, int *centers_int, int index_x, int k, size_t dim)
{
	double min = BIG_double;
	double dist;
	int j;
	for (j = 0; j < k; j++)
	{
		dist = calc_distance(dim, &X[centers_int[j] * dim], &X[index_x * dim]);
		if (dist < min)
		{
			min = dist;
		}
	}
	return min;
}

double *kmeans_init_plusplus(double *X, size_t N, size_t dim, size_t k)
{
	double *centers = (double *)malloc(sizeof(double) * k * dim);
	double *distance_cur_center = (double *)malloc(sizeof(double) * N);
	int *centers_int = (int *)malloc(sizeof(int) * k);
	double sum = 0;
	int first = rand() % N;
	// int first = chunk_ind;

	int i, j, best;
	centers_int[0] = first;
	for (i = 1; i < k; i++)
	{

		for (j = 0; j < N; j++)
		{
			distance_cur_center[j] = affect_x(X, centers_int, j, i, dim);
			sum += distance_cur_center[j];
		}

		centers_int[i] = find_best_center(distance_cur_center, centers_int, N, dim, i, sum);
	}
	for (i = 0; i < k; i++)
	{
		for (j = 0; j < dim; j++)
		{
			centers[i * dim + j] = X[centers_int[i] * dim + j];
			//	// printf ("%lf\t", X[centers_int[i]*dim+j]);
		}
		// // printf ("\n");
	}
	free(centers_int);
	centers_int = NULL;
	free(distance_cur_center);
	distance_cur_center = NULL;

	return centers;
}

// partial kmeans function
double *partial_mean(double *X, size_t dim, size_t chunk_len, size_t chunk_index)
{
	int i, j, k;
	double *sum = malloc(sizeof(double) * dim);

	for (i = 0; i < dim; i++)
		sum[i] = 0;
	for (i = 0; i < chunk_len; i++)
	{
		for (j = 0; j < dim; j++)
		{
			k = (chunk_index * dim) + (i * dim) + j;
			sum[j] = sum[j] + X[k];
		}
	}

	for (i = 0; i < dim; i++)
	{
		sum[i] /= chunk_len;
		// ("\t %lf", sum[i]);
	}
	// // printf("\n");

	return sum;
}

double sse_calculate(double *X, int dim, double *centroids, int *assignement, size_t n)
{
	double sse = 0;
	int i;
	for (i = 0; i < n; i++)
	{
		sse += calc_distance(dim, &centroids[assignement[i] * dim], &X[i * dim]);
	}
	return sse;
}

double var_calculate(double *X, int dim, double *centroids, int center, int *assignement, size_t n)
{
	double sse = 0;
	int i;
	for (i = 0; i < n; i++)
	{
		if (assignement[i] == center)
		{
			sse += calc_distance(dim, &centroids[center * dim], &X[i * dim]);
		}
	}
	return sse;
}

double inter_centroid_distance(int dim, int k, double *centroids)
{
	double sum = 0;
	int i, j;
	for (i = 0; i < k; i++)
	{
		for (j = i + 1; j < k; j++)
		{
			sum += calc_distance(dim, &centroids[i * dim], &centroids[j * dim]);
		}
	}
	return sum;
}

double *getsubmatrix(double *X, int dim, int offset, int len)
{
	double *Y = (double *)malloc((dim * len) * sizeof(double));
	int i;
	for (i = offset * dim; i < dim * len + offset * dim; i++)
	{
		Y[i - offset * dim] = X[i];
	}
	return Y;
}

int existJ(int k, int *tab, int j, int max)
{
	int i = 0;
	while (i < k && tab[i] != j && i < max)
	{
		i++;
	}
	if (i < k && i != max)
		return 1;
	else
		return 0;
}

int maximize(double *partial_means, int *jpris, int pris, int dim, int k, int chunks)
{
	int i, m, res;
	double sum, maxsum;
	maxsum = 0;
	res = 0;
	for (i = 0; i < chunks; i++)
	{ // pour tout chunk
		sum = 0;

		if (!existJ(k, jpris, i, pris))
		{ // s'il fait pas partie des chunks sélectionnés
			// for (j=0; j<chunks;j++){
			for (m = 0; m < pris; m++)
			{
				sum += calc_distance(dim, &partial_means[i * dim], &partial_means[jpris[m] * dim]);
			}
			//}
			if (sum > maxsum)
			{
				maxsum = sum;
				res = i;
			}
		}
	}

	// printf("stop %d\n", res);
	return res;
}

int existG(groupe *groupes, int nb, int g)
{
	int found = 0, i = 0, j;
	// iterate over the groupes
	while (i < nb && !found)
	{
		j = 0;
		// iterate over the chunks of the group
		while (j < groupes[i].nb && !found)
		{
			if (groupes[i].means[j] == g)
			{
				found = 1;
			}
			j++;
		}
		i++;
	}
	return found;
}

int existGi(groupe groupes, int g)
{
	int found = 0, j = 0;
	while (j < groupes.nb && !found)
	{
		if (groupes.means[j] == g)
		{
			found = 1;
		}
		j++;
	}
	return found;
}

// version 1: pour séléctionner aléatoirement un groupe parcourir toute la structure cluster_assignement jusqu'à tomber sur le point qui correspond au cluster
void form_index(int *cluster_assignement, int nb_groupes, groupe *grp, size_t N)
{
	int i = 0, j = 0, affected, *index_grp;
	index_grp = (int *)malloc(sizeof(int) * nb_groupes + 1);
	for (i = 0; i < nb_groupes; i++)
	{
		grp[i].members = (int *)malloc(sizeof(int) * grp[i].nb_members);
		index_grp[i] = 0;
	}

	// iterate over points
	for (i = 0; i < N; i++)
	{
		affected = 0;
		j = 0;
		// affect the point to a group
		while (!affected && j < nb_groupes)
		{
			// if the point cluster is exist in the group j
			if (existGi(grp[j], cluster_assignement[i]))
			{
				grp[j].members[index_grp[j]] = i;
				index_grp[j]++;
				affected = 1;
			}
			j++;
		}
	}

	free(index_grp);
}

double *form_chunk(groupe *grp, /*double *X*/ char *source, int *marks, int *cluster_assignment, int nb_groupes, size_t N, size_t dim, int k, size_t taille)
{
	double *chunk = (double *)malloc(sizeof(double) * (dim * taille));
	int nb_samples;
	int j, index_groupe, l, size = 0, i,*int_chunk; // found,
	int tmp;
	double *Y;
	float c;
	char **buf;

	index_chunk_elem index[((N / taille) - skip_chunk)];
	index_chunk_elem *cur[((N / taille) - skip_chunk)];
	index_chunk_elem *elem;
	for (i = 0; i < ((N / taille) - skip_chunk); i++)
	{
		index[i].element = -1;
		index[i].next = NULL;
		cur[i] = NULL;
	}
	// form groups index
	// // printf("BEGIN BUILD INDEX\n") ;
	form_index(cluster_assignment, nb_groupes, grp, N - (taille * skip_chunk));

	// form chunk
	while (size < taille)
	{
		// iterate over groups
		for (i = 0; i < nb_groupes; i++)
		{
			// get the elements ratio that can be in chunk for the group
			c = (float)grp[i].nb_members / (float)N * taille;
			nb_samples = (int)c;
			j = 0;
			// choose the elements in the chunk randomly
			for (j = 0; j < nb_samples; j++)
			{

				tmp = rand() % grp[i].nb_members;
				// tmp = (j + i * nb_samples) % grp[i].nb_members;

				if (size == taille)
					break;
				// if the slected point from the group does not exist
				// in the index , add it
				if (index[grp[i].members[tmp] / taille].element == -1)
				{
					index[grp[i].members[tmp] / taille].element = grp[i].members[tmp] % taille;
				}
				else
				{

					elem = (index_chunk_elem *)malloc(sizeof(index_chunk_elem));
					elem->element = grp[i].members[tmp] % taille;

					elem->next = NULL;
					if (cur[grp[i].members[tmp] / taille] != NULL)
						cur[grp[i].members[tmp] / taille]->next = elem;
					else
					{
						cur[grp[i].members[tmp] / taille] = (index_chunk_elem *)malloc(sizeof(index_chunk_elem));
						index[grp[i].members[tmp] / taille].next = elem;
					}
					cur[grp[i].members[tmp] / taille] = elem;
				}
				// free(elem);

				//}
				size++;
			}
		}
	}

	int found = 0;
	double estimated_chunk_delay = 0, rem_time = 0;
	long freq = get_optimal_freq(N, taille,&estimated_chunk_delay, &rem_time);

	kmlio_chunks_stats[chunk_ind].chunk_estimated_delay = estimated_chunk_delay;
	kmlio_chunks_stats[chunk_ind].chunk_rem_checkpoint = rem_time;
	kmlio_chunks_stats[chunk_ind].freq = freq;

	set_frequency(freq);

	size = 0;
	// iterate over chunk points
	for (i = 0; i < ((N / taille) - skip_chunk); i++)
	{
		if (index[i].element != -1)
		{
			// Y = getmatrix(source, dim, taille,taille,i*taille, marks);
			Y = getmatrix_buf(source, dim, k, taille, taille, i * taille, marks);
			// // printf ("chunk %d\n", i);
			for (j = 0; j < dim; j++)
			{
				chunk[size * dim + j] = Y[index[i].element * dim + j];
				// // printf("%lf\n", Y[index[i].element*dim+j]);
			}
			size++;
			// goto next point in index
			cur[i] = index[i].next;
			while (cur[i] != NULL)
			{
				// // printf("ELEMENT %ld\n", cur[i]->element);
				for (j = 0; j < dim; j++)
				{
					chunk[size * dim + j] = Y[cur[i]->element * dim + j];
				}
				cur[i] = cur[i]->next;
				size++;
			}
			free(cur[i]);
			free(Y);
		}
	}
	// // printf ("size %d\n", size);
	return chunk;
}

void kmeans_by_chunk(char *source, size_t dim, int taille, int N, int k, double **cluster_centroid)
{
	int i, j, n, kmeans_iterations = 0, kmeans_change_count = 0;
	double *Y = NULL;
	int *cluster_assignment_final = NULL, *cluster_assignment_final_Y = NULL;
	double *cluster_centroid_by_chunk = NULL, *chunk_centroid = NULL;
	double *var = NULL;
	int cluster_member_count[MAX_CLUSTERS]; // à modifier pour une allocation dynamique
	int nb_groupes = 0;
	groupe *groupes = NULL;

	// //////////////////////////////////////////
	struct timeval tv_1, tv_2, chunk_start, chunk_end;
	long freq;
	// //////////////////////////////////////////
	// *cluster_centroid = NULL ;

	/****** PHASE 1 : PARTIELS CHUNKS KMEANS ******/
	set_frequency(available_frequencies[0]);
	gettimeofday(&tv_1, NULL);
	// mark the chunks in dataset
	int *marks = mark(source, dim, N, taille);
	gettimeofday(&tv_2, NULL);
	real_time.mark_time = calc_delai_time(tv_1, tv_2);

	set_frequency(available_frequencies[nb_available_frequencies - 1]);

	// apply kmeans on each chunk
	if (((N / taille) - skip_chunk) > 1)
	{
		cluster_assignment_final = (int *)malloc(N * sizeof(int));
		cluster_centroid_by_chunk = (double *)malloc(sizeof(double) * (k * dim * ((N / taille) - skip_chunk)));
		var = (double *)malloc(sizeof(double) * k * ((N / taille) - skip_chunk));

		//////////////////////////////////////

		kmeans_iterations_durations = calloc(MAX_ITERATIONS * ((N / taille) + 1), sizeof(double));
		kmlio_chunks_stats = calloc((N / taille) + 1, sizeof(chunk_stats));

		//////////////////////////////////////
		chunk_ind = 0;
		while (chunk_ind < (N / taille) - skip_chunk)
		{
			gettimeofday(&chunk_start, NULL);
			//////////////////////////////////////

			// lecture du chunk
			// kmlio_chunks_stats[chunk_ind] = (chunk_stats*) calloc(1,sizeof(chunk_stats)) ;

			i = chunk_ind * taille;
			Y = getmatrix_buf(source, dim, k, taille, taille, i, marks);

			if (chunk_ind == 0) //SAVEPOINT:
			{
				gettimeofday(&tv_1, NULL);
				real_time.getmat_time = calc_delai_time(chunk_start, tv_1);
			}

			//////////////////////////////////////

			// chunk chunk_ind init kmeans++
			(*cluster_centroid) = kmeans_init_plusplus(Y, taille, dim, k);

			if (chunk_ind == 0)	 //SAVEPOINT:
			{
				gettimeofday(&tv_2, NULL);
				real_time.init_time = calc_delai_time(tv_1, tv_2);
				real_time.km_1_iteration_time = 0;
			}

			//////////////////////////////////////

			// apply kmeans on the chunk chunk_ind
			cluster_assignment_final_Y = (int *)malloc(taille * sizeof(int));
			kmeans(dim, Y, N, taille, k, (*cluster_centroid), cluster_assignment_final_Y, &kmeans_iterations, &kmeans_change_count);

			//////////////////////////////////////

			//////////////////////////////////////

			gettimeofday(&tv_1, NULL);

			// variance calculation on the chunk chunk_ind
			for (j = 0; j < k; j++)
			{
				var[chunk_ind * k + j] = var_calculate(Y, dim, (*cluster_centroid), j, cluster_assignment_final_Y, taille) / taille;
			}
			// copy the centroids to global array
			for (j = 0; j < k * dim; j++)
			{
				cluster_centroid_by_chunk[chunk_ind * k * dim + j] = (*cluster_centroid)[j];
			}
			// copy points clustering to global array
			for (n = 0; n < taille; n++)
			{
				cluster_assignment_final[chunk_ind * taille + n] = cluster_assignment_final_Y[n] + chunk_ind * k;
			}

			free(Y);
			Y = NULL;

			free((*cluster_centroid));
			(*cluster_centroid) = NULL;

			free(cluster_assignment_final_Y);
			cluster_assignment_final_Y = NULL;

			gettimeofday(&chunk_end, NULL);
			if (chunk_ind == 0)	 //SAVEPOINT:
				real_time.var_copy_time = calc_delai_time(tv_1, chunk_end);


			//SAVEPOINT: svae chunk real elapsed time
			kmlio_chunks_stats[chunk_ind].chunk_real_delay = calc_delai_time(chunk_start, chunk_end);

			////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////

			if( chunk_ind == 0 ){	//CHECKPOINT : DECIDE THE NUMBER OF CHUNKS TO DROP
				// calculate the observation times at each opration of kmlio
				calc_cycles_per_elem(N, n);

				// get optimal frequency and skip chunk for the first chunk execution
				decide_skip_chunk(N, n);

				// get optimal freq
				printf("chunks to drop = %d \n", skip_chunk);

				// ////////////////////////////////////
			}

			//check for the next chunk
			chunk_ind++;

			if (((N / taille) - chunk_ind - skip_chunk) <= 0) //CHECKPOINT : if processing the remaining partial chunks is complete
			{
				if(chunk_ind==1){ // if we performed only the first chunk , kmlio return without composing the final chunk
					// free all computations structures
					free(Y);
					Y = NULL;
					free(cluster_assignment_final_Y);
					cluster_assignment_final_Y = NULL;
					free(marks);
					marks = NULL;
					
					// KMLIO IS TERMINATED
					return ;
				}
				//else , we STOP PROECESSING PARTIAL CHUNKS and go the kmlio's next phase
				break;
			}
			else
			{
				printf("same skip chunk=%d get optimal freq for chunk %d\n", skip_chunk, chunk_ind);
				int found = 0;
				double estimated_chunk_delay = 0, rem_time = 0;

				freq = get_optimal_freq(N, taille, &estimated_chunk_delay, &rem_time);

				kmlio_chunks_stats[chunk_ind].chunk_estimated_delay = estimated_chunk_delay ;
				kmlio_chunks_stats[chunk_ind].chunk_rem_checkpoint = rem_time ;
				kmlio_chunks_stats[chunk_ind].freq = freq ;
				kmlio_chunks_stats[chunk_ind].skp_stat = skip_chunk ;

				set_frequency(freq) ;
			}

			////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////

		}// PARTIAL CHUNKS CLUSTERING END

		set_frequency(available_frequencies[nb_available_frequencies - 1]);

		/****** PHASE 2 : PARTIELS CLUSTERS GROUPING ******/

		groupes = (groupe *)malloc(sizeof(groupe) * ((N / taille) - skip_chunk) * k);
		nb_groupes = 0;
		for (j = 0; j < k * ((N / taille) - skip_chunk); j++)
		{
			if (!existG(groupes, nb_groupes, j))
			{
				groupes[nb_groupes].means = (int *)malloc(sizeof(int));
				groupes[nb_groupes].nb = 1;
				groupes[nb_groupes].means[0] = j;
				for (n = j + 1; n < k * ((N / taille) - skip_chunk); n++)
				{
					// condition de non existence
					if (calc_distance(dim, &cluster_centroid_by_chunk[j * dim], &cluster_centroid_by_chunk[n * dim]) < var[j] && !existG(groupes, nb_groupes, n))
					{
						groupes[nb_groupes].means = (int *)realloc(groupes[nb_groupes].means, (groupes[nb_groupes].nb + 1) * sizeof(int));
						groupes[nb_groupes].means[groupes[nb_groupes].nb] = n;
						// // printf ("next\n");
						groupes[nb_groupes].nb++;
					}
				}
				nb_groupes++;
			}
		}

		free(cluster_centroid_by_chunk);
		cluster_centroid_by_chunk = NULL;

		// groups members count update
		float count = k * ((N / taille) - skip_chunk);
		get_cluster_member_count(((int)((N / taille) - skip_chunk)) * taille, (int)count, cluster_assignment_final, cluster_member_count);

		for (j = 0; j < nb_groupes; j++)
		{
			groupes[j].nb_members = 0;
			for (n = 0; n < groupes[j].nb; n++)
			{
				groupes[j].nb_members += cluster_member_count[groupes[j].means[n]];
			}
		}

		real_time.nb_group = nb_groupes ; //SVAEPOINT : GROUPS NUMBER

		/****** PHASE 3 : FINAL CHUNK BUILDING ******/

		double *chunk;
		chunk = form_chunk(groupes, source, marks, cluster_assignment_final, nb_groupes, N, dim, k, taille);

		free(cluster_assignment_final);
		cluster_assignment_final = NULL;
		for (i = 0; i < nb_groupes; i++)
		{
			free(groupes[i].means);
			free(groupes[i].members);
		}
		free(groupes);
		groupes = NULL;

		/****** PHASE 4 : FINAL CHUNK KMEANS ******/

		// // printf("begin final chunk\n") ;

		(*cluster_centroid) = kmeans_init_plusplus(chunk, taille, dim, k);

		cluster_assignment_final_Y = (int *)malloc(taille * sizeof(int));
		kmeans(dim, chunk, N, taille, k, (*cluster_centroid), cluster_assignment_final_Y, &kmeans_iterations, &kmeans_change_count);

		free(chunk);
		chunk = NULL;
	}
	else
	{
		double *X;

		X = getmatrix_buf(source, dim, k, N, N, 0, marks);

		(*cluster_centroid) = kmeans_init_plusplus(X, N, dim, k);

		cluster_assignment_final = (int *)malloc(N * sizeof(int));
		kmeans(dim, X, N, N, k, (*cluster_centroid), cluster_assignment_final, &kmeans_iterations, &kmeans_change_count);

		free(X);
		X = NULL;
	}

	free(cluster_assignment_final_Y);
	cluster_assignment_final_Y = NULL;

	free(marks);
	marks = NULL;
}

int main(int argc, char **argv)
{
	cpu_set_t set;
	// clear cpu mask
	CPU_ZERO(&set);
	// set cpu 4
	CPU_SET(4, &set);
	// 0 is the calling process
	sched_setaffinity(0, sizeof(cpu_set_t), &set);
	// set priority
	setpriority(PRIO_PROCESS, 0, -20);

	clock_t begin, end;
	int i, j;
	size_t tot;
	int taille;
	size_t dim, N;
	double *X;

	char cmd[200];
	sprintf(cmd, "echo %d > /sys/fs/cgroup/kmeans/cgroup.procs", getpid()); // /cgroups/mem/kmeans/tasks
	system(cmd);

	sleep(1);

	char *source = argv[1];
	size_t k = atoi(argv[2]);
	N = atoi(argv[3]);
	taille = atoi(argv[4]);
	dim = atoi(argv[5]);
	D_max = atoi(argv[6]);	// kmlio maximal delay time in seconds

	// nb_it_sample = (argc > 7 ? atoi(argv[7]) : 1); // number of samples that enter in calculating iteration avg time
	// beta = (argc > 8 ? atof(argv[8]) : 0);		   // precison - energy tendency factor for skip chunk strategy

	skip_chunk = 0;

	srand(time(NULL));

	double *centroid = NULL;

	set_governor("userspace");

	gettimeofday(&kmlio_start, NULL);

	kmeans_by_chunk(source, dim, taille, N, k, &centroid);

	gettimeofday(&kmlio_end, NULL);

	set_governor("conservative");
	
	printf("KMLIO EXECUTION COMPLETE , remaining time = %f \n",calc_remaining_Time()) ;

	printf("\n----------------------------------------------------------------------------\n") ;
	printf("\n----------------------------------------------------------------------------\n") ;
	printf("\n----------------------------------------------------------------------------\n") ;

	

	kmlio_diag(k,dim,N,taille,D_max,centroid) ;

	printf("KMLIO EXECUTION REPORT READY\n") ;

	return 0;
}