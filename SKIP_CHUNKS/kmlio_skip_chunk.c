// kmeans.c
// Ethan Brodsky
// modified by Camélia SLIMANI

/*

DELL : cd /home/hafsa/Documents/@K-MLIO_Analysis/
BeagleBone : cd /home/debian/@K-MLIO_Analysis/
scripts/prog_script_cgroup 1155
scripts/prog_script_reset

(scripts/prog_script_cache) && (clear && gcc -g SKIP_CHUNKS/kmlio_skip_chunk.c -o SKIP_CHUNKS/kmlio_skip_chunk -lm -D _GNU_SOURCE) && (SKIP_CHUNKS/kmlio_skip_chunk generator/10D/13421800N/SEP-0.6/points.csv 10 13421800 6710900 10 0)

(scripts/prog_script_reset) && (clear && gcc -g kmlio.c -o program -lm -D _GNU_SOURCE) && (./program generator/CM13,4M_2400MO_SEP0,2/points.csv 10 13421800 6710900 10)
(scripts/prog_script_reset) && (clear && gcc -g kmlio.c -o program -lm -D _GNU_SOURCE) && ((./program generator/CM13,4M_2400MO_SEP0,2/points.csv 10 13421800 6710900 10) & (taskset -c 1 scripts/prog_script_launch))

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

// //////////////////////////////////////////////////
// //////////////////////////////////////////////////
// //////////////////////////////////////////////////

// double estimate_required_time(kmlio_time_estimation t_est,int freq){
	// double fre_req = D_max - calc_past_time() - (2 * (N/M) - chunk_ind - skip_chunk ) * M * real_time.C_1_it_elem ;
// 	return -1 ;
// }

struct kmlio_time_estimation
{
	double mark_time ;
	double getmat_time ;
	double init_time ;
	double form_chunk_time ;
	double km_1_iteration_time ;
	double var_copy_time ;

	double C_1_it_elem ;
	double C_init_elem ;
	double C_get_matrix_on_elem ;
	double T_1_read ;
};
typedef struct kmlio_time_estimation kmlio_time_estimation;

int chunk_ind = 0;

int skp_chk = 0 ;
int D_max ;

struct timeval kmlio_start, kmlio_end;

double* kmeans_iterations_durations ;
int* kmeans_nb_iteration ;
double* chunks_elapsed_durations ;

long MAX_FREQ = 1800000000 , MIN_FREQ = 400000000 , BASE_FREQ = 1600000000 ;

long FREQ_STEP = 100000000 ;

long set_speed ;

kmlio_time_estimation real_time ;



long * km_freq = NULL ;

double calc_delai_time(struct timeval tv_start,struct timeval tv_end){
	return ((double)(tv_end.tv_usec - tv_start.tv_usec) / 1000000 + (double)(tv_end.tv_sec - tv_start.tv_sec) ) ;
}

double calc_past_time(){
	struct timeval tv_now ;
	gettimeofday(&tv_now, NULL);
	return calc_delai_time(kmlio_start,tv_now) ;
}

double calc_remaining_Time(){
	return ( D_max - calc_past_time() ) ;
}

double estimate_T_get_matrix_time(size_t M,long freq){
	return ( (real_time.T_1_read + real_time.C_get_matrix_on_elem/freq) * M) ;
}

double estimate_T_init_time(size_t M,long freq){
	return (real_time.C_init_elem * M / freq ) ;
}

double estimate_T_1_iteration_time(size_t M,long freq){
	return (real_time.C_1_it_elem * M / freq ) ;
}

double estimate_tcp_max(size_t M,long freq){
	return (estimate_T_get_matrix_time(M,freq) + estimate_T_init_time(M,freq) + ((MAX_ITERATIONS+1) * estimate_T_1_iteration_time(M,freq)) ) ;
}

double estimate_tcf_max(size_t N,size_t M,long freq){
	return (estimate_T_get_matrix_time(M,freq) * (N/M) + estimate_T_init_time(M,freq) + ((MAX_ITERATIONS+1) * estimate_T_1_iteration_time(M,freq) ) ) ;
}

double calc_cycles_per_elem(size_t N,size_t M){
	real_time.T_1_read = real_time.mark_time / N ;
	real_time.C_get_matrix_on_elem = (real_time.getmat_time - real_time.T_1_read * M ) * set_speed / M ;
	real_time.C_init_elem = real_time.init_time * set_speed / M  ;
	real_time.C_1_it_elem = real_time.km_1_iteration_time * set_speed / M  ;
}

long get_optimal_freq(size_t N, size_t M,int skip_chunk, int* found){
	if(found)
		(*found) = 0 ;
	long freq = BASE_FREQ ;

	double first_chunk_km = (  chunk_ind == 0  ? ((MAX_ITERATIONS-1) * real_time.C_1_it_elem * M / freq)  : 0 ) ;
		
	double rem_time = calc_remaining_Time() ;
	double tcp = estimate_tcp_max(M,freq) ;
	double tcf = estimate_tcf_max(N,M,freq) ;
	int exec_chunk_nb = (N/M) - (chunk_ind == 0?chunk_ind+1:chunk_ind)- skip_chunk ;
	double req_time = exec_chunk_nb * tcp + first_chunk_km + tcf ;
	
	printf("chunk_id=%d, test skip=%d : freq %ld , required time = (%f)1st chunk + (%d)*tcp(%f)+ tcf(%f) = %f remaining time = %f\n",chunk_ind,skip_chunk,freq,first_chunk_km,exec_chunk_nb,tcp,tcf,req_time,rem_time) ;

	if ( req_time < rem_time){
		if (found)	
			(*found) = 1 ;
		return freq ;
	}

	return freq ;
}

void decide_skip_chunk(size_t N, size_t M,long * freq){
	(*freq) = BASE_FREQ  ;
	int skp = 0 , found = 0 ;

	while(skp < (N/M) - chunk_ind){
		(*freq) = get_optimal_freq(N,M,skp,&found) ;
		if( found == 1 ){
			skp_chk = skp ;
			printf("found optimal frequency = %ld\n",(*freq)) ;
			return ;
		}
		else {
			skp++ ;
		}
	}
	skp_chk = skp ;
}

void kmlio_diag(size_t k,size_t dim, size_t N, size_t taille,double D_max,double * centroid){
	char*  dirname ;
	asprintf(&dirname,"SKIP_CHUNKS/reports/%ldN_%ldM_%ldD_%ldK_%dL",N,taille,dim,k,(int)D_max) ;
	int check = mkdir(dirname,0777);
	sleep(1);

	FILE * fl = fopen("logs/log_skip_chunk.csv","at") ;

	// log kmlio dataset properties : N, DIM , K , M , DMAX , TOTAL kmlio time
	fprintf(fl,"%ld,%ld,%ld,%ld,%lf,%d,%lf,{",N,taille,dim,k,D_max,skp_chk,calc_delai_time(kmlio_start,kmlio_end));

	// LOG PROBLEM CONSTANT :

	fprintf(fl,"%lf,%lf,%lf,%lf,%lf,{",real_time.T_1_read,real_time.getmat_time,real_time.init_time,real_time.km_1_iteration_time,real_time.var_copy_time);

	// log each skip chunk checkpoint decision : CURRENT CHUNK , tcp , tcf , REMAINING TIME, REQUIERED TIME , OPTIMAL FREQ , SKIP_CHUNK , LOG EACH CHUNK REAL TIME

	FILE * fp = fopen(strcat(dirname,"/log_chunks_iterations.csv"),"wt") ;

	int T_all_iteration ;
	for(int m = 0 ; m<(N/taille)+1;m++){
		T_all_iteration = 0 ;
		for(int it=0;it<kmeans_nb_iteration[m];it++){
			fprintf(fp,"%d,%d,%lf\n",m,it,kmeans_iterations_durations[m*MAX_ITERATIONS+it]);
			T_all_iteration += kmeans_iterations_durations[it] ;
		}
		double estimated_chunk_delay = ( m==0 ? ((MAX_ITERATIONS-1) * real_time.C_1_it_elem * taille / kmeans_nb_iteration[m]) : ( m==(N/taille) ? estimate_tcf_max(N,taille,km_freq[m]) : estimate_tcp_max(taille,km_freq[m]) ) ) ;
		fprintf(fl,"(%d,%ld,%f,%f,%f),",kmeans_nb_iteration[m],km_freq[m]/1000000,T_all_iteration,estimated_chunk_delay,chunks_elapsed_durations[m]);
	}
	fprintf(fl,"\b}");

	// LOG_CENTERS
	
	r8mat_write (strcat(dirname,"/result_centers.csv"),dim,k,centroid) ;
}

// //////////////////////////////////////////////////
// //////////////////////////////////////////////////
// //////////////////////////////////////////////////


void fail(char *str)
{
	printf("%s", str);
	exit(-1);
}

/**
* Invokes the command source("foo.R").
*/
// void source(const char *name)
// {
//     SEXP e;
 
//     PROTECT(e = lang2(install("source"), mkString(name)));
//     R_tryEval(e, R_GlobalEnv, NULL);
//     UNPROTECT(1);
// }

// <pre>cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_setspeed</pre>

// os.system("echo "+governor+" | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor")
// # os.system("cpupower -c all frequency-set -f "+str(freq))
// # os.system("cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_setspeed")


void set_frequency(long freq)
{
	set_speed = freq ;
	// printf("set speed in hz %ld in khz %ld\n",set_speed,(set_speed/1000));
	system("echo userspace | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor") ;
	char *cmd ;
	asprintf(&cmd, "echo %ld | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_setspeed",(set_speed/1000)) ;
	system(cmd) ;

	// free(cmd) ;
	cmd = NULL ;
}

void set_governor()
{
	system("echo conservative | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor") ;
}

void save_phase_time(int phase, int step, int loop, int iteration)
{
	char *cmd;
	// printf("iterations=%d",iteration);
	asprintf(&cmd, "taskset -c $(($(nproc) - 1)) ./scripts/prog_script_phase %d %d %d %d -1 -1", phase, step, loop, iteration);
	system(cmd);
	free(cmd);
	cmd = NULL;
}

void save_kmeans_iterations(int phase, int step, int loop, int iteration, int change_count, double tot_D)
{
	char *cmd;
	// printf("iterations=%d",iteration);
	asprintf(&cmd, "taskset -c $(($(nproc) - 1)) ./scripts/prog_script_phase %d %d %d %d %d %f", phase, step, loop, iteration, change_count, tot_D);
	system(cmd);
	free(cmd);
	cmd = NULL;
}

double *r8mat_data_read(char *input_filename, int m, int n)
/******************************************************************************/
/*
  Purpose:

	R8MAT_DATA_READ reads the data from an R8MAT file.

  Discussion:

	An R8MAT is an array of R8's.

	The file is assumed to contain one record per line.

	Records beginning with the '#' character are comments, and are ignored.
	Blank lines are also ignored.

	Each line that is not ignored is assumed to contain exactly (or at least)
	M real numbers, representing the coordinates of a point.

	There are assumed to be exactly (or at least) N such records.

  Licensing:

	This code is distributed under the GNU LGPL license.

  Modified:

	27 January 2005

  Author:

	John Burkardt

  Parameters:

	Input, char *INPUT_FILENAME, the name of the input file.

	Input, int M, the number of spatial dimensions.

	Input, int N, the number of points.  The program
	will stop reading data once N values have been read.

	Output, double R8MAT_DATA_READ[M*N], the data.
*/
{
#define MY_LINE_MAX 255

	int error;
	char *got_string;
	FILE *input;
	int i;
	int j;
	char line[255];
	double *table;
	double *x;

	input = fopen(input_filename, "r");

	if (!input)
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "R8MAT_DATA_READ - Fatal error!\n");
		fprintf(stderr, "  Could not open the input file: \"%s\"\n", input_filename);
		exit(1);
	}

	table = (double *)malloc(m * n * sizeof(double));

	x = (double *)malloc(m * sizeof(double));

	j = 0;

	while (j < n)
	{
		got_string = fgets(line, MY_LINE_MAX, input);

		if (!got_string)
		{
			break;
		}

		// if ( line[0] == '#' || s_len_trim ( line ) == 0 ){continue;}

		char delim[3] = "\t";

		char *token = strtok(line, delim);
		i = 0;
		while (token != NULL)
		{
			x[i] = atof(token);
			i++;
			token = strtok(NULL, delim);
		}

		// error = s_to_r8vec ( line, m, x );

		// if ( error == 1 ){continue;}

		for (i = 0; i < m; i++)
		{
			table[i + j * m] = x[i];
		}
		j = j + 1;
	}

	fclose(input);

	free(x);

	return table;

#undef MY_LINE_MAX
}
/******************************************************************************/

void r8mat_write(char *output_filename, int m, int n, double table[])

/******************************************************************************/
/*
  Purpose:

	R8MAT_WRITE writes an R8MAT file.

  Discussion:

	An R8MAT is an array of R8's.

  Licensing:

	This code is distributed under the GNU LGPL license.

  Modified:

	01 June 2009

  Author:

	John Burkardt

  Parameters:

	Input, char *OUTPUT_FILENAME, the output filename.

	Input, int M, the spatial dimension.

	Input, int N, the number of points.

	Input, double TABLE[M*N], the data.
*/
{
	int i;
	int j;
	FILE *output;
	/*
	  Open the file.
	*/
	output = fopen(output_filename, "at");

	if (!output)
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "R8MAT_WRITE - Fatal error!\n");
		fprintf(stderr, "  Could not open the file '%s'.\n", output_filename);
		exit(1);
	}
	/*
	  Write the data.
	*/
	for (j = 0; j < n; j++)
	{
		for (i = 0; i < m; i++)
		{
			fprintf(output, "%2.16f\t", table[i + j * m]); //%24.16g
		}
		fprintf(output, "\n");
	}
	/*
	  Close the file.
	*/
	fclose(output);

	return;
}
/******************************************************************************/

void i4mat_write(char *output_filename, int m, int n, int table[])

/******************************************************************************/
/*
  Purpose:

	I4MAT_WRITE writes an I4MAT file.

  Discussion:

	An I4MAT is an array of I4's.

  Licensing:

	This code is distributed under the GNU LGPL license.

  Modified:

	01 June 2009

  Author:

	John Burkardt

  Parameters:

	Input, char *OUTPUT_FILENAME, the output filename.

	Input, int M, the spatial dimension.

	Input, int N, the number of points.

	Input, int TABLE[M*N], the data.
*/
{
	int i;
	int j;
	FILE *output;
	/*
	  Open the file.
	*/
	output = fopen(output_filename, "wt");

	if (!output)
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "I4MAT_WRITE - Fatal error!\n");
		fprintf(stderr, "  Could not open the output file '%s'\n", output_filename);
		exit(1);
	}
	/*
	  Write the data.
	*/
	for (j = 0; j < n; j++)
	{
		for (i = 0; i < m; i++)
		{
			fprintf(output, "  %d", table[i + j * m]);
		}
		fprintf(output, "\n");
	}
	/*
	  Close the file.
	*/
	fclose(output);

	return;
}
/******************************************************************************/

void rbinmat_write(char *output_filename, int m, int n, double table[])

/******************************************************************************/
/*
 */
{
	int i;
	int j;
	FILE *output;
	/*
	  Open the file.
	*/
	output = fopen(output_filename, "at");

	if (!output)
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "R8MAT_WRITE - Fatal error!\n");
		fprintf(stderr, "  Could not open the file '%s'.\n", output_filename);
		exit(1);
	}
	/*
	  Write the data.
	*/
	for (j = 0; j < n; j++)
	{
		for (i = 0; i < m; i++)
		{
			fwrite(&table[i + j * m], sizeof(double), 1, output); //%24.16g
		}
	}
	/*
	  Close the file.
	*/
	fclose(output);

	return;
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
			printf("WARNING: Empty cluster %d! \n", ii);
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
		printf("WARNING: Can't handle single-member clusters! \n");

	// update centroid array
	for (ii = 0; ii < dim; ii++)
	{
		cluster_centroid[cluster_old * dim + ii] -= (X[move_point * dim + ii] - cluster_centroid[cluster_old * dim + ii]) / cluster_member_count[cluster_old];
		cluster_centroid[cluster_new * dim + ii] += (X[move_point * dim + ii] - cluster_centroid[cluster_new * dim + ii]) / cluster_member_count[cluster_new];
	}
}

void cluster_diag(int dim, int n, int k, double *X, int *cluster_assignment_index, double *cluster_centroid)
{
	int cluster_member_count[MAX_CLUSTERS];

	// get_cluster_member_count(n, k, cluster_assignment_index, cluster_member_count);
	int ii, j;
	// printf("  Final clusters %d\n", k);
	for (ii = 0; ii < k; ii++)
	{
		// printf("    cluster %d:     members: %8d, centroid ", ii, cluster_member_count[ii]);
		for (j = 0; j < dim; j++)
			printf("%lf  ", cluster_centroid[ii * dim + j]);

		printf("\n");
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
	int N ,
	int n,	   // number of elements

	int k,						   // number of clusters
	double *cluster_centroid,	   // initial cluster centroids
	int *cluster_assignment_final, // output
	int *batch_iteration,
	int *change_count)
{
	double T_1_iteration , T_all_iteration ;
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
	
	while ((*batch_iteration) < MAX_ITERATIONS)
	{
		gettimeofday(&km_start, NULL); 
		
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
			printf("  negative progress made on this step - iteration completed (%.2f) \n", totD - prev_totD);
			// done with this phase
			break;
		}

		// save previous step
		copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);
		// move all points to nearest cluster
		calc_all_distances(dim, n, k, X, cluster_centroid, dist);
		//reassign points
		choose_all_clusters_from_distances(dim, n, k, dist, cluster_assignment_cur);
		//calculate change count
		int change_count = assignment_change_count(n, cluster_assignment_cur, cluster_assignment_prev);

		// done with this phase if nothing has changed
		if (change_count == 0){break;}

		prev_totD = totD;
		(*batch_iteration)++;

		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////

		gettimeofday(&km_end, NULL);
		kmeans_iterations_durations[(chunk_ind)*MAX_ITERATIONS + (*batch_iteration)-1] = calc_delai_time(km_start,km_end);

		// CHECKPOINT : at the first kmeans iteration of the first chunk
		if((*batch_iteration) == 1){ 
			if(chunk_ind == 0){
				// save iteration  time at etimation structure
				real_time.km_1_iteration_time = kmeans_iterations_durations[(chunk_ind)*MAX_ITERATIONS + (*batch_iteration)-1] ;
				// calculate the observation times at each opration of kmlio
				calc_cycles_per_elem(N,n) ;

				//get optimal frequency and skip chunk for the first chunk execution
				
				long freq ;
				decide_skip_chunk(N,n,&freq) ;
				set_frequency(freq) ;
				printf("in km , chunk %d : frequency set at %ld , skip chunk is %d \n",chunk_ind,freq,skp_chk) ;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////
	}

	kmeans_nb_iteration[chunk_ind] = (*batch_iteration);

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
	int *marks = (int *)malloc(sizeof(int) * (N / chunk_size)+1);
	char line[100000];
	FILE *src = fopen(source, "r+");
	char delim[3] = "\t";
	char *token;
	int j = 0, mark_index = 1;
	marks[0] = ftell(src);
	// printf ("marks[0] = %d\n", marks[0]);
	while (!feof(src))
	{
		fgets(line, 100000, src);
		// printf("j = %d , ",j) ;
		j = (j + 1) % chunk_size;
		// if(strlen(line)>max_line_size) max_line_size = strlen(line) ;
		if (j == 0)
		{
			marks[mark_index] = ftell(src);
			printf ("marks[%d] = %ld , %d\n", mark_index, ftell(src),j);
			mark_index++;
		}
	}
	// printf("max line size = %ld\n",max_line_size) ;

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

double *getmatrix_bin(char *source, size_t dim, size_t N, int sub, int offset)
{
	FILE *src = fopen(source, "rb");
	printf("binary file opened\n");
	double *X = NULL;
	size_t i = 0, j = 0;
	int stop = 0;

	X = (double *)malloc(sizeof(double) * sub * dim);
	fseek(src, sizeof(double) * offset, SEEK_CUR);
	printf("binary file fseek done\n");

	j = 0;
	while (!feof(src) && !stop)
	{
		fread(&X[j], sizeof(double), 1, src);
		j++;
		if (sub != 0)
			if ((j / (dim)) == sub)
				stop = 1;
	}
	fclose(src);

	return X;
}

double *getmatrix_mmap(char *source, size_t dim, int k, size_t N, int sub, int offset, int *marks)
{

	int src = open(source, O_RDONLY);

	double *X = NULL;
	size_t i = 0, j = 0;
	int stop = 0;

	// decide best buffer size
	int MAX_DOUBLE_LEN = 24;
	int MAX_LINE_SIZE = (MAX_DOUBLE_LEN + 1) * dim;
	char line[MAX_LINE_SIZE];

	// calculate the rest size allowed for the chunk kmeans
	int kmeans_size = (sizeof(double) * k + 2 * sizeof(int) + sizeof(int)) * sub;
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
	// sleep(20) ;
	// begin reading the file
	lseek(src, marks[offset / sub], SEEK_CUR);
	char delim[3] = "\t";
	// var to handle strtok tokens
	char *token;
	// store the buffered lines read form files at once
	char *mapped;

	// load the chunk
	j = 0;
	int L = 0;
	int to_read_lines = sub;
	int s;

	while (!stop)
	{
		// allocate buffer size
		if (L_opt > to_read_lines)
			L_opt = to_read_lines;

		// fill the buffer as much as possible
		/* Memory-map the file. */
		mapped = mmap(0, chunk_text_size, PROT_READ, MAP_PRIVATE, src, 0);

		// check (mapped == MAP_FAILED, "mmap failed\n");
		printf("stop after buf malloc @ = %p\n", mapped);
		sleep(10);

		to_read_lines -= L;

		// convert data strings to matrix values
		for (int l = 0; l < sub; l++)
		{
			for (int d = 0; d < 10; d++)
			{
				X[j] = atof("0");
				j++;
			}
		}

		// if we have read all the chunk lines
		if (sub != 0)
		{
			if ((j / (dim)) == sub)
			{
				stop = 1;
				munmap(mapped, chunk_text_size);
				s = malloc_trim(0);
				printf("stop after buf free @ = %p\n", mapped);
				sleep(10);
			}
		}
	}

	close(src);

	s = malloc_trim(0);
	printf("stop after buffer free , trimming status = %d\n", s);
	sleep(10);

	printf("GET MATRIX %d COMPLETE\n", offset / sub);

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
	// printf("stop after X malloc\n") ;

	// sleep(20) ;
	// begin reading the file
	fseek(src, marks[offset / sub], SEEK_CUR);
	char delim[3] = "\t";
	// var to handle strtok tokens
	char *token;
	// store the buffered lines read form files at once
	char **buf = (char **)calloc(L_opt, sizeof(char *));
	// printf("stop after buf malloc\n") ;

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

		// printf("stop after read\n") ;
		// sleep(15) ;

		// printf("Lmax = %d \tLopt = %d \tL = %d \tread lines = %ld \tmax line size = %ld \tstrings real size %f MB + buffer pointers size %f MB\n",L_max,L_opt, L,(j/(dim)),max_line_size+1+sizeof(char *),buf_len/pow(1024,2),(L_opt*sizeof(char*))/pow(1024,2)) ;

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
		// printf("stop after fill\n") ;
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

	// printf("stop after free\n") ;
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
	// int first = rand()%N ;
	int first = chunk_ind;

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
			//	printf ("%lf\t", X[centers_int[i]*dim+j]);
		}
		// printf ("\n");
	}
	free(centers_int);
	centers_int = NULL;
	free(distance_cur_center);
	distance_cur_center = NULL;

	return centers;
}

void dispersion(double *X, size_t dim, size_t N)
{
	int i, j;
	double *min = malloc(sizeof(size_t) * dim);
	double *max = malloc(sizeof(size_t) * dim);

	for (i = 0; i < dim; i++)
	{
		min[i] = X[i];
	}
	for (i = 0; i < dim; i++)
	{
		max[i] = X[i];
	}

	// iterate over point from the second to find the max
	double d = 0;
	for (i = dim; i < N; i++)
	{
		if (calc_distance(dim, min, X + i * dim) > d)
		{
			for (j = i; j < i + dim; j++)
			{
				max[j - i] = X[j];
			}
			d = calc_distance(dim, min, X + i * dim);
		}
	}
	/*for (i=0; i<dim; i++){
	printf ("%lf\t", max[i]);
	}
	printf("\n");
	for (i=0; i<dim; i++){
	printf ("%lf\t", min[i]);
	}*/
	printf("distance min-max %lf\n", calc_distance(dim, min, max));
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
	// printf("\n");

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

	printf("stop %d\n", res);
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
	int j, index_groupe, l, size = 0, i, found, m, *int_chunk;
	int tmp;
	double *Y;
	float c;
	char **buf;

	index_chunk_elem index[(( N / taille ) -skp_chk)];
	index_chunk_elem *cur[(( N / taille ) -skp_chk)];
	index_chunk_elem *elem;
	for (i = 0; i < (N / taille) - skp_chk ; i++)
	{
		index[i].element = -1;
		index[i].next = NULL;
		cur[i] = NULL;
	}
	// form groups index
	form_index(cluster_assignment, nb_groupes, grp, N);
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

				// tmp = rand()%grp[i].nb_members;
				tmp = (j + i * nb_samples) % grp[i].nb_members;

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
	size = 0;
	// iterate over chunk points
	for (i = 0; i < (N / taille) - skp_chk; i++)
	{
		if (index[i].element != -1)
		{
			// Y = getmatrix(source, dim, taille,taille,i*taille, marks);
			Y = getmatrix_buf(source, dim, k, taille, taille, i * taille, marks);
			// printf ("chunk %d\n", i);
			for (j = 0; j < dim; j++)
			{
				chunk[size * dim + j] = Y[index[i].element * dim + j];
				// printf("%lf\n", Y[index[i].element*dim+j]);
			}
			size++;
			// goto next point in index
			cur[i] = index[i].next;
			while (cur[i] != NULL)
			{
				// printf("ELEMENT %ld\n", cur[i]->element);
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
	// printf ("size %d\n", size);
	return chunk;
}


void compare_solution ( double * result_centroid, char * solution_centers_file_name , /*,char * solution_assignment_file_name*/ size_t dim, size_t k, size_t n){
	char *token; 
	char tmp[1000000]; 
	char delim[3] = "\t"; 
	int index=0, index2;  
	double delta=0; 
	
	FILE *c_solution = fopen (solution_centers_file_name, "r"); 
	double * solution_centroid = (double *) malloc (sizeof(double) * k * dim);
	printf("Real clusters %ld\n", k);
	while(!feof(c_solution)){
		fgets (tmp, 1000000, c_solution); 
		token = strtok(tmp, delim); 
		while (token != NULL){
			solution_centroid[index] = atof(token); 
			printf("%f\t",solution_centroid[index]) ;
			index++;
			token = strtok(NULL, delim); 	
		}
		printf("\n") ;
	}
	fclose(c_solution); 

	int ii, j;
	int min_c ;
	printf("Real clusters %ld\n", k);

	// FILE *a_solution = fopen (solution_assignment_file_name, "r"); 
	// int * assignment_solution = (int *) malloc (sizeof(int) * n);
	// index=0 ;
	// while(!feof(a_solution)){
	// 	fgets (tmp, 1000000, a_solution); 
	// 	token = strtok(tmp, delim); 
	// 	while (token != NULL){
	// 		assignment_solution[index] = atoi(token); 
	// 		printf("%f\t",assignment_solution[index]) ;
	// 		index++;
	// 		token = strtok(NULL, delim); 	
	// 	}
	// 	printf("\n") ;
	// }
	// fclose(a_solution); 
	// int * result_assignement = malloc(sizeof(int) * n) ;

	double * dist= malloc(sizeof(double)*k*k) ;

	calc_all_distances(dim,k,k,result_centroid,solution_centroid,dist) ;
	
	/////////////////////////////////////////

	// // Load the R function (source function defined above)
	// source("generator/cluster_evaluation.R") ;
 
	// // Allocate an R vector and copy the C array into it.
	// SEXP arg;
	// PROTECT(arg = allocVector(INTSXP, k*k));
	// memcpy(INTEGER(arg), dist, k * k * sizeof(int));

	free(dist);







	////////////////////////



	// for (index = 0; index<k; index++){
	// 	printf ("centres %d\n", index); 
	// 	for (index2=0; index2<dim; index2++){
	// 		printf ("sol: %lf  sol obt : %lf\n", solution_centroid[index*dim+index2], centroid[index*dim+index2]); 
	// 	}
		
	// 	delta +=calc_distance(dim,&solution_centroid[index*dim], &centroid[index*dim]);  	
	// }
	// printf ("delta solution_type solution obtenue : %lf\n", sqrt(delta)); 

	free(solution_centroid); 
}


void kmeans_by_chunk(char *source, size_t dim, int taille, int N, int k,double ** cluster_centroid)
{
	int i, j, n, m = 0, kmeans_iterations = 0, kmeans_change_count = 0;
	double *Y = NULL;
	int *cluster_assignment_final = NULL, *cluster_assignment_final_Y = NULL;
	double *cluster_centroid_by_chunk = NULL, *chunk_centroid = NULL;
	double *var = NULL;
	int cluster_member_count[MAX_CLUSTERS]; // à modifier pour une allocation dynamique
	int nb_groupes = 0;
	groupe *groupes = NULL;

	// //////////////////////////////////////////
	set_frequency(MIN_FREQ) ;
	struct timeval tv_1, tv_2 , chunk_start, chunk_end ;
	long freq ;
	
	// //////////////////////////////////////////
	
	// *cluster_centroid = NULL ;

	/****** PHASE 1 : PARTIELS CHUNKS KMEANS ******/
	gettimeofday(&tv_1, NULL);
	// mark the chunks in dataset
	int *marks = mark(source, dim, N, taille);
	gettimeofday(&tv_2, NULL);
	real_time.mark_time = calc_delai_time(tv_1,tv_2) ;
	
	set_frequency(BASE_FREQ) ;

	// apply kmeans on each chunk
	if ((( N / taille ) - skp_chk) > 1)
	{
		cluster_assignment_final = (int *)malloc(N * sizeof(int));
		cluster_centroid_by_chunk = (double *)malloc( sizeof(double) * (k * dim * (( N / taille ) -skp_chk) ) );
		var = (double *)malloc(sizeof(double) * k * (( N / taille ) -skp_chk) );
		
		//////////////////////////////////////
		kmeans_iterations_durations = malloc( sizeof(double) * MAX_ITERATIONS * (( N / taille )+1) ) ;
		kmeans_nb_iteration = malloc( sizeof(int) * (( N / taille )+1) ) ;
		chunks_elapsed_durations = malloc( sizeof(double) * (( N / taille )+1) ) ;
		//////////////////////////////////////
		m=0 ;
		while(m<(N / taille) - skp_chk)
		{
			gettimeofday(&chunk_start, NULL);
			//////////////////////////////////////

			// lecture du chunk
			i = m * taille;
			Y = getmatrix_buf(source, dim, k, taille, taille, i, marks);

			if(m==0){
				gettimeofday(&tv_1, NULL);
				real_time.getmat_time = calc_delai_time(chunk_start,tv_1) ;
			}

			//////////////////////////////////////
		
			// chunk m init kmeans++
			(*cluster_centroid) = kmeans_init_plusplus(Y, taille, dim, k);
			if(m==0){
				gettimeofday(&tv_2, NULL);
				real_time.init_time = calc_delai_time(tv_1,tv_2) ;
			}

			//////////////////////////////////////

			// apply kmeans on the chunk m
			cluster_assignment_final_Y = (int *)malloc(taille * sizeof(int));
			kmeans(dim, Y,N,taille, k, (*cluster_centroid), cluster_assignment_final_Y, &kmeans_iterations, &kmeans_change_count);
			
			//////////////////////////////////////

			//get optimal frequency and skip chunk for the first chunk execution
				
			printf("chunk id %d completed\n",m) ;
			chunk_ind ++ ;
			decide_skip_chunk(N,taille,&freq) ;

			// ////////////////////////////////////
			
			//if the first chunk is  the sole one to be treated , end kmlio with its centroids results
			if(m==0){
				if( ( ( N / taille ) - m - skp_chk ) <= 1 ){
					printf("stopped after the first chunk excecution 1/%d\n",N/taille) ;
					free(Y);
					Y = NULL;
					gettimeofday(&chunk_end, NULL);
					chunks_elapsed_durations[m]=calc_delai_time(chunk_start,chunk_end) ;
					return ;
				}
			}

			//////////////////////////////////////
			
			gettimeofday(&tv_1, NULL);
			
			// variance calculation on the chunk m
			for (j = 0; j < k; j++)
			{
				var[m * k + j] = var_calculate(Y, dim, (*cluster_centroid), j, cluster_assignment_final_Y, taille) / taille;
			}
			// copy the centroids to global array
			for (j = 0; j < k * dim; j++)
			{
				cluster_centroid_by_chunk[m * k * dim + j] = (*cluster_centroid)[j];
			}
			// copy points clustering to global array
			for (n = 0; n < taille; n++)
			{
				cluster_assignment_final[m * taille + n] = cluster_assignment_final_Y[n] + m * k;
			}

			free(Y);
			Y = NULL;

			free((*cluster_centroid));
			(*cluster_centroid) = NULL;

			free(cluster_assignment_final_Y);
			cluster_assignment_final_Y = NULL;

			gettimeofday(&chunk_end, NULL);
			if(m==0)
				real_time.var_copy_time = calc_delai_time(tv_1,chunk_end) ;

			////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////

			chunks_elapsed_durations[m]=calc_delai_time(chunk_start,chunk_end) ;

			if( ( ( N / taille ) - m - skp_chk ) <= 1 ){
				printf("skipped chunk reached after %d of %d\n",m+1,N/taille) ;
				break ;
			}else{
				km_freq[m] = freq ;
				set_frequency(freq) ;
				printf("chunk %d estimation : frequency set at %ld , skip chunk is %d \n",chunk_ind,freq,skp_chk) ;
			}

			////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////

			m++ ;
		}

		// decide_skip_chunk(N,taille,&freq) ;
		int found  = 0 ;
		freq = get_optimal_freq(N,taille,skp_chk,&found) ;
		set_frequency(freq) ;
		km_freq[chunk_ind] = freq ;
		
		printf("final chunk n %d : frequency set at %ld , skip chunk is %d \n",chunk_ind,freq,skp_chk) ;

		/****** PHASE 2 : PARTIELS CLUSTERS GROUPING ******/
		
		groupes = (groupe *)malloc(sizeof(groupe) * (( N / taille ) - skp_chk) * k) ;
		nb_groupes = 0 ;
		for (j = 0; j < k * (N / taille) - skp_chk ; j++)
		{
			if (!existG(groupes, nb_groupes, j))
			{
				groupes[nb_groupes].means = (int *)malloc(sizeof(int));
				groupes[nb_groupes].nb = 1;
				groupes[nb_groupes].means[0] = j;
				for (n = j + 1; n < k * ((N / taille) - skp_chk ) ; n++)
				{
					// condition de non existence
					if (calc_distance(dim, &cluster_centroid_by_chunk[j * dim], &cluster_centroid_by_chunk[n * dim]) < var[j] && !existG(groupes, nb_groupes, n))
					{
						groupes[nb_groupes].means = (int *)realloc(groupes[nb_groupes].means, (groupes[nb_groupes].nb + 1) * sizeof(int));
						groupes[nb_groupes].means[groupes[nb_groupes].nb] = n;
						// printf ("next\n");
						groupes[nb_groupes].nb++;
					}
				}
				nb_groupes++;
			}
		}

		free(cluster_centroid_by_chunk);
		cluster_centroid_by_chunk = NULL;

		// groups members count update
		float count = k *(( N / taille ) - skp_chk);
		get_cluster_member_count(((int)(( N / taille ) - skp_chk)) * taille, (int)count, cluster_assignment_final, cluster_member_count);

		for (j = 0; j < nb_groupes; j++)
		{
			groupes[j].nb_members = 0;
			for (n = 0; n < groupes[j].nb; n++)
			{
				groupes[j].nb_members += cluster_member_count[groupes[j].means[n]];
			}
		}

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

		(*cluster_centroid) = kmeans_init_plusplus(chunk, taille, dim, k);

		cluster_assignment_final_Y = (int *)malloc(taille * sizeof(int));
		kmeans(dim, chunk,N , taille, k, (*cluster_centroid), cluster_assignment_final_Y, &kmeans_iterations, &kmeans_change_count);
		
		free(chunk);
		chunk = NULL;
	}
	else
	{
		double *X;

		X = getmatrix_buf(source, dim, k, N, N, 0, marks);
	
		(*cluster_centroid) =kmeans_init_plusplus(X, N, dim, k);

		cluster_assignment_final=  (int *) malloc(N*sizeof(int));
		kmeans(dim,X,N,N, k, (*cluster_centroid), cluster_assignment_final,&kmeans_iterations,&kmeans_change_count);

		free(X);
		X = NULL;
	}

	free(cluster_assignment_final_Y);
	cluster_assignment_final_Y = NULL;

	// free(cluster_centroid);
	// cluster_centroid = NULL;

	free(marks);
	marks = NULL;
}


int main(int argc, char **argv)
{

	cpu_set_t set;
	// clear cpu mask
	CPU_ZERO(&set);
	// set cpu 0
	CPU_SET(0, &set);
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
	sprintf (cmd, "echo %d > /sys/fs/cgroup/memory/kmeans/cgroup.procs", getpid()); // /cgroups/mem/kmeans/tasks
	printf ("%s\n", cmd);
	system(cmd);

	sprintf (cmd, "echo %d > /sys/fs/cgroup/blkio/kmeans/cgroup.procs", getpid()); // /cgroups/mem/kmeans/tasks
	printf ("%s\n", cmd);
	system(cmd);

	sleep(1);

	char *source = argv[1];
	size_t k = atoi(argv[2]);
	N = atoi(argv[3]);
	taille = atoi(argv[4]);
	dim = atoi(argv[5]);
	D_max = atoi(argv[6]);	// kmlio maximal delay time in seconds
	skp_chk = 0 ;
	
	srand(time(NULL));

	//begin = clock();
	double * centroid = NULL ;

	km_freq = malloc(sizeof(long)*(N/taille+1)) ;
	
	gettimeofday(&kmlio_start, NULL); 
	
	kmeans_by_chunk(source,dim, taille, N,k,&centroid);
	
	gettimeofday(&kmlio_end, NULL); 

	kmlio_diag(k,dim,N,taille,D_max,centroid) ;

	return 0;
}
