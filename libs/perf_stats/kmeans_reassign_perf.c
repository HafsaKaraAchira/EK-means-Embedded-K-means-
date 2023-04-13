// kmeans.c
// Ethan Brodsky
// modified by Camélia SLIMANI

/*

DELL : cd /home/hafsa/Documents/@K-MLIO_Analysis/
BeagleBone : cd /home/debian/@K-MLIO_Analysis/
scripts/prog_script_cgroup 1104
scripts/prog_script_reset

(scripts/prog_script_reset) && (clear && gcc -g kmlio.c -o program -lm -D _GNU_SOURCE) && (./progra generator/CM13,4M_2400MO_SEP0,2/points.csv 10 13421800 6710900 10)
(scripts/prog_script_reset) && (clear && gcc -g kmlio.c -o program -lm -D _GNU_SOURCE) && ((./program generator/CM13,4M_2400MO_SEP0,2/points.csv 10 13421800 6710900 10) & (taskset -c 1 scripts/prog_script_launch))

(scripts/prog_script_cache) && (clear && gcc -g perf_stats/kmeans_reassign_perf.c -o perf_stats/kmeans_reassign_perf -lm -D _GNU_SOURCE)
./perf_stats/kmeans_reassign_perf generator/CM13,4M_2400MO_SEP0,2/points.csv 2 13421800 6710900 10

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

#include <signal.h>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <sched.h>

#define sqr(x) ((x) * (x))

#define MAX_ITERATIONS 150

#define BIG_double (INFINITY)

void fail(char *str)
{
	printf("%s", str);
	exit(-1);
}

int *mark(char *source, size_t dim, size_t N, int chunk_size)
{
	int *marks = (int *)malloc(sizeof(int) * (N / chunk_size));
	char line[100000];
	FILE *src = fopen(source, "r+");
	char delim[3] = "\t";
	char *token;
	int j = 0, mark_index = 1;
	marks[0] = ftell(src);
	while (!feof(src))
	{
		fgets(line, 100000, src);
		j = (j + 1) % chunk_size;
		if (j == 0)
		{
			marks[mark_index] = ftell(src);
			mark_index++;
		}
	}
	// printf("max line size = %ld\n",max_line_size) ;

	fclose(src);
	return marks;
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

	return X;
}

double calc_distance(int dim, double *p1, double *p2)
{
	double distance_sq_sum = 0;
	int ii;
	for (ii = 0; ii < dim; ii++)
		distance_sq_sum += sqr(p1[ii] - p2[ii]);

	return distance_sq_sum;
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

double *kmeans_init_plusplus(double *X, size_t N, size_t dim, size_t k)
{
	double *centers = (double *)malloc(sizeof(double) * k * dim);
	double *distance_cur_center = (double *)malloc(sizeof(double) * N);
	int *centers_int = (int *)malloc(sizeof(int) * k);
	double sum = 0;
	// int first = rand() % N;
	int first = 0; // chunk_ind++ ;

	int i, j, best;
	centers_int[0] = first;
	// centers_int[1] = 4386979;
	// centers_int[2] = 4801150;
	// centers_int[3] = 4811196;
	// centers_int[4] = 996475;
	// centers_int[5] = 2592191;
	// centers_int[6] = 6542895;
	// centers_int[7] = 6311008;
	// centers_int[8] = 4662352;
	// centers_int[9] = 3313892;

	////////////////////////////////////////////////

	for (i = 1; i < k; i++)
	{
		for (j = 0; j < N; j++)
		{
			distance_cur_center[j] = affect_x(X, centers_int, j, i, dim);
			sum += distance_cur_center[j];
		}
		centers_int[i] = find_best_center(distance_cur_center, centers_int, N, dim, i, sum);
		// printf("i = %d\n",centers_int[i]) ;
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

void calc_cluster_centroids(int dim, int n, int k, double *X, int *cluster_assignment_index, double *new_cluster_centroid)
{
	int cluster_member_count[k];
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
	int n,	   // number of elements

	int k,						   // number of clusters
	double *cluster_centroid,	   // initial cluster centroids
	int *cluster_assignment_final, // output
	int *batch_iteration,
	int *change_count)
{
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

	int cpid, pid;

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
		// done with this phase
	}

	// save previous step
	copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);

	calc_all_distances(dim, n, k, X, cluster_centroid, dist);

	// save_kmeans_iterations(-1,-1,chunk_ind,-1,-1,-1);

	// if((*batch_iteration)== 0){
	pid = getpid();
	cpid = fork();
	// }

	if (cpid == 0)
	{
		// child process .  Run your perf stat
		char msg[50];
		sprintf(msg, "perf stat -p %d   >> kmeans_stat.log 2>&1", pid);
		execl("/bin/sh", "sh", "-c", msg, NULL);
	}
	else
	{
		// set the child the leader of its process group
		setpgid(cpid, 0);

		while ((*batch_iteration) < MAX_ITERATIONS)
		{
			// update cluster centroids
			// calc_cluster_centroids(dim, n, k, X, cluster_assignment_cur, cluster_centroid);

			// see if we've failed to improve
			// double totD = calc_total_distance(dim, n, k, X, cluster_centroid, cluster_assignment_cur);
			
			// if (totD > prev_totD)
			// // failed to improve - currently solution worse than previous
			// {
			// 	printf("kmeans restore\n") ;
			// 	// restore old assignments
			// 	copy_assignment_array(n, cluster_assignment_prev, cluster_assignment_cur);
			// 	// recalc centroids
			// 	calc_cluster_centroids(dim, n, k, X, cluster_assignment_cur, cluster_centroid);
			// 	// done with this phase
			// 	break;
			// }

			// save previous step
			copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);

			// move all points to nearest cluster
			// calc_all_distances(dim, n, k, X, cluster_centroid, dist);
			
			//REASSIGN POINTS
			// choose_all_clusters_from_distances(dim, n, k, dist, cluster_assignment_cur);

			int change_count = assignment_change_count(n, cluster_assignment_cur, cluster_assignment_prev);

			// done with this phase if nothing has changed
			if (change_count == 0){	}//break;}

			prev_totD = totD;

			//  save_kmeans_iterations(-1,-1,chunk_ind,(*batch_iteration),change_count,totD) ;
			(*batch_iteration)++;
		}
		printf("kmeans iterations%d\n", (*batch_iteration));
		////////////////////////////////////////////////////////////////
		// stop perf stat by killing child process and all its descendants(sh, perf stat etc )
		kill(-cpid, SIGINT);
		////////////////////////////////////////////////////////////////////
		// rest of the program ...
	}

	// write to output array
	copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_final);
	free(dist);
	dist = NULL;
	free(cluster_assignment_cur);
	cluster_assignment_cur = NULL;
	free(cluster_assignment_prev);
	cluster_assignment_prev = NULL;
}

void kmeans_by_chunk(char *source, size_t dim, int taille, int N, int k)
{
	int i, j, n, m = 0, kmeans_iterations = 0, kmeans_change_count = 0;
	double *Y = NULL;
	int *cluster_assignment_final = NULL, *cluster_assignment_final_Y = NULL;
	double *cluster_centroid = NULL, *cluster_centroid_by_chunk = NULL, *chunk_centroid = NULL;
	double *var = NULL;
	// int cluster_member_count[k*taille]; //à modifier pour une allocation dynamique
	// int nb_groupes=0;
	// groupe * groupes=NULL;

	/****** PHASE 1 : PARTIELS CHUNKS KMEANS ******/

	// mark the chunks in dataset
	// save_phase_time(1,1,0,-1);
	int *marks = mark(source, dim, N, taille);
	// save_phase_time(1,1,1,-1);

	// apply kmeans on each chunk
	if (N / taille > 1)
	{
		// save_phase_time(1,2,0,-1);
		// cluster_assignment_final = (int *) malloc(N*sizeof(int));
		// cluster_centroid_by_chunk = (double *)malloc(sizeof(double)*(k*dim*N/taille));
		// var = (double *) malloc (sizeof(double)*k*N/taille);

		for (m = 0; m < N / taille; m++)
		{
			// lecture du chunk
			i = m * taille;
			Y = getmatrix_buf(source, dim, k, taille, taille, i, marks);
			// save_phase_time(1,2,(m+1),-1);	//get matrix complete

			// chunk m init kmeans++
			cluster_centroid = kmeans_init_plusplus(Y, taille, dim, k);
			// save_phase_time(1,3,(m+1),-1);

			// apply kmeans on the chunk m
			//  cluster_assignment_final_Y=  (int *) malloc(taille*sizeof(int));
			//  kmeans(dim,Y,taille, k, cluster_centroid, cluster_assignment_final_Y,&kmeans_iterations,&kmeans_change_count);
			//  save_kmeans_iterations(1,4,(m+1),kmeans_iterations,kmeans_change_count,-1);

			// variance calculation on the chunk m
			//  for ( j=0; j<k;j++){
			//  	var[m*k+j] = var_calculate (Y, dim, cluster_centroid, j, cluster_assignment_final_Y , taille)/taille;
			//  }
			//  //copy the centroids to global array
			//  for ( j=0; j<k*dim;j++){
			//  	cluster_centroid_by_chunk[m*k*dim+j] = cluster_centroid[j];
			//  }
			//  //copy points clustering to global array
			//  for ( n=0 ; n<taille; n++){
			//  	cluster_assignment_final[m*taille+n]=cluster_assignment_final_Y[n]+m*k;
			//  }

			free(Y);
			Y = NULL;

			free(cluster_centroid);
			cluster_centroid = NULL;

			// free(cluster_assignment_final_Y);
			// cluster_assignment_final_Y=NULL;
		}

		/****** PHASE 2 : PARTIELS CLUSTERS GROUPING ******/
		// save_phase_time(2,1,0,-1);
		// groupes = (groupe * ) malloc (sizeof(groupe )*N/taille*k);
		// nb_groupes =0;
		// for (j=0; j< k*N/taille; j++){
		// 	if (!existG(groupes, nb_groupes, j)){
		// 		groupes[nb_groupes].means = (int*) malloc(sizeof(int));
		// 		groupes[nb_groupes].nb =1;
		// 		groupes[nb_groupes].means[0]=j;
		// 		for ( n = j+1; n<k*N/taille; n++){
		// 			//condition de non existence
		// 			if (calc_distance(dim, &cluster_centroid_by_chunk[j*dim], &cluster_centroid_by_chunk[n*dim])<var[j] && !existG(groupes, nb_groupes, n)){
		// 				groupes[nb_groupes].means = (int *) realloc (groupes[nb_groupes].means, (groupes[nb_groupes].nb+1)*sizeof(int));
		// 				groupes[nb_groupes].means[groupes[nb_groupes].nb] = n;
		// 				//printf ("next\n");
		// 				groupes[nb_groupes].nb++;

		// 			}
		// 		}
		// 		nb_groupes++;
		// 	}
		// }

		// save_phase_time(2,1,1,j);

		// free(cluster_centroid_by_chunk);
		// cluster_centroid_by_chunk = NULL;

		// // save_phase_time(2,2,0,-1);
		// //groups members count update
		// float count = k*N/taille;
		// // get_cluster_member_count(((int) (N/taille)) * taille, (int) count, cluster_assignment_final, cluster_member_count);

		// for (j=0; j<nb_groupes; j++){
		// 	groupes[j].nb_members = 0;
		// 	for (n = 0; n< groupes[j].nb; n++){
		// 		groupes[j].nb_members+=cluster_member_count[groupes[j].means[n]];
		// 	}
		// }

		// save_phase_time(2,2,1,nb_groupes);

		/****** PHASE 3 : FINAL CHUNK BUILDING ******/
		// save_phase_time(3,1,0,-1);
		// double * chunk;
		// chunk = form_chunk(groupes,source, marks,cluster_assignment_final, nb_groupes, N, dim,k, taille);
		// // save_phase_time(3,1,1,-1);

		// free(cluster_assignment_final);
		// cluster_assignment_final =NULL;
		// for (i = 0; i<nb_groupes; i++){
		// 	free(groupes[i].means);
		// 	free(groupes[i].members);
		// }
		// free(groupes);
		// groupes = NULL;

		// /****** PHASE 4 : FINAL CHUNK KMEANS ******/

		// // save_phase_time(4,1,0,-1);
		// cluster_centroid =kmeans_init_plusplus(chunk, taille, dim, k);
		// // save_phase_time(4,1,1,-1);

		// // save_phase_time(4,2,0,-1);
		// cluster_assignment_final_Y=  (int *) malloc(taille*sizeof(int));
		// kmeans(dim,chunk,taille, k, cluster_centroid, cluster_assignment_final_Y,&kmeans_iterations,&kmeans_change_count);
		// // save_kmeans_iterations(4,2,1,kmeans_iterations,kmeans_change_count,-1);

		// // r8mat_write ("results/centers.csv",dim,k,cluster_centroid) ;
		// free(chunk);
		// chunk = NULL;
	}
	else
	{
		// sleep(20) ;
		//  save_phase_time(2,1,0,-1);
		double *X;

		X = getmatrix_buf(source, dim, k, N, N, 0, marks);
		// save_phase_time(2,1,1,-1);

		// save_phase_time(2,2,0,-1);
		// cluster_centroid =kmeans_init_plusplus(X, N, dim, k);
		// save_phase_time(2,2,1,-1);

		// save_phase_time(2,3,0,-1);
		// cluster_assignment_final=  (int *) malloc(N*sizeof(int));
		// kmeans(dim,X,N, k, cluster_centroid, cluster_assignment_final,&kmeans_iterations,&kmeans_change_count);
		// save_kmeans_iterations(2,3,1,kmeans_iterations,kmeans_change_count,-1);

		free(X);
		X = NULL;
	}

	free(cluster_assignment_final_Y);
	cluster_assignment_final_Y = NULL;
	free(cluster_centroid);
	cluster_centroid = NULL;
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
	// 2401015202

	char *source = argv[1];
	size_t k = atoi(argv[2]);
	size_t N = atoi(argv[3]);
	int M = atoi(argv[4]);
	size_t dim = atoi(argv[5]);

	// kmeans_by_chunk( source, dim, M, N, atoi(argv[2]));

	double *Y = NULL;
	double *cluster_centroid = NULL;
	int kmeans_iterations = 0, kmeans_change_count = 0;
	int *cluster_assignment_final_Y = NULL;

	// mark the chunks in dataset
	int *marks = mark(source, dim, N, M);
	int i = 0 * M;
	// read the chunk
	Y = getmatrix_buf(source, dim, k, M, M, i, marks);

	// chunk m init kmeans++
	cluster_centroid = kmeans_init_plusplus(Y, M, dim, k);

	cluster_assignment_final_Y = (int *)malloc(M * sizeof(int));
	kmeans(dim, Y, M, k, cluster_centroid, cluster_assignment_final_Y, &kmeans_iterations, &kmeans_change_count);

	free(Y);
	Y = NULL;

	free(cluster_centroid);
	cluster_centroid = NULL;

	return 0;
}
