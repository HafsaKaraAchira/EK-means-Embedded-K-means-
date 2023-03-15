// kmeans.c
// Ethan Brodsky
// modified by Camélia SLIMANI

/*

DELL : cd /home/hafsa/Documents/@K-MLIO_Analysis/
BeagleBone : cd /home/debian/@K-MLIO_Analysis/
scripts/prog_script_cgroup 1104
scripts/prog_script_reset

(scripts/prog_script_reset) && (clear && gcc -g kmlio.c -o program -lm -D _GNU_SOURCE) && (./program generator/CM13,4M_2400MO_SEP0,2/points.csv 10 13421800 6710900 10)
(scripts/prog_script_reset) && (clear && gcc -g kmlio.c -o program -lm -D _GNU_SOURCE) && ((./program generator/CM13,4M_2400MO_SEP0,2/points.csv 10 13421800 6710900 10) & (taskset -c 1 scripts/prog_script_launch))

(scripts/prog_script_cache) && (clear && gcc -g perf_stats/getmatrix_perf.c -o perf_stats/getmatrix_perf -lm -D _GNU_SOURCE)
./perf_stats/init_perf generator/CM13,4M_2400MO_SEP0,2/points.csv 2 13421800 6710900 10

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
	char **buf = (char **)calloc(sub, sizeof(char *));
	char line[] ="0.500040325395034\t0.491522145241467\t0.469243971882188\t0.653238734557954\t0.319925902731301\t0.544123893967396\t0.526092562503245\t0.467128339965845\t0.560311065512695\t0.575082213530508\n" ;
	// printf("stop after buf malloc\n") ;

	// load the chunk
	j = 0;
	int L = 0;
	int to_read_lines = sub;
	int s;
	size_t max_len;
	int cpid , pid ;

	// if(i==1){
		pid= getpid();
		cpid = fork();
	// }
		
	if( cpid == 0)
	{
		// child process .  Run your perf stat
		char buf[50];
		sprintf(buf, "perf stat -p %d   >> init_stat.log 2>&1",pid);
		execl("/bin/sh", "sh", "-c", buf, NULL);
	}
	else
	{
		// set the child the leader of its process group
		setpgid(cpid, 0);

		//////////////////////////////////////////////
		// part of program you wanted to perf stat ...
		
			while (!stop)
		{
			// allocate buffer size
			// if (L_opt > to_read_lines)
			// 	L_opt = to_read_lines;
			// fill the buffer as much as possible
			L = 0;
			size_t buf_len = 0;

			while (!feof(src) && L < sub)
			{
				// fgets(line, sizeof(line), src);
				
				buf[L] = calloc(strlen(line), sizeof(char));
				strncpy(buf[L], line, strlen(line));
				buf[L][strlen(line) - 1] = '\0';
				L++;
			}

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

			// if we have read all the chunk lines
			if (sub != 0)
			{
				if ((j / (dim)) == sub)
				{
					stop = 1;
				}
			}
		}
		////////////////////////////////////////////////////////////////
		// stop perf stat by killing child process and all its descendants(sh, perf stat etc )
		kill(-cpid, SIGINT);
		////////////////////////////////////////////////////////////////////

		// rest of the program ...
		if (buf != NULL)
		{
			free(buf);
			buf = NULL;
		}
		
	}

	

	fclose(src);

	// s = malloc_trim(0);


	return X;
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
	//set priority
	setpriority(PRIO_PROCESS, 0, -20);

	// 2401015202

	char *source = argv[1];
	size_t k = atoi(argv[2]);
	size_t N = atoi(argv[3]);
	int M = atoi(argv[4]);
	size_t dim = atoi(argv[5]);

	double *Y = NULL;
	double *cluster_centroid = NULL;

	// mark the chunks in dataset
	int *marks = mark(source, dim, N, M);
	int i = 0 * M;
	// read the chunk
	Y = getmatrix_buf(source, dim, k, M, M, i, marks);

	// chunk m init kmeans++
	// cluster_centroid = kmeans_init_plusplus(Y, M, dim, k);

	free(Y);
	Y = NULL;

	free(cluster_centroid);
	cluster_centroid = NULL;

	// kmeans_by_chunk( source, dim, M, N, atoi(argv[2]));



	// int *marks = (int *)malloc(sizeof(int) * 2);
	// marks[0] = 0;
	// marks[1] = 1200507295;
	// marks[2] = 2401015202;

	// X = getmatrix_buf(source, dim, k, M, M, 0, marks);
	// getmatrix_buf(source, dim,k, N,N,0, marks);
	// char *buf = malloc(sizeof(char)*strlen(line)) ;

	
	// char line[] ="0.500040325395034\t0.491522145241467\t0.469243971882188\t0.653238734557954\t0.319925902731301\t0.544123893967396\t0.526092562503245\t0.467128339965845\t0.560311065512695\t0.575082213530508\n" ;
	
	// double *X=malloc(sizeof(double)*dim*M);
	// // store the buffered lines read form files at once
	// char **buf = (char **)calloc(M, sizeof(char *));
	// int L = 0;
	// size_t buf_len = 0;
	// while (L < M)
	// {
	// 	buf[L] = calloc(strlen(line), sizeof(char));
	// 	strncpy(buf[L], line, strlen(line));
	// 	buf[L][strlen(line) - 1] = '\0';
	// 	L++;
	// }

	// for (int l = 0; l <M; l++)
	// {
	// 	int j = 0 ;
	// 	char * token = strtok(buf[l],"\t");
	// 	while (token != NULL)
	// 	{
	// 		X[j] = atof(token);
	// 		j++;
	// 		token = strtok(NULL,"\t");
	// 	}
	// 	free(buf[l]);
	// 	buf[l] = NULL;
	// }
	// free(buf);
	// buf = NULL;
	// free(X);
	// X = NULL;

	return 0;
}
