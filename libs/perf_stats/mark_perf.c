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

#ifndef _GNU_SOURCE
#define	_GNU_SOURCE 1
#endif

#include <sched.h>

#define sqr(x) ((x)*(x))

#define MAX_ITERATIONS 150

#define BIG_double (INFINITY)

void fail(char *str)
  {
    printf("%s",str);
    exit(-1);
  }

int *mark(char *source, size_t dim, size_t N, int chunk_size){
	int  * marks= (int *) malloc (sizeof(int)*(N/chunk_size));
	char line[100000];
	FILE *src = fopen(source, "r"); 
	char delim[3]="\t"; 
	char *token; 
	int j = 0, mark_index =1 ;
	marks[0] = ftell(src); 
	while (!feof(src)){
		fgets (line,100000, src); 
		j=(j+1)%chunk_size; 
		if (j==0){
			marks[mark_index]=ftell(src);
			mark_index++; 
		}
	}
	fclose(src);
	return marks; 
}

int main (int argc, char **argv){

	// cpu_set_t set;
	// // clear cpu mask
	// CPU_ZERO(&set);
	// // set cpu 0
	// CPU_SET(0, &set);
	// // 0 is the calling process
	// sched_setaffinity(0, sizeof(cpu_set_t), &set);
	// //set priority
	// setpriority(PRIO_PROCESS, 0, -20);

	size_t dim, N; 
	// double *X;

	// 2401015202
	
	char * source = argv[1]; 
	size_t k = atoi(argv[2]);
	N = atoi(argv[3]); 
	int M = atoi(argv[4]);
	dim = atoi(argv[5]);

	// kmeans_by_chunk( source, dim, M, N, atoi(argv[2])); 

	// int *marks = mark(source, dim, N, M);

	int  * marks= (int *) malloc (sizeof(int)*1);
	char line[100000];
	FILE *src = fopen(source, "r"); 
	char delim[3]="\t"; 
	char *token; 
	int j = 0, mark_index =1 ;
	marks[0] = ftell(src); 
	while (!feof(src)){
		fgets (line,100000, src); 
		j=(j+1)%1; 
		if (j==0){
			marks[mark_index]=ftell(src);
			mark_index++; 
			break;
		}
	}
	fclose(src);
		
	return 0;
}
