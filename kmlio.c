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
// #include <gc.h>
#include <malloc.h>
#include <sys/mman.h>

#ifndef _GNU_SOURCE
#define	_GNU_SOURCE 1
#endif

#include <sched.h>


#define sqr(x) ((x)*(x))

#define MAX_CLUSTERS 100000

#define MAX_ITERATIONS 25 // decrease iterations number for beaglebone overheating issue

#define BIG_double (INFINITY)

int chunk_index = 1 ;

size_t max_line_size = 0 ;



struct groupe{
int * means; 
int nb; 
int nb_members;
int * members;
};
typedef struct groupe groupe; 

struct index_chunk_elem{
int  element; 
struct index_chunk_elem *next;
};
typedef struct index_chunk_elem index_chunk_elem; 

void fail(char *str)
  {
    printf("%s",str);
    exit(-1);
  }


void save_phase_time(int phase, int step ,int loop , int iteration){
    char *cmd;
	//printf("iterations=%d",iteration);
    asprintf(&cmd,"taskset -c $(($(nproc) - 1)) ./scripts/prog_script_phase %d %d %d %d",phase,step,loop,iteration);
	system(cmd);
	free(cmd) ;
	cmd=NULL ;
}

double *r8mat_data_read ( char *input_filename, int m, int n )

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
# define MY_LINE_MAX 255

  int error;
  char *got_string;
  FILE *input;
  int i;
  int j;
  char line[255];
  double *table;
  double *x;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
  }

  table = ( double * ) malloc ( m * n * sizeof ( double ) );

  x = ( double * ) malloc ( m * sizeof ( double ) );

  j = 0;

  while ( j < n )
  {
    got_string = fgets ( line, MY_LINE_MAX, input );

    if ( !got_string ){break;}

    // if ( line[0] == '#' || s_len_trim ( line ) == 0 ){continue;}

	char delim[3] = "\t" ;

	char *token = strtok(line,delim);
	i = 0 ;
	while (token != NULL){
		x[i] = atof(token); 
		i++;
		token = strtok(NULL, delim);
	}

    //error = s_to_r8vec ( line, m, x );

    //if ( error == 1 ){continue;}

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  fclose ( input );

  free ( x );

  return table;

# undef MY_LINE_MAX
}
/******************************************************************************/
 
void r8mat_write ( char *output_filename, int m, int n, double table[] )

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
  output = fopen ( output_filename, "at" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the file '%s'.\n", output_filename );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "%2.16f\t", table[i+j*m] ); //%24.16g
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

void i4mat_write ( char *output_filename, int m, int n, int table[] )

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
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file '%s'\n", output_filename );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %d", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

void rbinmat_write ( char *output_filename, int m, int n, double table[] )

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
  output = fopen ( output_filename, "at" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the file '%s'.\n", output_filename );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fwrite ( &table[i+j*m],sizeof(double),1,output ); //%24.16g
    }
  }
/*
  Close the file.
*/
  fclose ( output );

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
    for ( ii = 0; ii < n; ii++) // for each point
      for ( jj = 0; jj < k; jj++) // for each cluster
        {
         // calculate distance between point and cluster centroid
          distance_output[ii*k + jj] = calc_distance(dim, &X[ii*dim], &centroid[jj*dim]);
        }
  }
  
double calc_total_distance(int dim, int n, int k, double *X, double *centroids, int *cluster_assignment_index)
 // NOTE: a point with cluster assignment -1 is ignored
  {
    double tot_D = 0;
    int ii; 
   // for every point
    for ( ii = 0; ii < n; ii++)
      {
       // which cluster is it in?
        int active_cluster = cluster_assignment_index[ii];
        
       // sum distance
        if (active_cluster != -1)
          tot_D += calc_distance(dim, &X[ii*dim], &centroids[active_cluster*dim]);
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
        for ( jj = 0; jj < k; jj++)
          {
           // distance between point and cluster centroid
           
            double cur_distance = distance_array[ii*k + jj];
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
    for ( ii = 0; ii < k; ii++) 
      {
        cluster_member_count[ii] = 0;
        
        for ( jj = 0; jj < dim; jj++)
          new_cluster_centroid[ii*dim + jj] = 0;
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
        for ( jj = 0; jj < dim; jj++)
          new_cluster_centroid[active_cluster*dim + jj] += X[ii*dim + jj];
      }
      
      
   // now divide each coordinate sum by number of members to find mean/centroid
   // for each cluster
    for ( ii = 0; ii < k; ii++) 
      {
     /*   if (cluster_member_count[ii] == 0)
         printf("WARNING: Empty cluster %d! \n", ii);
          */
       // for each dimension
        for (jj = 0; jj < dim; jj++)
          new_cluster_centroid[ii*dim + jj] /= cluster_member_count[ii];  /// XXXX will divide by zero here for any empty clusters!

      }
  }

void get_cluster_member_count(int n, int k, int *cluster_assignment_index, int *cluster_member_count)
  {
	int ii;  
   // initialize cluster member counts
    for ( ii = 0; ii < k; ii++) 
      cluster_member_count[ii] = 0;
   // count members of each cluster    
    for ( ii = 0; ii < n; ii++){
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
        for ( kk = 0; kk < dim; kk++)
          {
            double axis_dist = X[ii*dim + kk] - cluster_centroid[cc*dim + kk]; 
            dist_sum += sqr(axis_dist);
          }
          
        double mult = ((double)cluster_member_count[cc] / (cluster_member_count[cc] + ((cluster_assignment_cur[ii]==cc) ? -1 : +1)));

        point_move_score_table[ii*dim + cc] = dist_sum * mult;
      }
  }
  
void  perform_move(int dim, int n, int k, double *X, int *cluster_assignment, double *cluster_centroid, int *cluster_member_count, int move_point, int move_target_cluster)
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
    for ( ii = 0; ii < dim; ii++)
      {
        cluster_centroid[cluster_old*dim + ii] -= (X[move_point*dim + ii] - cluster_centroid[cluster_old*dim + ii]) / cluster_member_count[cluster_old];
        cluster_centroid[cluster_new*dim + ii] += (X[move_point*dim + ii] - cluster_centroid[cluster_new*dim + ii]) / cluster_member_count[cluster_new];
      }
  }  
  
void cluster_diag(int dim, int n, int k, double *X, int *cluster_assignment_index, double *cluster_centroid)
  {
    int cluster_member_count[MAX_CLUSTERS];
    
    get_cluster_member_count(n, k, cluster_assignment_index, cluster_member_count);
     int ii, j; 
    //printf("  Final clusters %d\n", k);
    for ( ii = 0; ii < k; ii++) {
     // printf("    cluster %d:     members: %8d, centroid ", ii, cluster_member_count[ii]);
    for ( j =0 ; j<dim; j++)
	printf("%lf  ", cluster_centroid[ii*dim + j]); 
    
	printf ("\n");
	}
}
  
void copy_assignment_array(int n, int *src, int *tgt)
  {
	int ii; 
    for ( ii = 0; ii < n; ii++)
      tgt[ii] = src[ii];
  }
  
int assignment_change_count(int n, int a[], int b[])
  {
    int change_count = 0;
    int ii; 
    for ( ii = 0; ii < n; ii++)
      if (a[ii] != b[ii])
        change_count++;
        
    return change_count;
  }

void kmeans(
            int  dim,		                     // dimension of data 

            double *X,                        // pointer to data
            int   n,                         // number of elements
            
            int   k,                         // number of clusters
            double *cluster_centroid,         // initial cluster centroids
            int   *cluster_assignment_final,  // output
			int *batch_iteration
           )
  {
   	
   double *dist                    = (double *)malloc(sizeof(double) * n * k);
int   *cluster_assignment_prev = NULL; 
   cluster_assignment_prev = (int *)malloc(sizeof(int) * n);
   int   *cluster_assignment_cur  = NULL; 
   cluster_assignment_cur = (int *)malloc(sizeof(int) * n);
   
 //   double *point_move_score        = NULL;//(double *)malloc(sizeof(double) * n * k);
    
    //if (!dist || !cluster_assignment_cur || !cluster_assignment_prev /*|| !point_move_score*/)
      //fail("Error allocating dist arrays");
    
   // initial setup  
    calc_all_distances(dim, n, k, X, cluster_centroid, dist);
    choose_all_clusters_from_distances(dim, n, k, dist, cluster_assignment_cur);
    copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);

   // BATCH UPDATE
    double prev_totD = BIG_double;
    (*batch_iteration) = 0;
    while ((*batch_iteration) < MAX_ITERATIONS)
      {
	 
//	printf("batch iteration %d \n", batch_iteration);
//        cluster_diag(dim, n, k, X, cluster_assignment_cur, cluster_centroid);
        
        // update cluster centroids
         calc_cluster_centroids(dim, n, k, X, cluster_assignment_cur, cluster_centroid);

        // deal with empty clusters
        // XXXXXXXXXXXXXX

        // see if we've failed to improve
         double totD = calc_total_distance(dim, n, k, X, cluster_centroid, cluster_assignment_cur);
         if (totD > prev_totD)
          // failed to improve - currently solution worse than previous
           {
            // restore old assignments
             copy_assignment_array(n, cluster_assignment_prev, cluster_assignment_cur);
             
            // recalc centroids
             calc_cluster_centroids(dim, n, k, X, cluster_assignment_cur, cluster_centroid);
             
           //  printf("  negative progress made on this step - iteration completed (%.2f) \n", totD - prev_totD);
             
            // done with this phase
             break;
           }
           
        // save previous step
         copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);
         
        // move all points to nearest cluster
         calc_all_distances(dim, n, k, X, cluster_centroid, dist);
         choose_all_clusters_from_distances(dim, n, k, dist, cluster_assignment_cur);
         
         int change_count = assignment_change_count(n, cluster_assignment_cur, cluster_assignment_prev);
         
        // printf("%3d   %u   %9d  %16.2f %17.2f\n", batch_iteration, 1, change_count, totD, totD - prev_totD);
         fflush(stdout);
         
        // done with this phase if nothing has changed
         if (change_count == 0)
           {
		
//		printf ("ook\n"); 
            // printf("  no change made on this step - iteration completed \n");
             break;
           }

         prev_totD = totD;
                        
         (*batch_iteration)++;
      }

//cluster_diag(dim, n, k, X, cluster_assignment_cur, cluster_centroid);


      
    //printf("iterations: %3d  \n", batch_iteration/*, online_iteration*/);
      
   // write to output array
    copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_final);    
    free(dist);
    dist = NULL; 
    free(cluster_assignment_cur);
    cluster_assignment_cur = NULL; 
    free(cluster_assignment_prev);
    cluster_assignment_prev = NULL;
   // free(point_move_score);
  }           

int *mark(char *source, size_t dim, size_t N, int chunk_size){
	int  * marks= (int *) malloc (sizeof(int)*(N/chunk_size));
	char line[100000];
	FILE *src = fopen(source, "r+"); char delim[3]="\t"; 
	char *token; 
	int j = 0, mark_index =1; 
	marks[0] = ftell(src); 
	//printf ("marks[0] = %d\n", marks[0]);
	while (!feof(src)){
		fgets (line,100000, src); 
		j=(j+1)%chunk_size; 
		//if(strlen(line)>max_line_size) max_line_size = strlen(line) ;	
		if (j==0){
			marks[mark_index]=ftell(src);
			printf ("marks[%d] = %ld\n", mark_index, ftell(src));
			mark_index++; 
		}
	}
	// printf("max line size = %ld\n",max_line_size) ;

	fclose(src);
	return marks; 
}


double * getmatrix( char * source, size_t dim, size_t N, int sub, int offset, int * marks){
	FILE *src = fopen(source, "r"); 
	char line[100000]; 
	double *X=NULL;
	size_t i =0, j=0; 
	int stop=0; 

	X = (double *) malloc(sizeof(double)*sub*dim); 
	fseek(src, marks[offset/sub], SEEK_CUR); 
	char delim[3]="\t"; 
	char *token; 

	j=0; 
	while (!feof(src) && !stop){
		fgets (line,100000, src) ;
		token = strtok(line, delim) ;
		while (token != NULL){
			X[j] = atof(token) ; 
			j++;
			token = strtok(NULL, delim) ;
		}

		if (sub!=0)
			if ((j/(dim))==sub)
				stop=1; 
	} 
	fclose(src);

	return X; 
}



//  double * getmatrix_buf( char * source, size_t dim, int k , size_t N, int sub, int offset, long * marks){
	
// 	FILE *src = fopen(source, "r"); 
// 	double *X=NULL;
// 	size_t i =0, j=0; 
// 	int stop=0; 
	
// 	//decide best buffer size
// 	int MAX_DOUBLE_LEN = 24 ;
// 	int MAX_LINE_SIZE = (MAX_DOUBLE_LEN+1) * dim ;
// 	char line[MAX_LINE_SIZE]; 
	
// 	// calculate the rest size allowed for the chunk kmeans
// 	int kmeans_size = (sizeof(double)*k + 2*sizeof(int) + sizeof(int))*sub ;
// 	// the size of chunk data in file
// 	long chunk_text_size = marks[(offset/sub)+1] - marks[offset/sub] ;
// 	// printf("chunk size = %ld\t begin = %ld \t end = %ld\n",chunk_text_size,marks[(offset/sub)+1],marks[offset/sub]) ;
// 	// the average size of chunk line in file
// 	int avg_line_size = ceil(chunk_text_size/(double)sub) ;

// 	// the maximum number of lines to allow to buffer
// 	long L_max = kmeans_size/(avg_line_size+1+sizeof(char *)) ;
// 	//int L_opt = ceil( (double)(sub*L_max) / (double)(sub+L_max) ) ;
	
// 	// ajust the number of lines in buffer 
// 	int it_num = ceil((double)sub/L_max) ;
// 	// balance the number of read lines in each iteration
// 	long L_opt = ceil((double)sub / it_num) ;
	
// 	//allocate matrix
// 	X = (double *) malloc(sizeof(double)*sub*dim); 
// 	//sleep(20) ;
// 	//begin reading the file
// 	fseek(src, marks[offset/sub], SEEK_CUR); 
// 	char delim[3]="\t"; 
// 	//var to handle strtok tokens
// 	char *token ; 
// 	//store the buffered lines read form files at once
	

// 	//load the chunk
// 	j=0; 
// 	long L=0 ;
// 	long to_read_lines = sub ;

// 	// printf("stop before buffer iteration\n") ;
// 	// sleep(20) ;
	
// 	// printf("Lmax = %ld \tLopt = %ld \tavg line size = %ld \t\n",L_max,L_opt,avg_line_size+1+sizeof(char *)) ;
	
	
// 	while (!feof(src) && !stop){
// 		//allocate buffer size
// 		if(L_opt>to_read_lines)	L_opt = to_read_lines ;

// 		buff = calloc( L_opt,sizeof(char *)) ;
// 		printf("stop after buf malloc\n") ;
// 		sleep(10) ;
// 		// fill the buffer as much as possible
// 		L = 0 ;
		
// 		// size_t buffer_len = 0 ;
// 		// while(!feof(src) && L < L_opt){	
// 		// 	fgets (line,sizeof(line), src);
// 		// 	//buff[L] = malloc(strlen(line)*sizeof(char)) ;
// 		// 	//strncpy(buff[L],line,strlen(line));
// 		// 	buff[L] = malloc(avg_line_size*sizeof(char)) ;
// 		// 	// strndup(line,strlen(line)) ;
// 		// 	buffer_len += strlen(line) ;
// 		// 	L++ ;
// 		// }

// 		L=2236967 ;
// 		// printf("Lmax = %ld \tLopt = %ld \tavg line size = %ld \tbuf_size %f MB + %f MB\n",L_max,L_opt,avg_line_size+1+sizeof(char *),buffer_len/pow(1024,2),(L_opt*sizeof(char*))/pow(1024,2)) ;
// 		to_read_lines -= L ;
// 		//printf("stop after buf fill\n") ;
// 		//sleep(20) ;

// 		// convert data strings to matrix values
// 		for(long l = 0 ; l<L ; l++){
// 			//lncpy = strdup(buff[l]) ;
// 			// token = strtok(buff[l], delim);
// 			// while (token != NULL){
// 			for(int d = 0 ; d<10 ; d++){
// 				X[j] = atof("0"); 
// 				j++;
// 				// token = strtok(NULL, delim);
// 			}
			
// 			// free(buff[l]) ;
// 			// buff[l] = NULL ;
// 		}

// 		if(buff != NULL){
// 			free(buff);
// 			buff = NULL ;
// 			printf("stop after buffer free\n") ;
// 			sleep(10) ;
// 		}

// 		if (sub!=0){
// 			if ((j/(dim))==sub) {
// 				stop=1;
// 				 printf("GET MATRIX COMPLETE\n") ;
// 			}
// 		}
			
// 	} 

// 	fclose(src);
	
// 	return X; 
// }


void free_buf(char ** buf, int L){
	if(buf != NULL){
		for(int l = 0 ; l<L ; l++)
			free(buf[l]) ;
		free(buf);
		buf = NULL ;
	}	
}

double * getmatrix_bin( char * source, size_t dim, size_t N, int sub, int offset){
	FILE *src = fopen(source, "rb"); 
	printf("binary file opened\n") ;
	double *X=NULL;
	size_t i =0, j=0; 
	int stop=0; 

	X = (double *) malloc(sizeof(double)*sub*dim); 
	fseek(src,sizeof(double)*offset, SEEK_CUR);  
	printf("binary file fseek done\n") ;

	j=0; 
	while (!feof(src) && !stop){
		fread(&X[j], sizeof(double),1,src) ;
		j++ ;
		if (sub!=0)
			if ((j/(dim))==sub)
				stop=1; 
	}
	fclose(src);

	return X; 
}


double * getmatrix_mmap( char * source, size_t dim, int k , size_t N, int sub, int offset, int * marks){
	
	int src = open (source, O_RDONLY);

	double *X=NULL;
	size_t i =0, j=0; 
	int stop=0; 
	
	//decide best buffer size
	int MAX_DOUBLE_LEN = 24 ;
	int MAX_LINE_SIZE = (MAX_DOUBLE_LEN+1) * dim ;
	char line[MAX_LINE_SIZE];
	
	// calculate the rest size allowed for the chunk kmeans
	int kmeans_size = (sizeof(double)*k + 2*sizeof(int) + sizeof(int))*sub ;
	// the size of chunk data in file
	int chunk_text_size = marks[(offset/sub)+1] - marks[offset/sub] ;
	// the average size of chunk line in file
	int avg_line_size = ceil(chunk_text_size/(double)sub) ;
	// the maximum number of lines to store in buffer
	int L_max = kmeans_size/(max_line_size+1+sizeof(char *)) ;
	// ajust the number of lines in buffer 
	int it_num = ceil((double)sub/L_max) ;
	// balance the number of read lines in each iteration
	int L_opt = ceil((double)sub / it_num) ;
	
	//allocate matrix
	X = (double *) malloc(sizeof(double)*sub*dim); 
	//sleep(20) ;
	//begin reading the file
	lseek(src, marks[offset/sub], SEEK_CUR); 
	char delim[3]="\t"; 
	//var to handle strtok tokens
	char *token ; 
	//store the buffered lines read form files at once
    char * mapped;
	

	//load the chunk
	j=0; 
	int L=0 ;
	int to_read_lines = sub ;
	int s ;
	
	
	while (!stop){
		//allocate buffer size
		if(L_opt>to_read_lines)	L_opt = to_read_lines ;

		// fill the buffer as much as possible
		/* Memory-map the file. */
    	mapped = mmap (0,chunk_text_size, PROT_READ, MAP_PRIVATE,src, 0);
			
    	// check (mapped == MAP_FAILED, "mmap failed\n");
		printf("stop after buf malloc @ = %p\n",mapped) ;
		sleep(10) ;
		
		to_read_lines -= L ;

		// convert data strings to matrix values
		for(int l = 0 ; l<sub ; l++){
			for(int d = 0 ; d<10 ; d++){
				X[j] = atof("0"); 
				j++;
			}
		}

		//if we have read all the chunk lines
		if (sub!=0){
			if ((j/(dim))==sub){
				stop=1;	
				munmap(mapped, chunk_text_size);
				s = malloc_trim(0);	
				printf("stop after buf free @ = %p\n",mapped) ;
				sleep(10) ;
			}
		}		
	} 

	close(src);

	s = malloc_trim(0);	
	printf("stop after buffer free , trimming status = %d\n",s) ;
	sleep(10) ;

	printf("GET MATRIX %d COMPLETE\n",offset/sub) ;
	
	return X; 
}

double * getmatrix_buf( char * source, size_t dim, int k , size_t N, int sub, int offset, int * marks){
	
	FILE *src = fopen(source, "r"); 

	double *X=NULL;
	size_t i =0, j=0; 
	int stop=0; 
	
	//decide best buffer size
	int MAX_DOUBLE_LEN = 24 ;
	int MAX_LINE_SIZE = (MAX_DOUBLE_LEN+1) * dim ;
	char line[MAX_LINE_SIZE+1];
	
	// calculate the rest size allowed for the chunk kmeans
	int kmeans_size = (sizeof(double)*k + 3*sizeof(int))*sub ;
	
	// the size of chunk data in file
	int chunk_text_size = marks[(offset/sub)+1] - marks[offset/sub] ;
	// the average size of chunk line in file
	int avg_line_size = ceil(chunk_text_size/(double)sub) ;
	
	// the maximum number of lines to store in buffer
	int L_max = kmeans_size/(avg_line_size+1+sizeof(char *)) ;
	// ajust the number of lines in buffer 
	int it_num = ceil((double)sub/L_max) ;
	// balance the number of read lines in each iteration
	int L_opt = ceil((double)sub / it_num) ;
	
	//allocate matrix
	X = (double *) malloc(sizeof(double)*sub*dim); 
	// printf("stop after X malloc\n") ;

	//sleep(20) ;
	//begin reading the file
	fseek(src, marks[offset/sub], SEEK_CUR); 
	char delim[3]="\t"; 
	//var to handle strtok tokens
	char *token ; 
	//store the buffered lines read form files at once
	char ** buf= (char **) calloc(L_opt , sizeof(char *)) ;
	// printf("stop after buf malloc\n") ;
	

	//load the chunk
	j=0; 
	int L=0 ;
	int to_read_lines = sub ;
	int s ;
	size_t max_len ;
	
	while (!feof(src) &&!stop){
		//allocate buffer size
		if(L_opt>to_read_lines)	L_opt = to_read_lines ;
		// fill the buffer as much as possible
		L = 0 ;
		size_t buf_len = 0 ; 
		while(!feof(src) && L < L_opt){	
			fgets (line,sizeof(line), src) ;
			buf[L] = strndup(line,strlen(line)) ;
			// buf[L] = calloc((strlen(line)+1),sizeof(char)) ;
			// strncpy(buf[L],line,strlen(line)) ;
			buf[L][strlen(line)-1] = '\0' ;
			buf_len += strlen(line)+1 ;
			L++ ;
		}
		
		// printf("stop after read\n") ;
		// sleep(15) ;

		// printf("Lmax = %d \tLopt = %d \tL = %d \tread lines = %ld \tmax line size = %ld \tstrings real size %f MB + buffer pointers size %f MB\n",L_max,L_opt, L,(j/(dim)),max_line_size+1+sizeof(char *),buf_len/pow(1024,2),(L_opt*sizeof(char*))/pow(1024,2)) ;
		
		L=L_opt ;
		to_read_lines -= L ;

		// convert data strings to matrix values
		for(int l = 0 ; l<L ; l++){
			// for(int d = 0 ; d<dim ; d++){
			// 	X[j] = atof("0"); 
			// 	j++;
			// 	// token = strtok(NULL, delim);
			// }
			token = strtok(buf[l],delim);
			while (token != NULL){
				X[j] = atof(token); 
				j++;
				token = strtok(NULL,delim);
			}
			free(buf[l]) ;
			buf[l]=NULL ;
		}

		// printf("stop after fill\n") ;
		// sleep(15) ;

		//if we have read all the chunk lines
		if (sub!=0){
			if ((j/(dim))==sub){
				stop=1;	
				s = malloc_trim(0);	
			}
		}		
	
	} 

	if(buf != NULL){
		free(buf);
		buf = NULL ;
	}
	
	fclose(src);

	s = malloc_trim(0);	

	// printf("stop after free\n") ;
	// sleep(10) ;
	
	return X; 
}

double * random_center_init(double * X, size_t N, size_t dim, size_t k){
	int j=0, i =0;
	int r, new, l, regen;  
	int *gen =(int *) malloc(sizeof(int)*k); 

	double *c = (double *)malloc (sizeof(double)*(k*dim));
	//iterate over k cluster
	while (j<k){
		r = rand()%N;
		//if not the first cluster
		if (j!=0){
			//search a new point to assign to cluster j
			new =0; 
			while(!new ){
				l = 0; 
				regen =0; 
				//find if the randomly chosen point was generated
				// to past cluster 
				while(l<j && !regen){
					if (gen[l]==r)
						regen =1;
					l++; 
				}
				//if so , reassign a new point to cluster j
				if (regen)
					r = rand()%N; 
	    		else 
					new =1; 
			}
		} 
		//asign randomly chosen point r to cluster i center
		gen[j] = r; 
		int m; 
		for ( m=0; m<dim;m++){
		c[i+m] = X[r*dim+m];
		}
		//goto next cluster
		i=i+dim; 
		j++;  

	}
	free(gen);
	gen = NULL;
	return c; 
} 

double * random_center_init_one(double * X, size_t N, size_t dim){
	int  i =0;
	int r, m;  
	double *c = malloc (dim*sizeof(double)); 
		r = rand()%N;
		for ( m=0; m<dim;m++){
		c[i+m] = X[r*dim+m];
		}
	return c; 
} 

int token_center(int *centers, int nb, int center){
	int found=0, i=0; 
	while (i<nb && !found){
		if (centers[i]==center)
			found =1; 
		i++; 
	}
	return found; 
}

int find_best_center(double * distance_cur_center, int *centers, size_t N, size_t dim, size_t k, double sum){
	int best, i; 
	double max=0;
	//iterate over points
	for (i=0; i<N; i++){
		//if the distance point-center/sumof distance is bigger than max and the current center of the point is not taken
		if (distance_cur_center[i]/sum>max && !token_center(centers, k,i)){
			//reassign the best found center and the max distance
			best = i; 
			max = distance_cur_center[i]/sum; 
		}
	}
	return best; 
}

double affect_x(double * X, int *centers_int, int index_x, int k, size_t dim){
	double min = BIG_double; 
	int j; 
	for (j=0; j<k;j++){
		if (calc_distance(dim, &X[centers_int[j]*dim], &X[index_x*dim])<min){
			min = calc_distance(dim, &X[centers_int[j]*dim], &X[index_x*dim]); 
		}
	}
	return min; 
} 

double * kmeans_init_plusplus(double *X, size_t N, size_t dim, size_t k){
	double * centers = (double *) malloc (sizeof(double)*k*dim) ; 
	double * distance_cur_center = (double *) malloc (sizeof(double)*N) ; 
	int * centers_int = (int *) malloc (sizeof(int)*k) ; 
	double sum =0 ; 
	//int first = rand()%N ; 
	int first = chunk_index++ ;

	int i, j, best;  
	centers_int[0] = first; 
	for (i=1; i<k; i++){
		for (j=0; j<N; j++){
			distance_cur_center[j] = affect_x( X, centers_int, j, i, dim); 
			sum+=distance_cur_center[j]; 	
		}
		 centers_int[i]= find_best_center(distance_cur_center,centers_int, N, dim, i, sum); 
	}
	for (i=0; i<k;i++){
		for (j=0; j<dim; j++){
				centers[i*dim+j]=X[centers_int[i]*dim+j]; 
			//	printf ("%lf\t", X[centers_int[i]*dim+j]); 
		}	
		//printf ("\n"); 
		
	}
	free(centers_int); 
	centers_int =NULL; 
	free(distance_cur_center);
	distance_cur_center = NULL;  

	return centers; 
}

void dispersion (double * X, size_t dim, size_t N){
	int i,j; 
	double *min = malloc(sizeof(size_t)*dim);
	double *max = malloc(sizeof(size_t)*dim); 

	for (i=0; i<dim; i++){
	min[i]=X[i];  
	}
	for (i=0; i<dim; i++){
	max[i]=X[i];  
	}

	//iterate over point from the second to find the max
	double d=0; 
	for (i=dim; i<N; i++){
	   if (calc_distance (dim, min, X+i*dim) > d){
		for(j=i;j<i+dim;j++){
			max[j-i] = X[j]; 
		}
		d = calc_distance (dim, min, X+i*dim); 
		} 
	}
	/*for (i=0; i<dim; i++){
	printf ("%lf\t", max[i]); 
	}
	printf("\n"); 
	for (i=0; i<dim; i++){
	printf ("%lf\t", min[i]); 
	}*/
	printf ("distance min-max %lf\n", calc_distance(dim, min, max));
}

//partial kmeans function 
double * partial_mean (double * X, size_t dim, size_t chunk_len, size_t chunk_index){
	int i,j,k; 
	double *sum = malloc(sizeof(double )*dim);

	for (i=0;i<dim;i++)
	sum[i]=0;  
	for (i=0;i<chunk_len;i++){
		for (j=0;j<dim;j++){
			k = (chunk_index*dim)+(i*dim)+j; 
			sum[j] = sum[j]+X[k]; 
		}
		
	} 

	for (i=0;i<dim;i++){
	sum[i]/=chunk_len;
	// ("\t %lf", sum[i]); 
	}
	//printf("\n"); 

	return sum; 
}

double  sse_calculate (double * X, int dim, double *  centroids, int *  assignement , size_t n) {
	double sse = 0; 
	int i; 
	for (i=0; i<n; i++){
		sse+= calc_distance (dim, &centroids[assignement[i]*dim], &X[i*dim]);  
	}
	return sse; 
}

double  var_calculate (double * X, int dim, double *  centroids, int center, int *  assignement , size_t n) {
	double sse = 0; 
	int i; 
	for (i=0; i<n; i++){
		if (assignement[i]==center){
		sse+= calc_distance (dim, &centroids[center*dim], &X[i*dim]);  
		}
	}
	return sse; 
}

double inter_centroid_distance (int dim, int k, double *  centroids){
	double sum =0; 
	int i, j; 
	for (i=0; i< k; i++){
		for (j=i+1; j<k; j++){
			sum+=calc_distance(dim, &centroids[i*dim], &centroids[j*dim]); 
		}	
	}
	return sum; 	
}

double *getsubmatrix(double *X, int dim, int offset, int len){
	double *Y = (double *)malloc ((dim*len)*sizeof(double)); 
	int i;
	for ( i=offset*dim; i <dim*len+offset*dim;i++){
		Y[i-offset*dim] = X[i]; 
	}
	return Y; 
}

int existJ (int k, int *tab, int j, int max){
	int i=0; 
	while (i<k && tab[i]!=j && i<max){
		i++; 
	}
	if (i<k && i!=max)
		return 1; 
	else 
		return 0; 
}

int maximize(double *partial_means, int * jpris, int pris, int dim, int k, int chunks){
	int i,m, res; 
	double sum, maxsum;
	maxsum =0; 
	res =0; 
	for (i=0; i<chunks; i++){ //pour tout chunk
		sum =0; 
		
		if (!existJ(k, jpris,i, pris)){//s'il fait pas partie des chunks sélectionnés
			//for (j=0; j<chunks;j++){
				for (m=0;m<pris;m++){
					sum+= calc_distance(dim, &partial_means[i*dim], &partial_means[jpris[m]*dim]); 
				}	
			//}
				if (sum>maxsum){
					maxsum=sum; 
					res = i; 
				}
		}
		

	}

		printf ("stop %d\n", res); 
	return res; 
}

int existG (groupe * groupes, int nb, int g){
	int found =0, i=0, j; 
	//iterate over the groupes
	while (i<nb && !found){
		j=0; 	
		//iterate over the chunks of the group
		while (j<groupes[i].nb && !found){
			if (groupes[i].means[j]==g){
				found=1; 
			}
			j++; 
		}
		i++; 
	}
	return found; 
}

int existGi (groupe  groupes, int g){
	int found =0, j=0; 	
	while (j<groupes.nb && !found){
		if (groupes.means[j]==g){
				found=1; 
		}
		j++; 
	}
	return found; 
}

//version 1: pour séléctionner aléatoirement un groupe parcourir toute la structure cluster_assignement jusqu'à tomber sur le point qui correspond au cluster 
void form_index ( int *cluster_assignement, int nb_groupes, groupe *grp, size_t N){
	int i =0, j=0, affected, *index_grp;
	index_grp = (int *)malloc (sizeof(int) *nb_groupes +1);  
	for (i=0; i<nb_groupes; i++){
	    grp[i].members = (int *) malloc (sizeof(int)*grp[i].nb_members); 
	    index_grp[i]= 0; 
	}

	//iterate over points
	for (i=0; i< N; i++){
		affected = 0; 
		j =0;
		// affect the point to a group
		while (!affected && j<nb_groupes){
			//if the point cluster is exist in the group j
			if (existGi (grp[j], cluster_assignement[i])){
				grp[j].members[index_grp[j]]=i; 
				index_grp[j]++; 
				affected =1; 
			}
		j++; 
		}
	}

	free(index_grp); 
} 

double * form_chunk( groupe *grp, /*double *X*/char * source, int * marks, int * cluster_assignment, int nb_groupes, size_t N, size_t dim, int k , size_t taille){
	double * chunk =(double *) malloc (sizeof(double )*(dim*taille)); 
	int nb_samples; 
	int j, index_groupe, l, size =0,i, found, m, *int_chunk; 
	int tmp; 
	double *Y; 
	float c; 
	char ** buf ;

	index_chunk_elem  index[N/taille];
	index_chunk_elem  * cur [N/taille];
	index_chunk_elem * elem; 
	for (i=0; i<N/taille;i++){
	index[i].element = -1; 
	index[i].next = NULL;
	cur[i]=NULL;
	
	}
	//form groups index
	form_index (cluster_assignment, nb_groupes, grp,  N); 
	//form chunk
	while(size < taille){
	//iterate over groups
	for ( i=0; i<nb_groupes; i++){
			//get the elements ratio that can be in chunk for the group
				c = (float)grp[i].nb_members /(float) N * taille;
				nb_samples = (int)c; 
				j = 0; 
				//choose the elements in the chunk randomly
				for (j=0; j<nb_samples;j++){
				
				//tmp = rand()%grp[i].nb_members;
				tmp = (j+i*nb_samples)%grp[i].nb_members;
				
				if (size==taille) 
						break;
				//if the slected point from the group does not exist
				//in the index , add it
				if (index[grp[i].members[tmp]/taille].element == -1){
					index[grp[i].members[tmp]/taille].element=grp[i].members[tmp]%taille; 
				}
				else{
	 
				elem =(index_chunk_elem*) malloc(sizeof(index_chunk_elem));
				elem->element = grp[i].members[tmp]%taille; 
				
				elem->next = NULL; 
					if(cur[grp[i].members[tmp]/taille] !=NULL)
						cur[grp[i].members[tmp]/taille]->next = elem; 
					else{
						cur[grp[i].members[tmp]/taille] =(index_chunk_elem*) malloc(sizeof(index_chunk_elem));
						index[grp[i].members[tmp]/taille].next = elem; 
				}
				cur[grp[i].members[tmp]/taille]=elem; 
				}
				//free(elem);	
				
			//}
			size++; 
			}
			
	}
	}
	size=0; 
	//iterate over chunk points
	for (i=0; i<N/taille;i++){
		if (index[i].element!=-1){
			// Y = getmatrix(source, dim, taille,taille,i*taille, marks); 
			Y = getmatrix_buf(source, dim,k, taille,taille,i*taille, marks); 
			//printf ("chunk %d\n", i); 
			for (j=0; j<dim;j++){
				chunk[size*dim+j] = Y[index[i].element*dim+j]; 
				//printf("%lf\n", Y[index[i].element*dim+j]);
			}
			size++; 
			//goto next point in index
			cur [i]= index[i].next; 	
			while(cur[i]!=NULL){
				//printf("ELEMENT %ld\n", cur[i]->element);
				for (j=0; j<dim;j++){
					chunk[size*dim+j] = Y[cur[i]->element*dim+j]; 
				}
				cur[i]=cur[i]->next; 
				size++; 
			}
			free(cur[i]);
			free(Y);
		}
	}
	//printf ("size %d\n", size);
	return chunk; 
}

void sort_centroid(size_t dim, size_t k, double * centroid){
	int i, j, l, index=0; 
	double *tmp = (double *) malloc (sizeof(double)*k*dim); 

	//iterate over clusters
	for (i =0; i< k; i++){
		//initialize
		for (j=0; j< dim; j++){
			tmp[j] = 0; 
		}
		//iterate over next clusters
		for (j=i;j<k; j++){
			if (tmp[0]<centroid[dim*j]){
				for (l=0; l< dim; l++){
					tmp[l] = centroid[j*dim+l];
					index = j;  
				}	
			}
		}
		//iterate over dimensions and swap the centers
		for (l=0; l< dim; l++){
			centroid[index*dim+l] = centroid[i*dim+l]; 
			centroid[i*dim+l] = tmp[l];  
		}
	}

	free(tmp); 
}

void kmeans_by_chunk( char * source, size_t dim, int taille, int N, int k){
	int i, j , n , m=0 , kmeans_iterations = 0 ;
	double * Y=NULL; 
	int * cluster_assignment_final=NULL, * cluster_assignment_final_Y=NULL; 
	double * cluster_centroid=NULL, *cluster_centroid_by_chunk=NULL , *chunk_centroid = NULL ; 
	double * var=NULL; 
	int cluster_member_count[MAX_CLUSTERS]; //à modifier pour une allocation dynamique 
	int nb_groupes=0; 
	groupe * groupes=NULL; 
	
	/****** PHASE 1 : PARTIELS CHUNKS KMEANS ******/

	//mark the chunks in dataset
	save_phase_time(1,1,0,-1);
	int *marks = mark(source, dim, N, taille);
	save_phase_time(1,1,1,-1);
	
	//apply kmeans on each chunk
	if (N/taille >1){
		save_phase_time(1,2,0,0);
		cluster_assignment_final = (int *) malloc(N*sizeof(int));
		cluster_centroid_by_chunk = (double *)malloc(sizeof(double)*(k*dim*N/taille)); 
		var = (double *) malloc (sizeof(double)*k*N/taille); 

		for( m = 0; m<N/taille; m++){
			//lecture du chunk 
			i = m*taille;
			Y = getmatrix_buf(source, dim,k, taille,taille,i, marks); 
			// Y = getmatrix(source, dim,taille,taille,i, marks); 
			// Y = getmatrix_bin(source, dim,taille,taille,i); 
			// Y = getmatrix_mmap(source, dim,k, taille,taille,i, marks);

			save_phase_time(1,2,(m+1),1);	//get matrix
			
			// chunk m init kmeans++ 
			cluster_centroid = kmeans_init_plusplus(Y, taille, dim, k);
			save_phase_time(1,3,(m+1),1);

			//apply kmeans on the chunk m
			cluster_assignment_final_Y=  (int *) malloc(taille*sizeof(int));
			kmeans(dim,Y,taille, k, cluster_centroid, cluster_assignment_final_Y,&kmeans_iterations);
			save_phase_time(1,4,(m+1),kmeans_iterations);

			//variance calculation on the chunk m
			for ( j=0; j<k;j++){
				var[m*k+j] = var_calculate (Y, dim, cluster_centroid, j, cluster_assignment_final_Y , taille)/taille; 
			}
			//copy the centroids to global array
			for ( j=0; j<k*dim;j++){
				cluster_centroid_by_chunk[m*k*dim+j] = cluster_centroid[j];   
			}
			//copy points clustering to global array
			for ( n=0 ; n<taille; n++){		
				cluster_assignment_final[m*taille+n]=cluster_assignment_final_Y[n]+m*k; 
			}

			
			
			free(Y);
			Y = NULL;

			free(cluster_centroid); 
			cluster_centroid=NULL; 
			
			free(cluster_assignment_final_Y); 
			cluster_assignment_final_Y=NULL; 

			save_phase_time(1,5,(m+1),1);

			// printf("stop after get matrix complete\n") ;
			// sleep(10) ;
			
		}

		/****** PHASE 2 : PARTIELS CLUSTERS GROUPING ******/
		save_phase_time(2,1,0,0);
		groupes = (groupe * ) malloc (sizeof(groupe )*N/taille*k);
		nb_groupes =0; 
		for (j=0; j< k*N/taille; j++){
			if (!existG(groupes, nb_groupes, j)){
				groupes[nb_groupes].means = (int*) malloc(sizeof(int)); 
				groupes[nb_groupes].nb =1; 		
				groupes[nb_groupes].means[0]=j; 
				for ( n = j+1; n<k*N/taille; n++){
					//condition de non existence
					if (calc_distance(dim, &cluster_centroid_by_chunk[j*dim], &cluster_centroid_by_chunk[n*dim])<var[j] && !existG(groupes, nb_groupes, n)){
						groupes[nb_groupes].means = (int *) realloc (groupes[nb_groupes].means, (groupes[nb_groupes].nb+1)*sizeof(int));  	
						groupes[nb_groupes].means[groupes[nb_groupes].nb] = n; 
						//printf ("next\n");
						groupes[nb_groupes].nb++; 
					
					}
				}
				nb_groupes++; 
			}
		}

		save_phase_time(2,1,1,j);
		
		free(cluster_centroid_by_chunk); 
		cluster_centroid_by_chunk = NULL;

		save_phase_time(2,2,0,0);
		//groups members count update
		float count = k*N/taille; 
		get_cluster_member_count(((int) (N/taille)) * taille, (int) count, cluster_assignment_final, cluster_member_count);

		for (j=0; j<nb_groupes; j++){
			groupes[j].nb_members = 0; 
			for (n = 0; n< groupes[j].nb; n++){
				groupes[j].nb_members+=cluster_member_count[groupes[j].means[n]]; 
			}
		}

		save_phase_time(2,2,1,nb_groupes);

		/****** PHASE 3 : FINAL CHUNK BUILDING ******/
		save_phase_time(3,1,0,-1);
		double * chunk; 
		chunk = form_chunk(groupes,source, marks,cluster_assignment_final, nb_groupes, N, dim,k, taille); 
		save_phase_time(3,1,1,-1);

		free(cluster_assignment_final); 
		cluster_assignment_final =NULL; 
		for (i = 0; i<nb_groupes; i++){
			free(groupes[i].means);
			free(groupes[i].members); 
		}
		free(groupes); 
		groupes = NULL;
		
		/****** PHASE 4 : FINAL CHUNK KMEANS ******/

		save_phase_time(4,1,0,-1);
		cluster_centroid =kmeans_init_plusplus(chunk, taille, dim, k);
		save_phase_time(4,1,1,-1);
		
		save_phase_time(4,2,0,0);
		cluster_assignment_final_Y=  (int *) malloc(taille*sizeof(int));
		kmeans(dim,chunk,taille, k, cluster_centroid, cluster_assignment_final_Y,&kmeans_iterations);
		save_phase_time(4,2,1,kmeans_iterations);

  		r8mat_write ("results/centers.csv",dim,k,cluster_centroid) ;
		free(chunk);
		chunk = NULL;
	}
	else{
		//sleep(20) ;
		save_phase_time(2,1,0,-1);
		double* X; 
		// X = getmatrix(source, dim, N,N,0, marks );
		// rbinmat_write ("points_bin.csv",dim,N,X) ;
		// X = getmatrix_mmap(source, dim,k,N,N,0,marks);
		// X = getmatrix_bin(source,dim,N,N,0); 

		X = getmatrix_buf(source, dim,k, N,N,0, marks);
		save_phase_time(2,1,1,-1);

		save_phase_time(2,2,0,0);
		cluster_centroid =kmeans_init_plusplus(X, N, dim, k);
		save_phase_time(2,2,1,1);

		save_phase_time(2,3,0,0);	
		cluster_assignment_final=  (int *) malloc(N*sizeof(int)); 
		kmeans(dim,X,N, k, cluster_centroid, cluster_assignment_final,&kmeans_iterations);
		save_phase_time(2,3,1,kmeans_iterations);

		free(X);
		X = NULL;
	}

	free(cluster_assignment_final_Y);
	cluster_assignment_final_Y = NULL; 
	free(cluster_centroid); 	
	cluster_centroid =NULL;
	free(marks) ;
	marks=NULL ;

}

void kmeans_kmeans(double * X,int dim, int taille, int N, int k){ //expérimentation de la méthode de kmeans de kmeans sans formation de chunk représentatif
	int i, * cluster_assignment_final, *cluster_assignment_final_Y, m, j, n; 
	double * Y; 
	double * cluster_centroid=NULL, *cluster_centroid_by_chunk=NULL; 
	double *var; 
	int nb_groupes; 
	int cluster_member_count[MAX_CLUSTERS]; //à modifier pour une allocation dynamique 
	int kmeans_iteration = 0;
	    
	groupe *groupes; 
	cluster_assignment_final=  (int *) malloc(N*sizeof(int));
	cluster_assignment_final_Y=  (int *) malloc(taille*sizeof(int));
	cluster_centroid = malloc(sizeof(double)*k*dim);
	cluster_centroid_by_chunk = malloc(sizeof(double)*k*dim*N/taille); 
	var = malloc (sizeof(double)*k*N/taille); 
	groupes = malloc (sizeof(groupe )*N/taille*k); 
	for( m = 0; m<N/taille; m++){
		i = m*taille;
		Y= getsubmatrix(X, dim,i, taille);
		//cluster_centroid = random_center_init(Y, taille, dim, k);
		cluster_centroid =kmeans_init_plusplus(Y, taille, dim, k);	        
		kmeans(dim,Y,taille, k, cluster_centroid, cluster_assignment_final_Y,&kmeans_iteration);
		for ( j=0; j<k;j++){
			var[m*k+j] = var_calculate (Y, dim, cluster_centroid, j, cluster_assignment_final_Y , taille)/taille; 
		}
		for ( j=0; j<k*dim;j++){
			cluster_centroid_by_chunk[m*k*dim+j] = cluster_centroid[j];   
		}
		for ( n=0 ; n<taille; n++){		cluster_assignment_final[m*taille+n]=cluster_assignment_final_Y[n]+m*k; 
		}
		
	}
	cluster_centroid = random_center_init(cluster_centroid_by_chunk, taille, dim, k);
	//cluster_centroid =kmeans_init_plusplus(cluster_centroid_by_chunk, taille, dim, k);
		
	kmeans(dim,cluster_centroid_by_chunk,k*N/taille, k, cluster_centroid, cluster_assignment_final_Y,&kmeans_iteration);

}

int main (int argc, char **argv){

	cpu_set_t set;
	// clear cpu mask
	CPU_ZERO(&set);
	// set cpu 0
	CPU_SET(0, &set);
	// 0 is the calling process
	sched_setaffinity(0, sizeof(cpu_set_t), &set);
	//set priority
	setpriority(PRIO_PROCESS, 0, -20);

	clock_t begin, end; 
	int i,j; 
	size_t  tot;
	int taille_chunk, taille; 
	size_t dim, N; 
	double *X;

	char cmd[200];  
	sprintf (cmd, "echo %d > /sys/fs/cgroup/memory/kmeans/cgroup.procs", getpid()); // /cgroups/mem/kmeans/tasks
	//printf ("%s\n", cmd);
	system(cmd); 

	sprintf (cmd, "echo %d > /sys/fs/cgroup/blkio/kmeans/cgroup.procs", getpid()); // /cgroups/mem/kmeans/tasks
	system(cmd); 

	sleep(1);
	
	char * source = argv[1]; 
	size_t k = atoi(argv[2]);
	N = atoi(argv[3]); 
	taille = atoi(argv[4]);
	dim = atoi(argv[5]); 
	taille_chunk = taille; 
	srand(time(NULL));
	//begin = clock();  

	struct timeval tv_begin, tv_end;
	//gettimeofday(&tv_begin, NULL); 
	kmeans_by_chunk( source, dim, taille, N, atoi(argv[2])); 
	//X= getmatrix(source,dim,N,0,0,NULL) ; 
	//kmeans_kmeans( X, dim, taille,N, atoi(argv[2])); 
	//gettimeofday(&tv_end, NULL); 

	//end = clock(); 
	//printf ("KMEANS TIME = %lf min\n", (float)(tv_end.tv_sec - tv_begin.tv_sec)/60);
	
	system("pkill watch");
	system("killall -9 prog_script_cpu_io");
	
	return 0;
}
