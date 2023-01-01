// kmeans.c
// Ethan Brodsky
// modified by Camélia SLIMANI

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <sys/time.h>

#define sqr(x) ((x)*(x))

#define MAX_CLUSTERS 100000

#define MAX_ITERATIONS 150

#define BIG_double (INFINITY)


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
    printf(str);
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
            int   *cluster_assignment_final  // output
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
    int batch_iteration = 0;
    while (batch_iteration < MAX_ITERATIONS)
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
                        
         batch_iteration++;
      }

//cluster_diag(dim, n, k, X, cluster_assignment_cur, cluster_centroid);

/* THe online update prtion of this code has never worked properly, but batch update has been adequate for our projects so far
   // ONLINE UPDATE
    int online_iteration = 0;
    int last_point_moved = 0;
    
    int cluster_changed[MAX_CLUSTERS];
    for (int ii = 0; ii < k; ii++)
      cluster_changed[ii] = 1;
    
    int cluster_member_count[MAX_CLUSTERS];
    get_cluster_member_count(n, k, cluster_assignment_cur, cluster_member_count);
    
    while (online_iteration < MAX_ITERATIONS)
      {
//        printf("online iteration %d \n", online_iteration);

       // for each cluster
        for (int ii = 0; ii < k; ii++)
          if (cluster_changed[ii])
            update_delta_score_table(dim, n, k, X, cluster_assignment_cur, cluster_centroid, cluster_member_count, point_move_score, ii);
            
       // pick a point to move
       // look at points in sequence starting at one after previously moved point
        int make_move = 0;
        int point_to_move = -1;
        int target_cluster = -1;
        for (int ii = 0; ii < n; ii++)
          {
            int point_to_consider = (last_point_moved + 1 + ii) % n;
              
           // find the best target for it
            int best_target_cluster = -1;
            int best_match_count    = 0;
            double best_delta        = BIG_double;
            
           // for each possible target
            for (int jj = 0; jj < k; jj++)
              {
                double cur_delta = point_move_score[point_to_consider*k + jj];

               // is this the best move so far?
                if (cur_delta < best_delta)
                 // yes - record it
                  {
                    best_target_cluster = jj;
                    best_delta = cur_delta;
                    best_match_count = 1;
                  }
                else if (cur_delta == best_delta)
                 // no, but it's tied with the best one
                 best_match_count++;
              }

           // is the best cluster for this point its current cluster?
            if (best_target_cluster == cluster_assignment_cur[point_to_consider])
             // yes - don't move this point
               continue;

           // do we have a unique best move?
            if (best_match_count > 1)
             // no - don't move this point (ignore ties)
              continue;
            else
             // yes - we've found a good point to move
              {
                point_to_move = point_to_consider;
                target_cluster = best_target_cluster;
                make_move = 1;
                break;
              }
          }

        if (make_move)
          {
           // where should we move it to?            
            printf("  %10d: moved %d to %d \n", point_to_move, cluster_assignment_cur[point_to_move], target_cluster);

           // mark which clusters have been modified          
            for (int ii = 0; ii < k; ii++)
              cluster_changed[ii] = 0;
            cluster_changed[cluster_assignment_cur[point_to_move]] = 1;
            cluster_changed[target_cluster] = 1;

           // perform move
            perform_move(dim, n, k, X, cluster_assignment_cur, cluster_centroid, cluster_member_count, point_to_move, target_cluster);

           // count an iteration every time we've cycled through all the points
            if (point_to_move < last_point_moved)
              online_iteration++;

            last_point_moved = point_to_move;
          }

      }

*/
      
    printf("iterations: %3d  \n", batch_iteration/*, online_iteration*/);
      
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
int  * marks= (int *) malloc (sizeof(int)*(N/chunk_size+1));
char line[100000];
FILE *src = fopen(source, "r+"); char delim[3]=" "; 

char *token; 
int j = 0, mark_index =1; 
marks[0] = ftell(src); 
while (!feof(src)){
	fgets (line,100000, src); 
	j=(j+1)%chunk_size; 
	if (j==0){
		marks[mark_index]=ftell(src);
		printf ("marks[%d] = %d\n", mark_index, ftell(src)); 
		mark_index++; 
	}
}

fclose(src);
return marks; 
}

double * getmatrix( char * source, size_t dim, size_t N, int sub, int offset, int * marks){
	FILE *src = fopen(source, "r+"); 
	char line[100000]; 
	double *X=NULL;
	size_t i =0, j=0; 
	int stop=0; 

	X = (double *) malloc(sizeof(double)*sub*dim); 
	fseek(src, marks[offset/sub], SEEK_CUR); 
	char delim[3]=" "; 
	char *token; 

	j=0; 
	while (!feof(src) && !stop){
		fgets (line,100000, src); 
		token = strtok(line, delim);
			while (token != NULL){
					X[j] = atof(token); 
					j++;
			token = strtok(NULL, delim);
		}

	if (sub!=0)
		if ((j/(dim))==sub)
			stop=1; 
	} 
	fclose(src);
	//*N=j/(*dim);
	//*d = *dim; 
	//free(dim);
	//dim = NULL;
	return X; 
}

double * random_center_init(double * X, size_t N, size_t dim, size_t k){
	int j=0, i =0;
	int r, new, l, regen;  
	int *gen =(int *) malloc(sizeof(int)*k); 

	double *c = (double *)malloc (sizeof(double)*(k*dim));
	while (j<k){
		r = rand()%N;
		if (j!=0){
			new =0; 
			while(!new ){
				l = 0; 
				regen =0; 
				while(l<j && !regen){
					if (gen[l]==r)
						regen =1;
					l++; 
				}
				if (regen)
					r = rand()%N; 
	    			else 
					new =1; 
			}
		} 
		gen[j] = r; 
		int m; 
		for ( m=0; m<dim;m++){
		c[i+m] = X[r*dim+m];
		}
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
	for (i=0; i<N; i++){
		if (distance_cur_center[i]/sum>max && !token_center(centers, k,i)){
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
	double * centers = (double *) malloc (sizeof(double)*k*dim); 
	double * distance_cur_center = (double *) malloc (sizeof(double)*N); 
	int * centers_int = (int *) malloc (sizeof(int)*k); 
	double sum =0; 
	int first = rand()%N; 
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
	printf ("\t %lf", sum[i]); 
	}
	printf("\n"); 

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
	while (i<nb && !found){
		j=0; 	
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

	for (i=0; i< N; i++){
		affected = 0; 
		j =0;
		while (!affected && j<nb_groupes){
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


	double * form_chunk( groupe *grp, /*double *X*/char * source, int * marks, int * cluster_assignment, int nb_groupes, size_t N, size_t dim, size_t taille){
	double * chunk =(double *) malloc (sizeof(double )*(dim*taille)); 
	int nb_samples; 
	int j, index_groupe, l, size =0,i, found, m, *int_chunk; 
	int tmp; 
	double *Y; 
	float c; 
	index_chunk_elem  index[N/taille];
	index_chunk_elem  * cur [N/taille];
	index_chunk_elem * elem; 
	for (i=0; i<N/taille;i++){
	index[i].element = -1; 
	index[i].next = NULL;
	cur[i]=NULL;
	}
	form_index (cluster_assignment, nb_groupes, grp,  N); 
	while(size < taille){
	for ( i=0; i<nb_groupes; i++){
			
				c = (float)grp[i].nb_members /(float) N * taille;
				nb_samples = (int)c; 
				j = 0; 
				for (j=0; j<nb_samples;j++){
				
				tmp = rand()%grp[i].nb_members;
				
				if (size==taille) 
						break;
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
	for (i=0; i<N/taille;i++){
	if (index[i].element!=-1){
	Y = getmatrix(source, dim, taille,taille,i*taille, marks); 
	//printf ("chunk %d\n", i); 
	for (j=0; j<dim;j++){
	chunk[size*dim+j] = Y[index[i].element*dim+j]; 
	//printf("%lf\n", Y[index[i].element*dim+j]);
	}
	size++; 
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
	printf ("size %d\n", size);
	return chunk; 
}

void sort_centroid(size_t dim, size_t k, double * centroid){
int i, j, l, index=0; 
double *tmp = (double *) malloc (sizeof(double)*k*dim); 

for (i =0; i< k; i++){
for (j=0; j< dim; j++){
tmp[j] = 0; 
}
	for (j=i;j<k; j++){
		if (tmp[0]<centroid[dim*j]){
			for (l=0; l< dim; l++){
				tmp[l] = centroid[j*dim+l];
				index = j;  
			}	
		}
	
	}
	for (l=0; l< dim; l++){
		centroid[index*dim+l] = centroid[i*dim+l]; 
		centroid[i*dim+l] = tmp[l];  
	}

}

free(tmp); 
}

void compare_solution ( double * centroid, char * solution_file_name, size_t dim, size_t k){

FILE *solution = fopen (solution_file_name, "r+"); 
char *token; 
char tmp[1000000]; 
char delim[3] = " \n"; 
double * centroid_solution = (double *) malloc (sizeof(double) *k*dim);
int index=0, index2;  
double delta=0; 
while(!feof(solution)){
	fgets (tmp, 1000000, solution); 
	token = strtok(tmp, delim); 
	while (token != NULL){
		centroid_solution[index] = atof(token); 
		index++;
		token = strtok(NULL, delim); 	
	}
}
fclose(solution); 
for (index = 0; index<k; index++){
	printf ("centres %d\n", index); 
	for (index2=0; index2<dim; index2++){
	 printf ("sol: %lf  sol obt : %lf\n", centroid_solution[index*dim+index2], centroid[index*dim+index2]); 
	}
	
	delta +=calc_distance(dim,&centroid_solution[index*dim], &centroid[index*dim]);  	
}
printf ("delta solution_type solution obtenue : %lf\n", sqrt(delta)); 

free(centroid_solution); 
}

void kmeans_by_chunk(/*double * X9*/ char * source, size_t dim, int taille, int N, int k){
	int i, * cluster_assignment_final=NULL, *cluster_assignment_final_Y=NULL, j,n; 
	double * Y=NULL; 
	double * cluster_centroid=NULL, *cluster_centroid_by_chunk=NULL; 
	double *var=NULL; 
	int nb_groupes,m; 
	int cluster_member_count[MAX_CLUSTERS]; //à modifier pour une allocation dynamique 
	double* X; 
	struct timeval tv_begin, tv_end;
	groupe *groupes=NULL; 
	size_t tot; 
	int *marks = mark(source, dim, N, taille); 
	//appliquer le kmeans sur chaque chunk 
	if (N/taille >1){
		cluster_assignment_final=  (int *) malloc(N*sizeof(int));
		cluster_centroid_by_chunk = (double *)malloc(sizeof(double)*(k*dim*N/taille)); 
		var = (double *) malloc (sizeof(double)*k*N/taille); 
		groupes = (groupe * ) malloc (sizeof(groupe )*N/taille*k); 
		for( m = 0; m<N/taille; m++){
			printf ("chunk %d/%d\n", m, N/taille); 
			i = m*taille;
			//lecture du chunk 
			Y = getmatrix(source, dim, taille,taille,i, marks); 
			printf ("unmap error"); 
			//Y= getsubmatrix(X, dim,i, taille);
			//initialisation aléatoire des centres pour ce chunk  
			//cluster_centroid = random_center_init(Y, taille, dim, k);
			
			gettimeofday(&tv_begin, NULL); 	
			cluster_centroid =kmeans_init_plusplus(Y, taille, dim, k);
			gettimeofday(&tv_end, NULL); 
			printf ("init time: %lf\n", (float)(tv_end.tv_sec - tv_begin.tv_sec)/60);
			cluster_assignment_final_Y=  (int *) malloc(taille*sizeof(int));
			//printf ("chunk %d\n", m); 
			kmeans(dim,Y,taille, k, cluster_centroid, cluster_assignment_final_Y);
			
			for ( j=0; j<k;j++){
				var[m*k+j] = var_calculate (Y, dim, cluster_centroid, j, cluster_assignment_final_Y , taille)/taille; 
			}
			for ( j=0; j<k*dim;j++){
				cluster_centroid_by_chunk[m*k*dim+j] = cluster_centroid[j];   
			}
			for ( n=0 ; n<taille; n++){		cluster_assignment_final[m*taille+n]=cluster_assignment_final_Y[n]+m*k; 
			}
			free(Y);
			Y = NULL;
			free(cluster_centroid); 
			cluster_centroid=NULL; 
			free(cluster_assignment_final_Y); 
			cluster_assignment_final_Y=NULL; 
		}
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
		printf ("group number %d\n", nb_groupes);
		free(cluster_centroid_by_chunk); 
		cluster_centroid_by_chunk = NULL;
		float count = k*N/taille; 

		get_cluster_member_count(((int) (N/taille)) * taille, (int) count, cluster_assignment_final, cluster_member_count);

		for (j=0; j<nb_groupes; j++){
			groupes[j].nb_members = 0; 
			for (n = 0; n< groupes[j].nb; n++){
			groupes[j].nb_members+=cluster_member_count[groupes[j].means[n]]; 
		}
		}

		double * chunk; 
		//X = getmatrix(source, dim, N,N,0, marks ); 

		chunk = form_chunk(groupes, /*X*/ source, marks,cluster_assignment_final, nb_groupes, N, dim, taille); 

		free(cluster_assignment_final); 
		cluster_assignment_final =NULL; 
		for (i = 0; i<nb_groupes; i++){
			free(groupes[i].means);
			free(groupes[i].members); 
		}
		free(groupes); 
		groupes = NULL;
		//printf ("cluster representatif formé\n"); 
		//cluster_centroid = random_center_init(chunk, taille, dim, k);

		gettimeofday(&tv_begin, NULL);
		cluster_centroid =kmeans_init_plusplus(chunk, taille, dim, k);
		gettimeofday(&tv_end, NULL);
		printf ("init time: %lf\n", (float)(tv_end.tv_sec - tv_begin.tv_sec)/60);
		
		//printf ("centres initialisés\n");

		cluster_assignment_final_Y=  (int *) malloc(taille*sizeof(int));
		kmeans(dim,chunk,taille, k, cluster_centroid, cluster_assignment_final_Y);
		//sort_centroid (dim, k, cluster_centroid);
		cluster_diag(dim, taille, k, chunk, cluster_assignment_final_Y, cluster_centroid);
		free(chunk);
		chunk = NULL;
		}
	else{
		printf ("lecture matrice\n"); 
		X = getmatrix(source, dim, N,N,0, marks );
		printf ("initialisation des centres\n"); 

		gettimeofday(&tv_begin, NULL);
		cluster_centroid =kmeans_init_plusplus(X, N, dim, k);
		gettimeofday(&tv_end, NULL);
		printf ("init time: %lf\n", (float)(tv_end.tv_sec - tv_begin.tv_sec)/60);
		//cluster_centroid = random_center_init(X, N, dim, k);
			
		cluster_assignment_final=  (int *) malloc(N*sizeof(int)); 
		printf ("kmeans\n");
		kmeans(dim,X,N, k, cluster_centroid, cluster_assignment_final);
		//sort_centroid (dim, k, cluster_centroid);
		cluster_diag(dim, N, k, X, cluster_assignment_final, cluster_centroid);
		free(X);
		X = NULL;
	}

	//compare_solution ( cluster_centroid, "a1-ga-cb.txt", dim, k); 
	free(var);
	var = NULL; 

	free(cluster_assignment_final_Y);
	cluster_assignment_final_Y = NULL; 
	free(cluster_centroid); 	
	cluster_centroid =NULL;
}
void kmeans_kmeans(double * X,int dim, int taille, int N, int k){ //expérimentation de la méthode de kmeans de kmeans sans formation de chunk représentatif
	int i, * cluster_assignment_final, *cluster_assignment_final_Y, m, j, n; 
	double * Y; 
	double * cluster_centroid=NULL, *cluster_centroid_by_chunk=NULL; 
	double *var; 
	int nb_groupes; 
	int cluster_member_count[MAX_CLUSTERS]; //à modifier pour une allocation dynamique 
	    
	    
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
		kmeans(dim,Y,taille, k, cluster_centroid, cluster_assignment_final_Y);
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
		
	kmeans(dim,cluster_centroid_by_chunk,k*N/taille, k, cluster_centroid, cluster_assignment_final_Y);

}

int main (int argc, char **argv){
	clock_t begin, end; 
	int i,j; 
	size_t  tot;
	char cmd[200];  
	sprintf (cmd, "echo %d > /cgroups/mem/kmeans/tasks", getpid()); 
	printf ("%s\n", cmd);
	system(cmd); 
	sleep(1);
	char * source = argv[1]; 
	size_t k = atoi(argv[2]);
	int taille_chunk, taille; 
	size_t dim, N; 
	double *X;
	struct timeval tv_begin, tv_end;
	N = atoi(argv[3]); 
	taille = atoi(argv[4]);
	dim = atoi(argv[5]); 
	taille_chunk = taille; 
	srand(time(NULL));
	//X= getmatrix(source, &dim, &tot,NULL,0 ); 
	//begin = clock();  
	gettimeofday(&tv_begin, NULL); 
	kmeans_by_chunk( source, dim, taille, N, atoi(argv[2])); 
	gettimeofday(&tv_end, NULL); 

	//kmeans_kmeans( X, dim, taille, tot, atoi(argv[2])); 

	//end = clock(); 
	printf ("%lf\n", (float)(tv_end.tv_sec - tv_begin.tv_sec)/60);
	return 0;
}
