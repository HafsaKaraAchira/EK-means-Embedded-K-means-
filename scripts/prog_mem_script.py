#!/usr/bin/env python3

import os
import math
import sys




N = int(sys.argv[1])
MC = int(sys.argv[2])
k = int(sys.argv[3])
dim = int(sys.argv[4])





def chunk_memory_est(N,dim,MC) :
    double_size = 8
    M = int(N/MC)
    chunk_size = M * dim * double_size
    
    return ( (chunk_size) / (1024**2) )


def prog_memory_est(N,dim,k,MC) :
    double_size = 8
    int_size = 4
    M = int(N/MC)
    chunk_size = M * dim * double_size
    dist_mat_size = M * k * double_size
    kmeans_assign_size = 2 * M *int_size
    chunk_assign_size = M * int_size
    dataset_assign_size = N * int_size
    
    size = chunk_size+dist_mat_size+kmeans_assign_size+chunk_assign_size+dataset_assign_size
    size = math.ceil( ( (size) / (1024**2) ) ) + 50
    
    return size




#################################################################
#################################################################
#################################################################


S=prog_memory_est(N,dim,k,MC)
os.system("scripts/prog_script_cgroup " + str(int(S)))
print("N/M="+str(MC)+"      CHUNK_SIZE="+str(chunk_memory_est(N,dim,MC))+"    MEM_SIZE="+str(int(S)))


#################################################################
#################################################################