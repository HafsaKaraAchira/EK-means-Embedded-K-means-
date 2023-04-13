
import os
import time
import math

from pyJoules.device.rapl_device import RaplCoreDomain , RaplDramDomain
from pyJoules.handler.csv_handler import CSVHandler
from pyJoules.energy_meter import EnergyContext

csv_handler = CSVHandler('./logs/log_energy_io.csv')

file = "generator/10D/13421800N/SEP-0.9/points.csv"
N = 13421800
k = 10
dim = 10

#@measure_energy(handler=csv_handler, domains=[RaplCoreDomain(0),RaplDramDomain(0)])
def program_call(MC):
    os.system("./prog_energy_ref "+file+" "+str(k)+" "+str(N)+" "+str(int(N/MC))+" "+str(dim))

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
    size = math.ceil(( (size) / (1024**2)) + 1.5 ) +1
    
    return size


affinity_mask = {1}
pid = 0
os.sched_setaffinity(0, affinity_mask)

governor = "ondemand"
memory_constraint = [1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8]

os.system("echo "+governor+" | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor")

S_dataset=prog_memory_est(N,dim,k,1)    

for MC in memory_constraint :
    S=prog_memory_est(N,dim,k,MC)
    os.system("echo $(( "+str(int(S))+" * 1024 * 1024)) > /sys/fs/cgroup/memory/kmeans/memory.limit_in_bytes")
    print(str(MC)+"  "+str(int(1024/MC))+"    "+str(S))
    os.system("./scripts/prog_script_cache")
    time.sleep(5)
    with EnergyContext(handler=csv_handler, domains=[RaplCoreDomain(0),RaplDramDomain(0)], start_tag=str(N)+","+str(int(N))+","+str(MC)+","+str(S)) as ctx:
        # call the target program
        program_call(MC=1)
    csv_handler.save_data()



#################################################################