
import os
import time
import math

from pyJoules.device.rapl_device import RaplCoreDomain , RaplDramDomain
from pyJoules.handler.csv_handler import CSVHandler
# from pyJoules.handler.print_handler import PrintHandler

from pyJoules.energy_meter import EnergyContext
# from pyJoules.energy_meter import EnergyMeter
# from pyJoules.device import DeviceFactory
# from pyJoules.energy_meter import measure_energy
	
csv_handler = CSVHandler('./logs/log_energy_gov.csv')

file = "generator/10D/13421800N/SEP-0.9/points.csv"
N = 13421800
k = 10
dim = 10

#@measure_energy(handler=csv_handler, domains=[RaplCoreDomain(0),RaplDramDomain(0)])
def program_call(M):
    os.system("./prog_energy_ref "+file+" "+str(k)+" "+str(N)+" "+str(int(N/M))+" "+str(dim))

affinity_mask = {1}
pid = 0
os.sched_setaffinity(0, affinity_mask)

governors = ["ondemand","performance","conservative","schedutil"] #,"powersave"]
memory_constraint = [1,2,4,5,8,10]

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
    
    

for MC in reversed(memory_constraint) :
    # $1 is the governor name
    S=prog_memory_est(N,dim,k,MC)
    os.system("echo $(( "+str(int(S))+" * 1024 * 1024)) > /sys/fs/cgroup/memory/kmeans/memory.limit_in_bytes")
    print(str(MC)+"  "+str(int(1024/MC))+"    "+str(S))
    for governor in governors :
        os.system("echo "+governor+" | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor")
        os.system("./scripts/prog_script_cache")
        time.sleep(5)
        with EnergyContext(handler=csv_handler, domains=[RaplCoreDomain(0),RaplDramDomain(0)], start_tag=str(N)+","+str(int(N/MC))+","+str(MC)+","+governor) as ctx:
            # call the target program
            program_call(MC)
            
        csv_handler.save_data()


#################################################################