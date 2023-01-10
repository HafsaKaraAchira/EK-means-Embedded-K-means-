
import os
import time

from pyJoules.device.rapl_device import RaplCoreDomain , RaplDramDomain
from pyJoules.handler.csv_handler import CSVHandler
# from pyJoules.handler.print_handler import PrintHandler

from pyJoules.energy_meter import EnergyContext
# from pyJoules.energy_meter import EnergyMeter
# from pyJoules.device import DeviceFactory
# from pyJoules.energy_meter import measure_energy
	
csv_handler = CSVHandler('./logs/log_energy.csv')

file = "generator/CM5,8M_1000MO_SEP0,3/points.csv"
N = 5800000
k = 10
dim = 10

#@measure_energy(handler=csv_handler, domains=[RaplCoreDomain(0),RaplDramDomain(0)])
def program_call(M):
    os.system("./program "+file+" "+str(k)+" "+str(N)+" "+str(int(N/M))+" "+str(dim))

affinity_mask = {1}
pid = 0
os.sched_setaffinity(0, affinity_mask)

governors = ["ondemand"]  #["performance","conservative","schedutil","performance","powersave"];
memory_constraint = [1] #0.8,1,2,4,5,8,10]

for M in memory_constraint :
    # $1 is the governor name
    os.system("echo $(( "+str(int(1000/M))+" * 1024 * 1024)) > /sys/fs/cgroup/memory/kmeans/memory.limit_in_bytes")
    print(str(M)+"  "+str(int(1000/M))+"    "+str(int(N/M)))
    for governor in governors :
        os.system("echo "+governor+" | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor")
        os.system("./scripts/prog_script_cache")
        time.sleep(30)
        with EnergyContext(handler=csv_handler, domains=[RaplCoreDomain(0),RaplDramDomain(0)], start_tag=governor+","+str(M)) as ctx:
            # call the target program
            program_call(M)


csv_handler.save_data()

#################################################################