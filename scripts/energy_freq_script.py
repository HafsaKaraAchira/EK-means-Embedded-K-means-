
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



#@measure_energy(handler=csv_handler, domains=[RaplCoreDomain(0),RaplDramDomain(0)])
def program_call():
    os.system("./program")

affinity_mask = {1}
pid = 0
os.sched_setaffinity(0, affinity_mask)

governor = "userspace"
available_memory = [100,250,500,750,1000]
unit = 100000

os.system("echo "+governor+" | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor")

for M in reversed(available_memory) :
    # $1 is the governor name
    os.system("echo $(( "+str(M)+" * 1024 * 1024)) > /sys/fs/cgroup/memory/prog_cgroup/memory.limit_in_bytes")
    print(M)
    for freq in range(29 * unit, 7 * unit, -unit) :
        os.system("./scripts/prog_script_cache")
        os.system("cpupower -c all frequency-set -f "+str(freq))
        os.system("cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_setspeed")
        time.sleep(30)
        with EnergyContext(handler=csv_handler, domains=[RaplCoreDomain(0),RaplDramDomain(0)], start_tag=governor+","+str(M)) as ctx:
            # call the target program
            program_call()
            
    csv_handler.save_data()



#################################################################