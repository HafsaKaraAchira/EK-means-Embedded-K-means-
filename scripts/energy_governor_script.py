
import os

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

governors = ["ondemand","conservative","schedutil","performance","powersave"];
available_memory = [100,200,300,400,500,600,700,800,900,1000]

for M in reversed(available_memory) :
    # $1 is the governor name
    os.system("echo $(( "+str(M)+" * 1024 * 1024)) > /sys/fs/cgroup/memory/prog_cgroup/memory.limit_in_bytes")
    print(M)
    for governor in governors :
        os.system("./scripts/prog_script_cache")
        os.system("echo "+governor+" | tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor")
        with EnergyContext(handler=csv_handler, domains=[RaplCoreDomain(0),RaplDramDomain(0)], start_tag=governor+","+str(M)) as ctx:
            # call the target program
            program_call()


csv_handler.save_data()

#################################################################