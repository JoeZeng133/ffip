import ffip
import subprocess

# subprocess.run(['mpirun', '--version'])
# subprocess.run(['mpirun', '-np', '2', 'run_sim_json'], check=True, shell=True)
process = subprocess.Popen('mpirun -np 2 run_sim_json', shell=True)
process.wait()