import subprocess


# subprocess.run('pwd')
subprocess.run(['mpirun', '-np', '2', 'run_sim_json'])