import subprocess


# subprocess.run('pwd')
subprocess.run(['mpirun', '-np', '2', '/home/zhou/anaconda3/envs/ffip/bin/run_sim_json'])