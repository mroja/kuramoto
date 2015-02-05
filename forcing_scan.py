
import os
import sys

N_steps = 1500
dt = 0.01
noise = 0.0
coupling_type = 'kuramoto_2'
dump_interval = 30

name_i = 10
for ff in [1.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0]:
    for fs in [0.1, 1.0, 5.0, 10.0, 20.0, 40.0]:
        name = str(name_i)

        directory = 'dump_' + name
        if not os.path.exists(directory):
            os.makedirs(directory)

        #forcing_strength = fs # 20.5
        #forcing_freq = ff #30.0
        #os.system('python gen_presets.py ' + name)
        #conf = "%i\n%f\n%f\n%i\n%s\n%f\n%f" % (N_steps, dt, noise, dump_interval, coupling_type, forcing_strength, forcing_freq)
        #with open(name + '.conf.txt', 'w') as f:
        #    f.write(conf)
        #os.system('kuramoto_simulation.exe ' + name)

        os.system('python gen_plots.py ' + name)
        
        name_i += 1
