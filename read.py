import numpy as np
import math
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import sys

def readpos(gridsize,n_particles_1d,a_initial,delta_a,snapshot_delta_a,total_steps,dt=1):
    f = open('pos.dat','r')

    plt.clf()
    fig = plt.figure(1)
    sp = fig.add_subplot(111)
    
    n_steps = total_steps/dt
    colors = iter(cm.rainbow(np.linspace(0,1,n_steps)))
    length = (n_particles_1d)**3 - n_particles_1d

    PI = 3.14159
    a_cross = 10.0 * a_initial
    nump = float(n_particles_1d)
    l_box = float(gridsize)
    k = 2 * PI / l_box
    amp = 1.0 / (a_cross * k)
    OMEGA_M = 1.0
    OMEGA_L = 0.0

    a_half = a_initial - (delta_a/ 2.0)
    D_dot =(math.sqrt(OMEGA_M + OMEGA_L * a_half**3.0) / math.sqrt(a_half));
   
    a = a_initial
    for i in range(7):
            line = f.readline()
            
    for t in range(total_steps):
        x = []
        v = []
        for i in range(n_particles_1d):
            line = f.readline()
            data = line.split()
            x.append(data[0])
            v.append(data[3])
        for i in range(length):
            line = f.readline()
       
        xtest = []
        vtest = []
        for i in range(n_particles_1d):
            checkx = (float(i)*(l_box/nump) + a*amp*math.sin(k*float(i)*(l_box/nump)))
            xtest.append(checkx)
            temp_a_half = a - (delta_a/2.0)
            D_dot =(math.sqrt(OMEGA_M + OMEGA_L * temp_a_half**3.0) / math.sqrt(temp_a_half));
            checkv = temp_a_half*temp_a_half*D_dot*amp*math.sin(k*float(i)*(l_box/nump))
            vtest.append(checkv)

        if (t%dt == 0):
            c = next(colors)
            sp.plot(x,v,'o',color=c)
            sp.plot(xtest,vtest,'-',color=c)
            print a,t

        a = a + snapshot_delta_a
           
#     plt.show()
    plt.savefig("ZA_pos.jpeg")
    f.close()

def readdens(gridsize,a_initial,delta_a,snapshot_delta_a,total_steps,dt):
    f = open('density.dat','r')

    plt.clf()
    fig = plt.figure(1)
    sp = fig.add_subplot(111)
    
    n_steps = total_steps/dt
    colors = iter(cm.rainbow(np.linspace(0,1,n_steps)))
    length = gridsize*gridsize - 1

    a = a_initial
    line = f.readline()
    
    for t in range(total_steps):
        x = []
        dens = []
        line = f.readline()
        data = line.split()
        
        for i in range(gridsize):
            x.append(i)
            dens.append(data[i])
        for i in range(length):
            line = f.readline()
        
        if (t%dt == 0):
            c = next(colors)
            sp.plot(x,dens,'o',color=c)
            sp.plot(x,dens,'-',color=c)
            print a,t
            
        a = a + snapshot_delta_a
       
    plt.show()
    f.close()
    


if __name__ == '__main__':
    gridsize = int(sys.argv[1])
    n_p_1d = int(sys.argv[2])
    a_init = float(sys.argv[3])
    delta_a = float(sys.argv[4])
    snapshot_delta_a = float(sys.argv[5])
    total_steps = int(sys.argv[6])
    dt_steps = int(sys.argv[7])
    readpos(gridsize,n_p_1d,a_init,delta_a,snapshot_delta_a,total_steps,dt_steps)
    #readdens(gridsize,a_init,delta_a,total_steps,dt_steps)
