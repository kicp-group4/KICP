import numpy as np
import math
from matplotlib import pyplot as plt
import matplotlib.cm as cm

def readpos():
    f = open('pos.dat','r')

    plt.clf()
    fig = plt.figure(1)
    sp = fig.add_subplot(111)

    #testing123
    
    colors = iter(cm.rainbow(np.linspace(0,1,10)))
    length = 1000 - 8
    length = 32*32*32 - 30

    DELTA_A  = 0.0001
    PI = 3.14159
    H_0 = 72
    A_INITIAL = 0.2
    a_cross = 10.0 * A_INITIAL
    Np = 32
    L_BOX = float((10.0*3.086e16))
    L_NP = L_BOX / (Np-1)
    k = 2 * PI / L_BOX
    D_aini = 1
    amp = -1.0 / (a_cross * k / A_INITIAL)
    OMEGA_M = 0.23
    OMEGA_L = 0.77
    r_0 = (10.0*3.086e16)/32.0

    a_temp = A_INITIAL - (DELTA_A / 2.0)
    D_dot = (H_0 * 3.241e-20)*(math.sqrt(OMEGA_M + OMEGA_L * a_temp**3.0) / math.sqrt(a_temp)) / A_INITIAL;
   
    temp_a = A_INITIAL
    for dt in range(1):
        x = []
        v = []
        for i in range(7):
            a = f.readline()
            print a
        for i in range(32):
            a = f.readline()
            data = a.split()
            x.append(data[0])
            v.append(data[3])
        for i in range(length):
            a = f.readline()

        xtest = []
        vtest = []
        for i in range(32):
            checkx = (float(i)*L_NP + (temp_a/A_INITIAL)*amp*math.sin(k*float(i)*L_NP))/(3.086e16)
            xtest.append(checkx)
            checkv = temp_a*D_dot*amp*math.sin(k*float(i)*L_NP)
            vtest.append(checkv)
        
        lx = 9
        lv = max(v)
        
        c = next(colors)
        sp.plot(x,v,'o',color=c)
        sp.plot(xtest,vtest,'-',color=c)
        a = "Step #"+str(dt)
        plt.text(lx,lv,a,color=c)
        temp_a = temp_a + DELTA_A
        
    plt.show()
    f.close()

    xg = []
    f = open('density.dat','r')
    a = f.readline()
    a = f.readline()
    data = a.split()
    dens = []
    for i in range(32):
        step = (10.0/32.0)*float(i)
        xg.append(step)
        dens.append(data[i])

    print xg
    print x
    
    f.close()

def readdens():
    f = open('density.dat','r')

    plt.clf()
    fig = plt.figure(1)
    sp = fig.add_subplot(111)

    #testing123
    
    colors = iter(cm.rainbow(np.linspace(0,1,10)))
    length = 1000 - 1
    length = 32*32*32 - 30

    for dt in range(1):
        x = []
        a = f.readline()
        a = f.readline()
        data = a.split()
        dens = []
        for i in range(32):
            step = (10.0/32.0)*float(i)
            x.append(step)
            dens.append(data[i])
        
        lx = 9
        lv = max(dens)
        
        c = next(colors)
        sp.plot(x,dens,'o',color=c)
       
    plt.show()
    f.close()
    
