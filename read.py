import numpy as np
import math
from matplotlib import pyplot as plt
import matplotlib.cm as cm

def readpos():
    f = open('pos.dat','r')

    plt.clf()
    fig = plt.figure(1)
    sp = fig.add_subplot(111)
    
    colors = iter(cm.rainbow(np.linspace(0,1,20)))
    length = 1000 - 8
    length = 64*64*64 - 62

<<<<<<< HEAD
    DELTA_A  = 0.02
=======
    DELTA_A  = 0.01
>>>>>>> 27cdbc009f60e758e06aba7ecf990cc4003899a3
    PI = 3.14159
    H_0 = 72
    A_INITIAL = 0.2
    a_cross = 10.0 * A_INITIAL
    Np = 32
<<<<<<< HEAD
    L_BOX = 64.0
    k = 2 * PI / 64.0
=======
    L_BOX = 32.0
    k = 2 * PI / 32.0
>>>>>>> 27cdbc009f60e758e06aba7ecf990cc4003899a3
    D_aini = 1
    amp = 1.0 / (a_cross * k)
    OMEGA_M = 1.0
    OMEGA_L = 0.0
    r_0 = (10.0*3.086e16)/64.0

    a_temp = A_INITIAL - (DELTA_A / 2.0)
    D_dot =(math.sqrt(OMEGA_M + OMEGA_L * a_temp**3.0) / math.sqrt(a_temp));
   
    temp_a = A_INITIAL
    for dt in range(200):
        x = []
        v = []
        for i in range(7):
            a = f.readline()
        for i in range(64):
            a = f.readline()
            data = a.split()
            x.append(data[0])
            v.append(data[3])
        for i in range(length):
            a = f.readline()
       
        xtest = []
        vtest = []
        for i in range(64):
            checkx = (float(i) + temp_a*amp*math.sin(k*float(i)))
            xtest.append(checkx)
            temp_a_half = temp_a - (DELTA_A/2.0)
            D_dot =(math.sqrt(OMEGA_M + OMEGA_L * temp_a_half**3.0) / math.sqrt(temp_a_half));
            checkv = temp_a_half*temp_a_half*D_dot*amp*math.sin(k*float(i))
            vtest.append(checkv)
        

        lx = 63
        lv = max(v)
        
        if (dt%10 == 0):
            c = next(colors)
            sp.plot(x,v,'o',color=c)
            sp.plot(xtest,vtest,'-',color=c)
            print temp_a,dt
            a = "a = "+str(temp_a)
            plt.text(lx,lv,a,color=c)
        temp_a = temp_a + DELTA_A
        #print temp_a,dt

        
    #plt.show()
    plt.savefig("ZA_pos.jpeg")
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
    
    f.close()

def readdens():
    f = open('density.dat','r')

    plt.clf()
    fig = plt.figure(1)
    sp = fig.add_subplot(111)

    #testing123
    
    colors = iter(cm.rainbow(np.linspace(0,1,20)))
    length = 1000 - 1
    length = 32*32*32 - 30

    for dt in range(20):
        x = []
        a = f.readline()
        a = f.readline()
        data = a.split()
        dens = []
        for i in range(32):
            x.append(i)
            dens.append(data[i])
        for i in range(31*32):
            a = f.readline()
        
        lx = 9
        lv = max(dens)
        
        c = next(colors)
        sp.plot(x,dens,'o',color=c)
        sp.plot(x,dens,'-',color=c)
       
    plt.show()
    f.close()
    
def checkdensity():
    f = open('density.dat','r')

    length = 32*32*32 - 30

    for dt in range(1):
        x = []
        a = f.readline()
       
        dens = []
        total = 0.0
        for i in range(32):
            for j in range(32):
                 a = f.readline()
                 print a
                 data = a.split()
                 for k in range(32):
                    total = total + float(data[k])

            print total
            
        print 32*32*32
        f.close()
    
def testgrav():
    x = [0, 	1.09936 ,	2.1949, 	3.28295 ,	4.36013 ,	5.42346 ,	6.47053 ,	7.49951 ,	8.5093 ,	9.49951 ,	10.4705 ,	11.4235 ,	12.3601 ,	13.283, 	14.1949 ,	15.0994 ,	16 ,	16.9006 ,	17.8051 ,	18.7171 ,	19.6399 ,	20.5765 ,	21.5295 ,	22.5005 ,	23.4907 ,	24.5005 ,	25.5295 ,	26.5765 ,	27.6399 ,	28.717, 	29.8051 ,	30.9006]

    g = [-3.16758e-07 ,	0.445213 ,	1.23765 ,	1.96725 ,	2.6013 ,	3.11378, 	3.48665 ,	3.71026 ,	3.78302 ,	3.71041 ,	3.50347 ,	3.17716 ,	2.74871 ,	2.23623 ,	1.65776 ,	0.993836 ,	1.32012e-05, 	-0.993815 ,	-1.65774 ,	-2.23621 ,	-2.7487, 	-3.17715 ,	-3.50346 ,	-3.7104 ,	-3.78302 ,	-3.71027 ,	-3.48666 ,	-3.11379 ,	-2.60131 ,	-1.96727 ,	-1.23767 ,	-0.445227 ]

    plt.clf()
    fig = plt.figure(1)
    sp = fig.add_subplot(111)

    grav = []
    for i in range(32):
        grav.append(((3/(2 * 0.2))*(x[i] - float(i)))) 

    sp.scatter(x,g)
    sp.plot(x,grav)
    plt.show()

    
