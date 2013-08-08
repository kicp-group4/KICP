import numpy as np
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
    
    for dt in range(5):
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
            
        lx = 9
        lv = max(v)
        
        c = next(colors)
        sp.plot(x,v,'o',color=c)
        a = "Step #"+str(dt)
        plt.text(lx,lv,a,color=c)

    plt.show()
        
   
    
    f.close()
    
