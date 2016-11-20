# -*- coding: utf-8 -*-


import scipy as sc
import numpy
import matplotlib.pylab as plt

#time segment [0, TT]
TT = 3.0
#list of jumps
times = [0.0]
#list of values
values = [0.0]
#current time
TR = 0.0
#current value
N = 0.0
#parameter of the process
lambd = 4.0

#we modelise the intervals of jumps as the values having the exponential distribution
#we do it until the leaving our segment [0, TT]
while(TR < TT):
    TR = TR + numpy.random.exponential(1/lambd, 1)
    times = times + [TR]# we add the time of the next jump
    N = N + 1.0#we increase the value of process
    values = values + [N]#we add the value corresponding to time of jump
    
#we return the plot
plt.xlabel('t', fontsize = 20)
plt.ylabel('N',fontsize = 20)
plt.title('Processus de Poisson', fontsize = 20)
plt.step(times, values, where='post')
plt.savefig('poisson.png')#or 'poisson.png' 
        
