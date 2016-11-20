# -*- coding: utf-8 -*-
"""
Created on Sat Jun 04 21:56:49 2016

@author: Sergey
"""


import scipy as sc
import matplotlib.pylab as plt
from scipy.integrate import ode
import numpy

lambd = 2.0#parameter of the problem
n = 500#number of points that we obtain after each integration
TT = 4#time segment

#We define the function for integration and the jacobian(that we won't use here)
def f(t, y):
    return -lambd * numpy.sqrt(1 + 0.5 * numpy.sin(y))
def jac(t, y):
    return -0.5 * lambd * (numpy.cos(y))/numpy.sqrt(4 + 2 * numpy.sin(y))

#This code implements integration on the segment [t0, t1] for the initial value y0.
#This code is created according the example from the documentation of scipy.integrate.
#k is a number of points, that we are going to obtain after the integration.
#Attention!! the number of points of integration is set by the integrator 'isoda' and it's done automatically
#so k is only the number of points that we are going to have.
def integrate(y0, t0, t1, k):
    loc_list = []
    r = ode(f, jac).set_integrator('lsoda', with_jacobian=False)
    r.set_initial_value(y0, t0)
    step = (t1 - t0)/k
    while r.successful() and r.t <= t1:
        r.integrate(r.t + step)
        loc_list = loc_list + [[r.t, r.y]]
    return loc_list
 
#After having obtained all points of the trajectory of the process
#we need to calculate the average. 
#This code returns values for a given time t   
def search_in_array(P,t):
    for i in range(len(P)):
        if (P[i][0] >= t):
            return P[i][1]
    
#Arrays for check the hypothethis. A - for the moments of time: 0, 0.4, 0.8...
#B is for the average value in these points
A = [i * 0.4 for i in range(11)]
B = [0 for i in range(11)]

#Number of trajectories
num = 100

#The main loop
for i in range(num):
    
    L = []#total list of all points of the trajectory
    TR = 0.0
    #times - is list of the moments of jumps of the Poisson process
    times = [0.0]
    
    #Here we modelise the Poisson process as in S1
    while(TR < TT):
        TR = TR + numpy.random.exponential(1/lambd, 1)
        times = times + [TR]
        
    y_curr = 1#initial value
    L = L + [[0, y_curr]]#adds initial point
    
    #this loop builds the trajectory
    for j in range(len(times) - 1):
        c = integrate(y_curr, times[j], times[j + 1], n)#we integrate between two neghbouring points
        y_curr = c[len(c) - 1][1]# we upload the current value
        y_curr = y_curr + numpy.sqrt(1 + 0.5 * numpy.sin(y_curr))#we considering the next jump
        L = L + c[:len(c) - 1] + [[times[j + 1], y_curr]]#we add obtained points to the total list also considering the jump
        
    #we  add the values of the trajectory to calculate the average 
    for j in range(len(B)):
        B[j] = B[j] + search_in_array(L, A[j])
        
    #this code returns the graphs for the first 4 trajectories  
    
    if (i < 4):
        plt.clf()
        plt.cla()
        plt.xlabel('t', fontsize = 20)
        plt.ylabel('X',fontsize = 20)
        plt.title('X(t)',fontsize = 20)
        plt.plot([L[i][0] for i in xrange(len(L))], [L[i][1] for i in xrange(len(L))])
        s = 'graph' + str(i) + '.png'#or .eps
        plt.savefig(s)
    



#calculates the average for each point
B = [B[i]/num for i in range(len(B))]


#print(B)
#returns the plot of the average
t = numpy.arange(0., 4., 0.2)
plt.clf()
plt.cla()
plt.xlabel('t', fontsize = 20)
plt.ylabel('Moyenne(t)',fontsize = 20)
plt.title('Moyenne pour 100 trajectoires',fontsize = 20)
plt.plot(A, B, 'bs', t, [1.0 for i in range(len(t))],'r--')
#plt.savefig('Moyenne.eps')   















