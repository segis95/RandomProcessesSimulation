# -*- coding: utf-8 -*-
"""
Created on Sun May 29 00:16:06 2016

@author: Sergey
"""
import scipy as sc
import numpy
import matplotlib.pylab as plt

#Segment of time
TT = 3
#List of the jumps of the process of Poisson
L = [0.0]
#parameters of the problem
lambd = 3.0
sigma = 2.0
#current time
TR = 0.0
#initial value
x = 1.0



#As we can integrate here explicitely in this problem, we do it!
#This code makes this "integration" and inserts the points obtained between
#the moments of jumps(i.e. we EXPAND our array)
def expand(T):
    S = []#array for time
    k = 1000#number of point demanded between two jumps
    Y = []#array for values
    
    
    for i in range(len(T) - 1):
        P = [0 for q in range(k)]#local array of times between to jumps
        P[0] = T[i]#giving the initial value
        tau = (T[i + 1] - T[i])/k#step of "integration"
        #adding values to P
        for j in range(1,k):
            P[j] = P[0] + tau * j
        S = S + P#expanding the total array of time
    
    #We have already expanded the array of times
    #So now we expand the array of values by the formula obtained in T8
    for i in range(len(S)):
            Y = Y + [x * numpy.exp(-lambd*sigma*S[i])*(sigma + 1)**(i/k)]

    return (S,Y)
    

#********************************************
#Start of the implementing part

#we simulate the Poisson process as in S1
while(TR < TT):
    TR = TR + numpy.random.exponential(1/lambd, 1)
    L = L + [TR]
    
#We expand the arrays
X,Y = expand(L)

#We save the plot
plt.plot(X,Y)
plt.xlabel('t', fontsize = 20)
plt.ylabel('X',fontsize = 20)
plt.title('X(t)',fontsize = 20)

plt.savefig('graph.eps') #or 'graph.png'
