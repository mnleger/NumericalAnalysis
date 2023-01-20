#Matthew Noah Leger
#Final Project
 
import numpy as np
import math

# constants
k = 1/16
g = 9.81
loopcount = 11

#Question 8
 
#Question a
#function inspired by Dr. Lovekin's example
def rkSystem(a,b,m,N,alphas,functions):
    #set up step size
    h = (b - a) / N
 
    #set up t and w arrays
    t = np.linspace(a,b, N+1)
    w = np.zeros((m,len(t)))
 
    #initialize k arrays
    k1 = np.zeros(m)
    k2 = np.zeros(m)
    k3 = np.zeros(m)
    k4 = np.zeros(m)
 
    for i in range(0,m):
        #set the initial solutions to be the initial conditions
        w[i,0] = alphas[i]
    for i in range(1, N+1):
        for j in range(0,m):
            #calculate k1 for each function
            f = functions[j]
            k1[j] = h * f(t[i-1], w[:,i-1])
        for j in range(0,m):
            #calculate k2 for each function
            f = functions[j]
            k2[j] = h * f(t[i-1] + h/2, w[:,i-1] + k1/2)
        for j in range(0,m):
            #calculate k3 for each function
            f = functions[j]
            k3[j] = h * f(t[i-1] + h/3, w[:,i-1] + k2/2)
        for j in range(0,m):
            #calculate k4 for each function
            f = functions[j]
            k4[j] = h * f(t[i-1] + h, w[:,i-1] + k3)
        for j in range(0,m):
            #Calculate the solution at i + 1 for each j
            w[j,i] = w[j,i-1] + (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6
    return (t,w)
 
#rate of change of position (velocity)
def dpdt(t,I):
    return math.pi * I[2]**2 * I[0] / I[1] * drdt(t,I) * I[0] / I[1] + I[1] * g
 
#rate of change of mass
def dmdt(t,I):
    return math.pi * I[2]**2 * I[0] / I[1] * drdt(t,I)
 
#rate of change of radius (assume k = 1/16 for simplicity's sake)
def drdt(t,I):
    return k
 
#assume r(0) = a = 0.00001
t,w = rkSystem(0,10, 3, 10, [0,(math.pi * 4 * 0.00001**3)/3,k], [dpdt,dmdt,drdt])
position = w[0]
mass = w[1]
radius = w[2]
 
vel = []
 
#calculate velocity using the position values previously calculated
def velocity(p,m,r):
    return math.pi * r**2 * p**2 / m**2 * k + m * g
 
for i in range(0,loopcount):
    vel.append(velocity(position[i], mass[i], radius[i]))
 
#print the results in tabluar form
print("Tabluar form of the data in the format (time, position, velocity, acceleration)")

#built in python library used to format the data in a tabular form
#print(0, position[0], vel[0], 0)
for i in range(0,loopcount):
    print("{:>5,.0f}".format(i), "{:>12,.5f}".format(position[i]), "{:>12,.5f}".format(vel[i]), "{:>12,.5f}".format(g))
 
#Question b
print()
#calculate the rate of change of volume (which is equivalent to the rate of change of mass, after doing a unit conversion)
def dvdt(p,m,r):
    return math.pi * r**2 * p / m * k
 
vol = []
 
for i in range(0,loopcount):
    vol.append(dvdt(position[i], mass[i], radius[i]))
 
#output the rate of change of volume
print("Rate of change of volume (time, volume):")
for i in range(0, loopcount):
    print("{:>5,.0f}".format(i), "{:>12,.5f}".format(vol[i]))
#print(vol)
