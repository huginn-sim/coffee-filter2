# -*- coding: utf-8 -*-
"""
.. module:: coffee-filter
   :synopsis: Models the position, velocity, and acceleration of a falling coffee filter.

.. moduleauthor:: Huginn
"""

#~ Modules
from pylab import *

#~ ODE Solvers
def euler(f,t,dt,x):
    return x + f(t,x)*dt

def euler_richardson(f,t,dt,x):
    return x + f(t + dt/2., x + f(t,x)*dt/2.)*dt

def rk4(f,t,dt,x):
    k1 = f(t        , x        )*dt
    k2 = f(t + dt/2., x + k1/2.)*dt
    k3 = f(t + dt/2., x + k2/2.)*dt
    k4 = f(t + dt   , x + k3   )*dt

    return x + (1./6.)*(k1 + 2*k2 + 2*k3 + k4)

def predict_correct(f,t,dt,x):
    t[0:0] = [t[0] - dt]
    # Roll back by one time-step.
    pc_state = [rk4(f,t[0],-dt,x[0])]
    pc_state.append(x[0])
    # Roll forward.
    for t in times:
        xp = pc_state[-2] + f(t,pc_state[-1])*2*dt
        xc = pc_state[-1] + 0.5*(f(t,xp) + f(t,pc_state[-1]))*dt
        pc_state.append(xc)

    return times, array(pc_state)
#/~ ODE Solvers

#~ ODEs
def fall(t,x):
    return np.array([x,g])
#/~ ODEs

#-> Not quite sure this is 100% correct...
def finite_diff(f,times,dt,x):
    v_fd = lambda f,t,dt,x: (x+f(t,x) - x+f(t-dt,x))/(2.*dt)
    a_fd = lambda f,t,dt,x: (x+f(t,x) - 2*x + x+f(t-dt,x))/(dt**2)

    print len(times)
    v = [v_fd(f,times[0],dt,x[0])]
    a = [a_fd(f,times[0],dt,x[0])]
    for i in range(1, len(times)):
        t = times[i]
        dt = t - times[i-1] 
        v.append(v_fd(f,t,dt,x[i]))
        a.append(a_fd(f,t,dt,x[i]))

    v = array(v)
    a = array(a)
    return v[:,0], a[:,0]

def sample_data(data=None):
    """ Calculates thermal constants from a dataset.
        
        :param data: The dataset used to derive cooling constants. If 'None' use default 'data'.
    """

    if data == None:    # Use default 'data' if none is provided by user.
        data = np.array([[  .2055,.2302,.2550,.2797,.3045,.3292,.3539,.3786,.4033,
                            .4280,.4526,.4773,.5020,.5266,.5513,.5759,.6005,.6252,
                            .6498,.6744,.6990,.7236,.7482,.7728,.7974,.8220,.8466],
                         [  .4188,.4164,.4128,.4082,.4026,.3958,.3878,.3802,.3708,
                            .3609,.3505,.3400,.3297,.3181,.3051,.2913,.2788,.2667,
                            .2497,.2337, .2175,.2008,.1846,.1696,.1566,.1393,.1263]])

    times = data[0,:]; print times
    pos = data[1,:]; print pos
    vel, acc  = finite_diff(fall,times,dt,pos); print vel; print acc
    
    plt.plot(acc, vel, 'b.', lw=5)
    plt.show()

    #~> Fix this
    return None

#~ Global Variables
t0 = 0; tf = 10*pi; dt = .1
times = np.arange(t0, tf, dt)

g = 9.81
x0 = 1.; v0 = 0.; a0 = g
state = [np.array([x0,v0,a0])]
#/~ Global Variables

#~ Entry point of the script.
if __name__ == "__main__":
    sample_data()
    
#sliding average for smoothing
def slide_avg(x):
    avg = []
    for i in range(len(x)-2):
        avg.append((x[i]+x[i+1]+x[i+2])/3)
    return np.array(avg)
position_avg = slide_avg(position)
#can change the last line to velocity and then acc
#This is what I have been working on so far
t = time    
x = position
for i in range(1,len(x)-1):
    dt = ((t[i+1] + t[i-1])/2)
    v.append((x[i+1] - x[i-1])/(2*(dt)))
    a.append((x[i+1] - 2*(x[i]) + x[i-1])/(dt**2))

vT = -0.025

def falld1(t,x):
    return np.array([x[1],g(1-(x[1]/vT))])

def falld2(t,x):
    return np.array([x[1],g(1-(x[1]/vT)**2)])
# this is when I was going to run rk with the two falling models

