# -*- coding: utf-8 -*-
"""
.. module:: coffee-filter
   :synopsis: Models the position, velocity, and acceleration of a falling coffee filter.

.. moduleauthor:: Huginn
"""

#~ Modules
from pylab import *
import  sys
#/~ Modules

#~ Custom Modules
sys.path.append("C:\Users\Evin\Documents\Github")
from viz.display.plot import configure
#/~ Custom Modules

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

def finite_diff(times,x):
    v_fd = lambda dt,x: (x[0] - x[1])/(2.*dt)
    a_fd = lambda dt,x: (x[0] - 2*x[1] + x[2])/(dt**2)

    v = []; a = []; dts = []
    for i in range(1, len(times)-1):
        dts.append((times[i+1] - times[i-1]) /2.)
        v.append(v_fd(dts[i-1],[x[i+1], x[i-1]]))
        a.append(a_fd(dts[i-1],[x[i+1], x[i], x[i-1]]))

    v = array(v); a = array(a); dts = array(dts)
    return v, a, dts

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
                            .2497,.2337,.2175,.2008,.1846,.1696,.1566,.1393,.1263]])

    times = data[0,:]; print times
    pos = data[1,:]; print pos
    vel, acc, dts = finite_diff(times,pos); print vel; print acc
    
    return times, dts, pos, vel, acc

def plot_samples(times, pos, vel, acc):
    fig, axes = subplots(2, 2)

    ptmark, = axes[0,0].plot(times, pos, 'b--o')
    start_ptmark, = axes[0,0].plot(times[0], pos[0], 'go', ms=10)
    end_ptmark, = axes[0,0].plot(times[-1], pos[-1], 'ro', ms=10)
    
    axes[0,0].legend(  [ptmark, start_ptmark, end_ptmark],
                [r'$\left(t,p\right)$', r'First', r'Last'],
                numpoints=1,
                loc="upper right")

    configure(  ax=axes[0,0],
                title="Position vs Time",
                xlabel=r"Time $\left(seconds\right)$",
                ylabel=r"Position $\left(meters\right)$",
                xbounds=None, ybounds=None)

    vtmark, = axes[1,0].plot(times[1:-1], vel, 'b--o')
    start_vtmark, = axes[1,0].plot(times[1], vel[0], 'go', ms=10)
    end_vtmark, = axes[1,0].plot(times[-2], vel[-1], 'ro', ms=10)

    axes[1,0].legend(  [vtmark, start_vtmark, end_vtmark],
                [r'$\left(t,v\right)$', r'First', r'Last'],
                numpoints=1,
                loc="upper right")
    
    configure(  ax=axes[1,0],
                title="Velocity vs Time",
                xlabel=r"Time $\left(seconds\right)$",
                ylabel=r"Velocity $\left(\frac{meters}{second}\right)$",
                xbounds=None, ybounds=None)

    atmark, = axes[0,1].plot(times[1:-1], acc, 'b--o')
    start_atmark, = axes[0,1].plot(times[1], acc[0], 'go', ms=10)
    end_atmark, = axes[0,1].plot(times[-2], acc[-1], 'ro', ms=10)

    axes[0,1].legend(  [atmark, start_atmark, end_atmark],
                [r'$\left(t,a\right)$', r'First', r'Last'],
                numpoints=1,
                loc="lower right")
    
    configure(  ax=axes[0,1],
                title="Acceleration vs Time",
                xlabel=r"Time $\left(seconds\right)$",
                ylabel=r"Acceleration $\left(\frac{meters}{seconds^2}\right)$",
                xbounds=None, ybounds=(-10,10))

    avmark, = axes[1,1].plot(vel, acc, 'b--o')
    start_avmark, = axes[1,1].plot(vel[0], acc[0], 'go', ms=10)
    end_avmark, = axes[1,1].plot(vel[-1], acc[-1], 'ro', ms=10)

    axes[1,1].legend(  [avmark, start_avmark, end_avmark],
                [r'$\left(v,a\right)$', r'First', r'Last'],
                numpoints=1,
                loc="lower right")
    
    configure(  ax=axes[1,1],
                title="Acceleration vs Velocity",
                xlabel=r"Velocity $\left(\frac{meters}{second}\right)$",
                ylabel=r"Acceleration $\left(\frac{meters}{seconds^2}\right)$",
                xbounds=None, ybounds=(-10,10))

    fig.suptitle("Falling Coffee Filter", size=30)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.08)

    show()

def plot_acc_vs_vel(vel,acc):
    fig, ax = plt.subplots()

    avmark, = ax.plot(vel, acc, 'b--o')
    start_avmark, = ax.plot(vel[0], acc[0], 'go', ms=10)
    end_avmark, = ax.plot(vel[-1], acc[-1], 'ro', ms=10)

    ax.legend(  [avmark, start_avmark, end_avmark],
                [r'$\left(v,a\right)$', r'First', r'Last'],
                numpoints=1,
                loc="lower right")
    
    configure(  ax=ax,
                title="Acceleration vs Velocity",
                xlabel=r"Velocity $\left(\frac{meters}{second}\right)$",
                ylabel=r"Acceleration $\left(\frac{meters}{seconds^2}\right)$",
                xbounds=None, ybounds=(-10,10))

    fig.suptitle("Falling Coffee Filter", size=30)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.08)

    show()

def plot_samples_and_models(times, pos, vel, acc):
    F1 = [array([pos[0], vel[0]])]
    F2 = [array([pos[0], vel[0]])]
    for t,dt,p,v in zip(times,dts,pos,vel)[1:]:
        F1.append(rk4(falld1, t, dt, F1[-1]))
        F2.append(rk4(falld2, t, dt, F2[-1]))

    F1 = array(F1)
    F2 = array(F2)

    fig, axes = subplots(2)

    pmark, = axes[0].plot(times, pos, 'b--o')
    plmark, = axes[0].plot(times[:-2], F1[:,0], 'g--o')
    pqmark, = axes[0].plot(times[:-2], F2[:,0], 'r--o')

    axes[0].legend(  [pmark, plmark, pqmark],
                [r'$\left(t,p\right)$', r'Linear', r'Quadratic'],
                numpoints=1,
                loc="upper right")
    
    configure(  ax=axes[0],
                title=r'Position vs Time $\left(v_t='+ str(vT)+ r'\right)$',
                xlabel=r"Time $\left(seconds\right)$",
                ylabel=r"Position $\left(meters\right)$",
                xbounds=None, ybounds=(0,.5))

    vmark, = axes[1].plot(times[:-2], vel, 'b--o')
    vlmark, = axes[1].plot(times[:-2], F1[:,1], 'g--o')
    vqmark, = axes[1].plot(times[:-2], F2[:,1], 'r--o')

    axes[1].legend(  [vmark, vlmark, vqmark],
                [r'$\left(t,v\right)$', r'Linear', r'Quadratic'],
                numpoints=1,
                loc="upper right")
    
    configure(  ax=axes[1],
                title=r'Velocity vs Time $\left(v_t='+str(vT)+r'\right)$',
                xlabel=r"Time $\left(seconds\right)$",
                ylabel=r"Velocity $\left(\frac{meters}{second}\right)$",
                xbounds=None, ybounds=None)

    fig.suptitle("Falling Coffee Filter", size=30)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.08)

    show()

#sliding average for smoothing
def slide_avg(x):
    avg = []
    for i in range(len(x)-2):
        avg.append((x[i]+x[i+1]+x[i+2])/3)
    return np.array(avg)

def falld1(t,x):
    return np.array([x[1],g*(1-(x[1]/vT))])

def falld2(t,x):
    return np.array([x[1],g*(1-(x[1]/vT)**2)])
# this is when I was going to run rk with the two falling models

#~ Entry point of the script.
if __name__ == "__main__":
    #~ State Variables
    t0 = 0; tf = 10*pi; dt = .1
    times = np.arange(t0, tf, dt)

    g = -9.81
    x0 = 1.; v0 = 0.; vT = -0.425; a0 = g
    state = [np.array([x0,v0,a0])]
    #/~ State Variables

    times, dts, pos, vel, acc = sample_data()
    
    plot_samples(times, pos, vel, acc)
    plot_acc_vs_vel(vel,acc)
    plot_samples_and_models(times, pos, vel, acc)
