from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

beta=0.9
mu=5
phi=0

def equ(state, t):

    dydx = np.zeros_like(state)
    dydx[0] = state[1]
    dydx[1] = beta*mu*sin(mu*t+phi)*state[1]/(1+beta*cos(mu*t+phi))+\
    -sin(state[0])*cos(state[0])*state[1]**2/(1/3+sin(state[0])**2)\
    +sin(state[0])/(1/6+0.5*sin(state[0])**2)  

    return dydx

# create a time array
dt = 0.001
t = np.arange(0.0, 8, dt)

# th is the initial angles
# w is initial angular velocities
th = 0.00001
w = 0

# initial state
#state = np.radians([th, w])
state = [th, w]

# integrate
y = integrate.odeint(equ, state, t)

x1 = 0.5*sin(y[:, 0])
y1 = 0*y[:, 0]

x2 = -0.5*sin(y[:, 0])
y2 = cos(y[:, 0])

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1, 1), ylim=(0, 1))
ax.grid()
ax.set_title( r'$\beta$='+str(beta) +r', $\mu$='+str(round(mu)) +r', $\phi$='+str(round(phi/np.pi*180)) +' degrees')

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1f'
time_text = ax.text(0.05, 0.5, '', transform=ax.transAxes)


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

n=10

def animate(i):
    thisx = [x1[n*i], x2[n*i]]
    thisy = [y1[n*i], y2[n*i]]

    line.set_data(thisx, thisy) # see matplotlib.lines
    line.set_color('pink')
    line.set_linewidth(10)
    time_text.set_text(time_template % (n*i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, round(len(y)/n)-1),
                              interval=1, blit=True, init_func=init)

#ani.save('fall_2.gif')
plt.show()