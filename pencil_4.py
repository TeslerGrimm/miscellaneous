from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def y(y, t, beta, mu, phi):
	theta, dottheta = y
	return np.array([ dottheta, -1 + beta*mu*np.sin(mu*t+phi)*dottheta/(1+beta*np.cos(mu*t+phi)) ])

def x(y, t, beta, mu, phi):
	theta, dottheta = y
	return np.array([ dottheta, beta*mu*np.sin(mu*t+phi)*dottheta/(1+beta*np.cos(mu*t+phi)) ])

t = np.arange(0, 5, 0.001)

beta=0.3
mu=1
phi=np.pi
vx0=np.sqrt(2)/2
vy0=np.sqrt(1-vx0**2)

trackx = odeint(x, (0, vx0), t, args=(beta,mu,phi))
tracky = odeint(y, (0, vy0), t, args=(beta,mu,phi))

x, dotx=trackx[:,0], trackx[:,1]
y, doty=tracky[:,0], tracky[:,1]

fig = plt.figure()
ax = fig.add_subplot(131)
ax.plot(t, x, '-',label='m(t)')
ax.plot(t, vx0*t, '-',label='m=const')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.legend(loc='lower left')
ax.grid(True)
ax.set_title( r'$\beta$='+str(beta) +r', $\mu$='+str(round(mu)) +r', $\phi$='+str(round(phi/np.pi*180)) +' degrees')

ax = fig.add_subplot(132)
ax.plot(t, y, '-',label='m(t)')
ax.plot(t, vy0*t-0.5*t**2, '-',label='m=const')
ax.set_xlabel('t')
ax.set_ylabel('y')
ax.legend(loc='lower left')
ax.grid(True)

ax = fig.add_subplot(133)
ax.plot(x, y, '-',label='m(t)')
ax.plot(vx0*t, vy0*t-0.5*t**2, '-',label='m=const')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.ylim(0,0.5)
plt.xlim(0,2)
ax.legend(loc='upper left')
ax.grid(True)

plt.tight_layout()
plt.show()