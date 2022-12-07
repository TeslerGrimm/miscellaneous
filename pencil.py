from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def pencil(y, t, beta, mu, phi):
	theta, dottheta = y
	return np.array([ dottheta, 6*np.sin(theta) - beta*mu*np.sin(mu*t+phi)*dottheta/(1+beta*np.cos(mu*t+phi)) ])

t = np.arange(0, 10, 0.001)

beta=0.5
mu=100
phi=0

track = odeint(pencil, (0.00001, 0), t, args=(beta,mu,phi))

theta, dottheta=track[:,0], track[:,1]

N=1-3*np.sin(theta)**2-0.5*np.cos(theta)*dottheta**2

f=3*np.sin(theta)*np.cos(theta)-0.5*np.sin(theta)*dottheta**2

r=(9*np.sin(theta)*np.cos(theta)-6*np.sin(theta))/(1-3*np.sin(theta)**2+6*np.cos(theta)**2-6*np.cos(theta))

fig = plt.figure()
ax = fig.add_subplot(121)
ax.plot(theta*180/np.pi, f/N, '-',label='m(t)')
ax.plot(theta*180/np.pi, r, '-',label='m=const')
ax.set_xlabel(r'$\theta$')
ax.set_ylabel('f/N')
plt.xlim(0,30)
plt.ylim(0,50)
ax.legend(loc='upper left')
ax.grid(True)
ax.set_title( r'$\beta$='+str(beta) +r', $\mu$='+str(round(mu)) +r', $\phi$='+str(round(phi/np.pi*180)) +' degrees')

ax = fig.add_subplot(122)
ax.plot(theta*180/np.pi, dottheta**2/24+np.cos(theta)/2, '-',label='m(t)')
ax.plot(theta*180/np.pi, 0.5*np.ones(np.size(t)), '-',label='m=const')
ax.set_xlabel(r'$\theta$')
ax.set_ylabel('H/m')
plt.xlim(0,30)
plt.ylim(0.47,0.53)
ax.legend(loc='upper left')
ax.grid(True)

plt.tight_layout()
plt.show()