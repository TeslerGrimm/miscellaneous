from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def pencil(y, t, beta, mu, phi):
	theta, dottheta = y
	return np.array([ dottheta, 3/2*np.sin(theta) + beta*mu*np.sin(mu*t+phi)*dottheta/(1+beta*np.cos(mu*t+phi)) ])

t = np.arange(0, 15, 0.001)

beta=0.1
mu=100
phi=0

track = odeint(pencil, (0.00001, 0), t, args=(beta,mu,phi))

theta, dottheta=track[:,0], track[:,1]

N=1-3/4*np.sin(theta)**2-0.5*np.cos(theta)*dottheta**2

f=3/4*np.sin(theta)*np.cos(theta)-0.5*np.sin(theta)*dottheta**2

r=(9/4*np.sin(theta)*np.cos(theta)-3/2*np.sin(theta))/(1-3/4*np.sin(theta)**2+3/2*np.cos(theta)**2-3/2*np.cos(theta))

fig = plt.figure()
ax = fig.add_subplot(131)
ax.plot(theta*180/np.pi, f/N, '-',label='m(t)')
ax.plot(theta*180/np.pi, r, '-',label='m=const')
ax.set_xlabel(r'$\theta$')
ax.set_ylabel('f/N')
plt.xlim(0,72)
plt.ylim(-10,0.5)
plt.axhline(0,c="k")
ax.legend(loc='lower left')
ax.grid(True)
ax.set_title( r'$\beta$='+str(beta) +r', $\mu$='+str(round(mu)) +r', $\phi$='+str(round(phi/np.pi*180)) +' degrees')

ax = fig.add_subplot(132)
ax.plot(theta*180/np.pi, N, '-',label='m(t)')
ax.plot(theta*180/np.pi, 1-3/4*np.sin(theta)**2+3/2*np.cos(theta)**2-3/2*np.cos(theta), '-',label='m=const')
ax.set_xlabel(r'$\theta$')
ax.set_ylabel('N')
plt.xlim(0,72)
plt.ylim(-0.5,1.5)
ax.legend(loc='upper left')
ax.grid(True)

# ax = fig.add_subplot(133)
# ax.plot(theta*180/np.pi, dottheta**2/6+np.cos(theta)/2, '-',label='m(t)')
# ax.plot(theta*180/np.pi, 0.5*np.ones(np.size(t)), '-',label='m=const')
# ax.set_xlabel(r'$\theta$')
# ax.set_ylabel('H/m')
# plt.xlim(0,72)
# plt.ylim(0.47,0.53)
# ax.legend(loc='upper left')
# ax.grid(True)

ax = fig.add_subplot(133)
ax.plot(theta*180/np.pi, dottheta, '-',label='m(t)')
ax.plot(theta*180/np.pi, np.sqrt(3*(1-np.cos(theta))), '-',label='m=const')
ax.set_xlabel(r'$\theta$')
ax.set_ylabel(r'$\dot\theta$')
plt.xlim(0,72)
ax.legend(loc='upper left')
ax.grid(True)

print(theta[-1]*180/np.pi)

plt.tight_layout()
plt.show()