from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def pencil(y, t, beta, mu, phi):
	theta, dottheta = y
	return np.array([ dottheta, 10*np.cos(2*t)-5*dottheta -theta + beta*mu*np.sin(mu*t+phi)*dottheta/(1+beta*np.cos(mu*t+phi)) ])

t = np.arange(0, 60, 0.001)

beta=0.4
mu=2
phi=0
theta0=5*np.pi/180

track = odeint(pencil, (theta0, 0), t, args=(beta,mu,phi))

theta, dottheta=track[:,0], track[:,1]

fig = plt.figure()
ax = fig.add_subplot(131)
ax.plot(t, theta*180/np.pi, '-')
#ax.plot(t, theta0*180/np.pi*np.exp(beta*t/2), '-')
ax.set_xlabel('t')
ax.set_ylabel(r'$\theta$')
ax.grid(True)
ax.set_title( r'$\beta$='+str(beta) +r', $\mu$='+str(round(mu)) +r', $\phi$='+str(round(phi/np.pi*180)) +' degrees')

ax = fig.add_subplot(132)
ax.plot(t, dottheta, '-')
ax.set_xlabel('t')
ax.set_ylabel(r'$\dot\theta$')
ax.grid(True)

ax = fig.add_subplot(133)
ax.plot(theta*180/np.pi, dottheta, '-')
ax.set_xlabel(r'$\theta$')
ax.set_ylabel(r'$\dot\theta$')
ax.grid(True)

plt.tight_layout()
plt.show()