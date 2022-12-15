from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

#Jd=J-J_theta
def dumbbell(y, t, beta, mu, m0, I0, Jd):
	R, dotR = y
	return np.array([ dotR, beta*mu*np.sin(mu*t)*dotR/(1+beta*np.cos(mu*t))+ \
		((1+m0)*Jd/m0/(1+beta*np.cos(mu*t)))**2/R**3-\
		(1+m0)*(1+beta*np.cos(mu*t))/R**2 ])

t = np.arange(0, 400, 0.01)

R0=10
Omega=1/np.sqrt(R0**3)
beta=0.01
mu=Omega*2
m0=0.1
I0=0.02
k=1
Jd=np.sqrt(k*R0*(1+beta)**3/(1+m0))*m0

#epicycle period
print(2*np.pi*np.sqrt(R0**3/(1+m0)))
#circular orbit  period
print(2*np.pi/Omega)

track = odeint(dumbbell, (R0, 0), t, args=(beta, mu, m0, I0, Jd))

R, dotR=track[:,0], track[:,1]

fig = plt.figure()
ax = fig.add_subplot(121)
ax.plot(t, R, '-')
ax.set_xlabel('t')
ax.set_ylabel(r'$R$')
ax.grid(False)
ax.set_title( r'$\beta$='+str(round(beta,2)) +r', $\mu$='+str(round(mu, 2)) +r', $m_0$='+str(round(m0,2)) \
	+r', $I_0$='+str(round(I0, 2)))

ax = fig.add_subplot(122)
ax.plot(t, dotR, '-')
ax.set_xlabel('t')
ax.set_ylabel(r'$\dot R$')
ax.grid(False)

plt.tight_layout()
plt.show()