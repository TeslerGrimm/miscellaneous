import numpy as np
import matplotlib.pyplot as plt

K=0
N=6*K+4
r =2*np.sqrt(3)*2/3
I=1/N
p=np.zeros([N**2,2])
A=np.zeros([N**2,2])
B=np.zeros([N**2,2])
#rotation angles for B:
vb=np.arange(0,120,5)
#rotation angles for A:
va=np.array([0])
#output:
o1=np.zeros([len(vb),len(va)])
o2=np.zeros([len(vb),len(va)])
o3=np.zeros([len(vb),len(va)])
o4=np.zeros([len(vb),len(va)])

#build an ET array:
m=1
p[0,0]=-(1-I)/2
p[0,1]=-np.sqrt(3)*(1-I)/6
for i in range(0,N-1):
    for j in range(0,2*N-2*i-2):
        m=m+1
        if np.mod(j,2)==0:
            p[m-1,1]=p[m-2,1]+np.sqrt(3)*I/6
        else:
            p[m-1,1]=p[m-2,1]-np.sqrt(3)*I/6
        p[m-1,0]= p[m-2,0]+I/2
    m=m+1
    p[m-1,1]=-np.sqrt(3)*(1-I)/6+np.sqrt(3)*I*(i+1)/2
    p[m-1,0]= (i+1)*I/2-0.5*(1-I)
#-------------------------------------------
# plt.plot(p[:,0], p[:,1], 'o-')
# plt.show()
#-------------------------------------------
timerA=1
for tA in va:
    a=tA*np.pi/180
    timerB=1
    for tB in vb:
        b=tB*np.pi/180
        B[:,0]=p[:,0]*np.cos(b)-p[:,1]*np.sin(b)
        B[:,1]=p[:,1]*np.cos(b)+p[:,0]*np.sin(b)+r
        A[:,0]=p[:,0]*np.cos(a)-p[:,1]*np.sin(a)
        A[:,1]=p[:,1]*np.cos(a)+p[:,0]*np.sin(a)
        #force exerted on A:
        Fx=0
        Fy=0
        #torque exerted on A:
        Ta=0
        #torque exerted on B:
        Tb=0
        for i in range(0,N**2):
            for j in range(0,N**2):
                k=1
                R=((A[i,0]-B[j,0])**2+(A[i,1]-B[j,1])**2)**((k-1)/2)
                Fx=Fx+(-A[i,0]+B[j,0])*R
                Fy=Fy+(-A[i,1]+B[j,1])*R
                Ta=Ta+(-A[i,1]*(-A[i,0]+B[j,0]) + A[i,0]*(-A[i,1]+ B[j,1]))*R
                Tb=Tb+(-(B[j,1]-r)*(-B[j,0]+A[i,0]) + B[j,0]*(-B[j,1]+A[i,1]))*R
        Fx=Fx/N**4
        Fy=Fy/N**4
        Ta=Ta/N**4
        Tb=Tb/N**4      
        o1[timerB-1,timerA-1]=Fx
        o2[timerB-1,timerA-1]=Fy
        o3[timerB-1,timerA-1]=Ta
        o4[timerB-1,timerA-1]=Tb       
        timerB=timerB+1
    timerA=timerA+1

fig = plt.figure()
ax = fig.add_subplot(121)
ax.plot(vb, o1-o1[1], '-',label='Fx-Fx(0)')
ax.plot(vb, o2-o2[1], '-',label='Fy-Fy(0)')
ax.set_xlabel('angle')
ax.set_ylabel('Force on A')
ax.legend(loc='upper right')
ax.grid(True)
plt.gca().xaxis.set_major_locator(plt.MultipleLocator(30))
plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(10))

ax = fig.add_subplot(122)
ax.plot(vb, o3, '-',label='Torque on A')
ax.plot(vb, o4, '-',label='Torque on B')
ax.set_xlabel('angle')
ax.set_ylabel('Toeque')
ax.legend(loc='upper right')
ax.grid(True)
ax.set_title( '$k$='+str(k) +', $r/d$='+str(round(r/(np.sqrt(3)*2/3),2)) +', $N^2$='+str(N**2) )
plt.gca().xaxis.set_major_locator(plt.MultipleLocator(30))
plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(10))

plt.tight_layout()
#plt.subplots_adjust (left=None,bottom=None,right=None,top=None,wspace=0.5,hspace=0)
plt.show()