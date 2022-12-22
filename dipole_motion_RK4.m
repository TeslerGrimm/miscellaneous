%motion of an ideal electric dipole in a Coulomb 
% potential, solved with Runge-Kutta 4th Order Method
clc
clear

Q=1;
q=1;
K=1;
m=1;
d=1e-4;

%ODE:
f1=@(y1,y2,y3,y4,y5,y6,t)y2;
f2=@(y1,y2,y3,y4,y5,y6,t)((Q*q*d*K*cos(y3-y5)/y1^3/m)+y4^2*y1);
f3=@(y1,y2,y3,y4,y5,y6,t)y4;
f4=@(y1,y2,y3,y4,y5,y6,t)(((Q*q*d*K*sin(y5-y3)/m/y1^3)-2*y4*y2)/y1);
f5=@(y1,y2,y3,y4,y5,y6,t)y6;
f6=@(y1,y2,y3,y4,y5,y6,t)(Q*q*K*sin(y3-y5)/2/m/y1^2/d);

%initial condition:
r0=0.5;
y(1,1)=r0;
y(2,1)=-0.1;
%
theta0=30*pi/180;
y(3,1)=theta0;
y(4,1)=0;
%
phi0=pi/2;
y(5,1)=phi0;
y(6,1)=0;

%time and step:
t0=0;
tend=2;
h=1e-4;
N=round((tend-t0)/h);
t=zeros(N,1);

%4th order RK:
for n=1:N
    %
    K(1,1)=f1(y(1,n),y(2,n),y(3,n),y(4,n),y(5,n),y(6,n),t(n));
    K(2,1)=f2(y(1,n),y(2,n),y(3,n),y(4,n),y(5,n),y(6,n),t(n));
    K(3,1)=f3(y(1,n),y(2,n),y(3,n),y(4,n),y(5,n),y(6,n),t(n));
    K(4,1)=f4(y(1,n),y(2,n),y(3,n),y(4,n),y(5,n),y(6,n),t(n));
    K(5,1)=f5(y(1,n),y(2,n),y(3,n),y(4,n),y(5,n),y(6,n),t(n));
    K(6,1)=f6(y(1,n),y(2,n),y(3,n),y(4,n),y(5,n),y(6,n),t(n));
    %
    K(1,2)=f1(y(1,n)+1/2*h*K(1,1),y(2,n)+1/2*h*K(2,1),y(3,n)+1/2*h*K(3,1),y(4,n)+1/2*h*K(4,1),y(5,n)+1/2*h*K(5,1),y(6,n)+1/2*h*K(6,1),t(n)+1/2*h);
    K(2,2)=f2(y(1,n)+1/2*h*K(1,1),y(2,n)+1/2*h*K(2,1),y(3,n)+1/2*h*K(3,1),y(4,n)+1/2*h*K(4,1),y(5,n)+1/2*h*K(5,1),y(6,n)+1/2*h*K(6,1),t(n)+1/2*h);
    K(3,2)=f3(y(1,n)+1/2*h*K(1,1),y(2,n)+1/2*h*K(2,1),y(3,n)+1/2*h*K(3,1),y(4,n)+1/2*h*K(4,1),y(5,n)+1/2*h*K(5,1),y(6,n)+1/2*h*K(6,1),t(n)+1/2*h);
    K(4,2)=f4(y(1,n)+1/2*h*K(1,1),y(2,n)+1/2*h*K(2,1),y(3,n)+1/2*h*K(3,1),y(4,n)+1/2*h*K(4,1),y(5,n)+1/2*h*K(5,1),y(6,n)+1/2*h*K(6,1),t(n)+1/2*h);
    K(5,2)=f5(y(1,n)+1/2*h*K(1,1),y(2,n)+1/2*h*K(2,1),y(3,n)+1/2*h*K(3,1),y(4,n)+1/2*h*K(4,1),y(5,n)+1/2*h*K(5,1),y(6,n)+1/2*h*K(6,1),t(n)+1/2*h);
    K(6,2)=f6(y(1,n)+1/2*h*K(1,1),y(2,n)+1/2*h*K(2,1),y(3,n)+1/2*h*K(3,1),y(4,n)+1/2*h*K(4,1),y(5,n)+1/2*h*K(5,1),y(6,n)+1/2*h*K(6,1),t(n)+1/2*h);
    %
    K(1,3)=f1(y(1,n)+1/2*h*K(1,2),y(2,n)+1/2*h*K(2,2),y(3,n)+1/2*h*K(3,2),y(4,n)+1/2*h*K(4,2),y(5,n)+1/2*h*K(5,2),y(6,n)+1/2*h*K(6,2),t(n)+1/2*h);
    K(2,3)=f2(y(1,n)+1/2*h*K(1,2),y(2,n)+1/2*h*K(2,2),y(3,n)+1/2*h*K(3,2),y(4,n)+1/2*h*K(4,2),y(5,n)+1/2*h*K(5,2),y(6,n)+1/2*h*K(6,2),t(n)+1/2*h);
    K(3,3)=f3(y(1,n)+1/2*h*K(1,2),y(2,n)+1/2*h*K(2,2),y(3,n)+1/2*h*K(3,2),y(4,n)+1/2*h*K(4,2),y(5,n)+1/2*h*K(5,2),y(6,n)+1/2*h*K(6,2),t(n)+1/2*h);
    K(4,3)=f4(y(1,n)+1/2*h*K(1,2),y(2,n)+1/2*h*K(2,2),y(3,n)+1/2*h*K(3,2),y(4,n)+1/2*h*K(4,2),y(5,n)+1/2*h*K(5,2),y(6,n)+1/2*h*K(6,2),t(n)+1/2*h);
    K(5,3)=f5(y(1,n)+1/2*h*K(1,2),y(2,n)+1/2*h*K(2,2),y(3,n)+1/2*h*K(3,2),y(4,n)+1/2*h*K(4,2),y(5,n)+1/2*h*K(5,2),y(6,n)+1/2*h*K(6,2),t(n)+1/2*h);
    K(6,3)=f6(y(1,n)+1/2*h*K(1,2),y(2,n)+1/2*h*K(2,2),y(3,n)+1/2*h*K(3,2),y(4,n)+1/2*h*K(4,2),y(5,n)+1/2*h*K(5,2),y(6,n)+1/2*h*K(6,2),t(n)+1/2*h);
    %
    K(1,4)=f1(y(1,n)+h*K(1,3),y(2,n)+h*K(2,3),y(3,n)+h*K(3,3),y(4,n)+h*K(4,3),y(5,n)+h*K(5,3),y(6,n)+h*K(6,3),t(n)+h);
    K(2,4)=f2(y(1,n)+h*K(1,3),y(2,n)+h*K(2,3),y(3,n)+h*K(3,3),y(4,n)+h*K(4,3),y(5,n)+h*K(5,3),y(6,n)+h*K(6,3),t(n)+h);
    K(3,4)=f3(y(1,n)+h*K(1,3),y(2,n)+h*K(2,3),y(3,n)+h*K(3,3),y(4,n)+h*K(4,3),y(5,n)+h*K(5,3),y(6,n)+h*K(6,3),t(n)+h);
    K(4,4)=f4(y(1,n)+h*K(1,3),y(2,n)+h*K(2,3),y(3,n)+h*K(3,3),y(4,n)+h*K(4,3),y(5,n)+h*K(5,3),y(6,n)+h*K(6,3),t(n)+h);
    K(5,4)=f5(y(1,n)+h*K(1,3),y(2,n)+h*K(2,3),y(3,n)+h*K(3,3),y(4,n)+h*K(4,3),y(5,n)+h*K(5,3),y(6,n)+h*K(6,3),t(n)+h);
    K(6,4)=f6(y(1,n)+h*K(1,3),y(2,n)+h*K(2,3),y(3,n)+h*K(3,3),y(4,n)+h*K(4,3),y(5,n)+h*K(5,3),y(6,n)+h*K(6,3),t(n)+h);
    %
    for i=1:6
        y(i,n+1)=y(i,n)+h/6*(K(i,1)+2*K(i,2)+2*K(i,3)+K(i,4)); 
    end
    t(n+1)=t(n)+h;
end

%observables:
Xc=y(1,:).*cos(y(3,:));
Yc=y(1,:).*sin(y(3,:));
X1=Xc+d*cos(y(5,:))/2;
X2=Xc-d*cos(y(5,:))/2;
Dy=d*sin(y(5,:))/2;
Y1=Yc+Dy;
Y2=Yc-Dy;
r1=sqrt(X1.^2+Y1.^2);
r2=sqrt(X2.^2+Y2.^2);

%plot:
subplot(2,2,1)
plot(Xc,Yc,'k-','LineWidth',1);
hold on
plot(Xc(1),Yc(1),'ko','LineWidth',1,'MarkerFaceColor','g');
hold on
plot(0,0,'ko','LineWidth',1,'MarkerFaceColor','k');
hold on
plot([Xc(1)-0.1*cos(phi0),Xc(1)+0.1*cos(phi0)],[Yc(1)-0.1*sin(phi0),Yc(1)+0.1*sin(phi0)],'ko-','LineWidth',1);
hold on
plot(Xc(1)-0.1*cos(phi0),Yc(1)-0.1*sin(phi0),'ko','markerfacecolor','b');
hold on
plot(Xc(1)+0.1*cos(phi0),Yc(1)+0.1*sin(phi0),'ko','markerfacecolor','r');
xlabel('x');
ylabel('y');
grid on

subplot(2,2,2)
plot(t,y(6,:),'k-','LineWidth',1)
xlabel('t');
ylabel('\Omega_{\phi}');
grid on

subplot(2,2,3)
plot(t,1./r1-1./r2,'k-','LineWidth',1)
xlabel('t');
ylabel('potnetial energy');
grid on

subplot(2,2,4)
plot(t,sin(y(3,:)-y(5,:)),'k-','LineWidth',1)
xlabel('t');
ylabel('sin(\theta-\phi)');
grid on




