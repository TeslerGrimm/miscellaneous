%Solve 2D Poisson equ with inner and outer boundaries using FDM.
%Boundaries are created via B&W image inputs:
%(region to solve should be painted white)
%outer boundary: Shape1.jpg
%inner boundary: Shape2.jpg
%Dirichlet boundary condition (value) is used.
%Warning: the boundry line should be closed!
%(some details are not explained)

clc
clear
%rectangular canvas: [Lx~Rx,Uy~Dy]
Lx=-1;
Rx=1;
Dy=-1;
Uy=1;
%a point source can be added (q the charge)
xs=0;
ys=0.7;
q=0;
%read the image
i=imread('Shape1.jpg');
%read only the B in RGB
i=i(:,:,1);
[ny,nx]=size(i);
%pixel coordinates of the point source:
nxs=round((xs-Lx)/(Rx-Lx)*nx);
nys=round((ys-Uy)/(Dy-Uy)*ny);
%similarly
j=imread('Shape2.jpg');
j=j(:,:,1);
%kind warning:
if max(size(i)-size(j))~=0 
    error('Sizes of Two Photos Must be Same!');
end
%array marking the boundary (boundary point: 1, otherwise: 0)
P1=edge(i,'canny');
%(canny is just one choice to identify the edges)
%sequence number of boundary points:
B1=find(P1==1);
%similarly
P2=edge(j,'canny');
B2=find(P2==1);

P=P1+P2;
[y,x]=find(P==1);

B=find(P==1);

%image binarization (black: 0, white: 1)
Q1=imbinarize(i,graythresh(i));
Q2=imbinarize(j,graythresh(j));
Q=Q1+Q2-2*P;
%in Q, common black boundary point takes value 0,
%uncommon black boundary point takes value -1,
%uncommon black non-boundary point takes value 1
OQ1=Q+3*Q1;
OQ2=Q+3*Q2;

%now the boundary points:
p=find(Q==2);
[Y,X]=find(Q==2);

%uncommon black non-boundary points:
p1=find(OQ1==1);
p2=find(OQ2==1);

%grid setup:
V1=zeros(size(Q));
Q=V1;%now Q is the grid
%boundary condition:
V1(B1)=1;
V1(B2)=3;
V1(p1)=1;
V1(p2)=3;
V2=V1;
%point source:
Q(nys,nxs)=q;
%iteration
E=1;%initial value of relative error
w=1.816;%relaxation parameter
er=1e-5;%precession goal
time=0;
while E>er
    time=time+1;
    for k=1:length(p)
        j=X(k);
        i=Y(k);
        V2(i,j)=V1(i,j)+w*(Q(i,j)+V2(i-1,j)+V2(i,j-1)-4*V1(i,j)+V1(i+1,j)+V1(i,j+1))/4;
    end
    E=max(max(abs(V1-V2)));
    V1=V2;
end
%done!

disp('times of iteration:');
disp(time);
disp('relative error:');
disp(E);

Rx=linspace(Lx,Rx,nx);
Ry=linspace(Dy,Uy,ny);
ry=fliplr(Ry);

%potential plot
figure(1)
colormap(jet(500));
surf(Rx,ry,V2,V2,'linestyle','none','FaceAlpha',1,'LineWidth',1);
shading interp
view([0,90]);
colorbar;
grid on
xlabel('x');ylabel('y');zlabel('V');
title({
     'number of grid nodes:';
     [num2str(nx),'×',num2str(ny),'=',num2str(nx*ny)];
     });
hold on

%check the value of potential in a measuring point (ax,ay)
ax=0;
ay=0;
h=Rx(2)-Rx(1);
Ox=round((ax-Lx)/h);
Oy=round((Uy-ay)/h);
disp('potential:');
U=V2(Oy,Ox);
disp(U);

%field lines plot
figure(2)
%subplot(1,2,1)
colormap(jet(500));
contour(Rx,ry,V2,30,'LineWidth',1.5);
grid on
hold on
[Ex,Ey]=gradient(-V2);
[RX,RY]=meshgrid(Rx,ry);
sB=B(1:5:end);%field lines start at boundary
ef=streamline(RX,RY,Ex,-Ey,RX(sB),RY(sB));
set(ef,'LineWidth',1,'Color','k');
axis([-1,1,-1,1]);
axis equal
xlabel('x');
ylabel('y');
hold on
%draw the measuring point:
plot(Rx(Ox),ry(Oy),'ko','Linewidth',1,'markerfacecolor','g');
hold on
%plot(RX(p),RY(p),'r.','Linewidth',1);%inner point
plot(RX(B),RY(B),'k.','Linewidth',1);%boundary point

%total electrostatic energy:
energy=sum(sum(E.^2))*h^2;
disp('energy:');
disp(energy);

%next steps: 
%computation of charge density on the bondary;
%other boundary conditions
