clc
clear
B=0;
for A=0:2:120 
    K=2;
    N=6*K+4;  
    h =1*sqrt(3)*2/3;
    I=1/N;
    q=zeros(N^2,2);
    p=q;
    m=1;
    p(1,1)=-(1-I)/2;
    p(1,2)=-sqrt(3)*(1-I)/6;
    for i=0:N-2
        for j = 0:2*N-1-2*i-2
            m=m+1;
            if mod(j,2)==0
                p(m,2)=p(m-1,2)+sqrt(3)*I/6;
            else
                p(m,2)=p(m-1,2)-sqrt(3)*I/6;
            end
            p(m,1)= p(m-1,1)+I/2;
        end
        m=m+1;
        p(m,2)=-sqrt(3)*(1-I)/6+sqrt(3)*I*(i+1)/2;
        p(m,1)= (i+1)*I/2-0.5*(1-I);
    end
    P=p;
    a=A*pi/180;%rotated
    b=B*pi/180;%fixed
    for m=1:N^2
        q(m,1)=P(m,1)*cos(a)-P(m,2)*sin(a);
        q(m,2)=P(m,2)*cos(a)+P(m,1)*sin(a)+h;
    end
    for m=1:N^2
        p(m,1)=P(m,1)*cos(b)-P(m,2)*sin(b);
        p(m,2)=P(m,2)*cos(b)+P(m,1)*sin(b);
    end
    [x,y]=meshgrid(-1:0.06:1,-0.5:0.06:h+0.5);
    [ny,nx]=size(x);
    V=zeros(ny,nx);
    for i=1:N^2
        V=V-1./sqrt((q(i,1)-x).^2+(q(i,2)-y).^2)-1./sqrt((p(i,1)-x).^2+(p(i,2)-y).^2);
    end
    %mesh(x,y,V,'linestyle','-','FaceAlpha',1,'LineWidth',1);
    contour(x,y,V,50,'LineWidth',1);axis equal
    hold on
    plot(q(:,1),q(:,2),'k.-',p(:,1),p(:,2),'b.-','LineWidth',1);
    axis([-1,1,-0.5,h+0.7 ]);
    hold off
    getframe;
end