%Modified Runga Kutta with TVD
dx=1;
dt=1;
IM=11;
NM=121;
x=0:dx:4;
t=0:dt:120;
y=1.4;
%Initial conditions
for i=1:5
    u(i,1)=0;
    d(i,1)=1;
    p(i,1)=1;
end
for i=6:IM
    u(i,1)=0;
    d(i,1)=0.125;
    p(i,1)=0.1;
    
end
for i=1:11
    et(i,1)=(0.5*u(i,1)*u(i,1))+p(i,1)/(d(i,1)*(y-1));
    %3D array with each column as vector Q at ith position and nth time, i
    %is the 2nd dimension and n is the 3rd dimension
    Q(:,i,1)=[d(i,1);d(i,1)*u(i,1);d(i,1)*et(i,1)];
    E(:,i,1)=[d(i,1)*u(i,1);d(i,1)*u(i,1)*u(i,1)+p(i,1);u(i,1)*(d(i,1)*et(i,1)+p(i,1))];
end
for n=1:121
for i=2:IM-1
K1(:,i,n)=(E(:,i+1,n)-E(:,i-1,n))/(2*dx);
end
K1(:,1,n)=K1(:,2,n);
K1(:,IM,n)=K1(:,IM-1,n);
for i=1:IM
Q2(:,i,n)=Q(:,i,n)-(dt/4)*(K1(:,i,n));
end
for i=1:IM
q=Q2(:,i,n)';
E2(:,i,n)=myfunc(q);
end
for i=2:IM-1
K2(:,i,n)=(E2(:,i+1,n)-E2(:,i-1,n))/(2*dx);
end
K2(:,1,n)=K2(:,2,n);
K2(:,IM,n)=K2(:,IM-1,n);
for i=1:IM
Q3(:,i,n)=Q(:,i,n)-(dt/3)*(K2(:,i,n));
end
E3(:,i,n)=myfunc(Q3(:,i,n));
for i=2:IM-1
K3(:,i,n)=(E3(:,i+1,n)-E3(:,i-1,n))/(2*dx);
end
K3(:,1,n)=K3(:,2,n);
K3(:,IM,n)=K3(:,IM-1,n);
for i=1:IM
Q4(:,i,n)=Q(:,i,n)-(dt/2)*(K3(:,i,n));
end
E4(:,i,n)=myfunc(Q4(:,i,n));
for i=2:IM-1
K4(:,i,n)=(E4(:,i+1,n)-E4(:,i-1,n))/(2*dx);
end
K4(:,1,n)=K4(:,2,n);
K4(:,IM,n)=K4(:,IM-1,n);
for i=1:IM
Q(:,i,n+1)=Q(:,i,n)-dt*(K4(:,i,n));
d(i,n+1)=Q[1,i,n+1];
end
end
A=[0,1,0;0.5*(y-3)*u*u,-(y-3)*u,(y-1);-y*p*u/(d*(y-1))+(0.5y-1)*u.^3,y*p/(d*(y-1))+(1.5-y)*u*u,y*u];
[X,D]=eig(A);
W=inv (X);
for n=1:NM
    plot(x,d(:,n));
    hold on;
end
