%Supersonic inflow and subsonic outflow in a wind tunnel
%1st order Steger and Warming
%Initial conditions
JM=101;
NM=9000;
y=1.4;
R=1716;
dx=0.1;
dt=0.00001;
for j=1:29
    u(j,1)=1676.55 ;
end
for j=29:JM
    u(j,1)=572.76 ;
end
for j=1:JM
%density
d(j,1)=0.002241;
%total specific energy
p(j,1)=2000;
%pressure
et(j,1)=p(j,1)/(d(j,1)*(y-1))+0.5*u(j,1)*u(j,1);
%momentum
m(j,1)=d(j,1)*u(j,1);
%total energy
E(j,1)=et(j,1)*d(j,1);
end
for j=1:JM
x(j)=(j-1)*dx;
%cross sectional area of wind tunnel at x(j)
S(j)=1.398+0.347*tanh(0.8*x(j)-4);
syms z;
f=@(z) (1.398+0.347*tanh(0.8*z-4));
%derivative of S(j)
dvt(j)= eval( (subs(diff(f,z,1),z,x(j))) );
end
for j=1:JM
Q(:,j,1)=[(d(j,1));(m(j,1));(E(j,1))];
%flux matrix
F(:,j,1)=S(j)*[d(j,1)*u(j,1);(d(j,1)*u(j,1)*u(j,1)+p(j,1));u(j,1)*(d(j,1)*et(j,1)+p(j,1))];
%positive splitted flux matrix
Fp(:,j,1)=F(:,j,1);
%Negative splitted flux matrix
Fm(:,j,1)=[0 0 0];
H(:,j,1)=[0;dvt(j)*p(j,1);0];
end
for n=1:NM    
Q(:,1,n+1)=(Q(:,1,n));
for j=2:JM-1
Q(:,j,n+1)=Q(:,j,n)-0.5*(dt/dx)*(1/S(j))*(Fp(:,j,n)-Fp(:,j-1,n))+((1/S(j))*dt*(H(:,j,n)));
end
Q(:,JM,n+1)=(Q(:,JM-1,n+1));
for j=1:JM
d(j,n+1)=(Q(1,j,n+1));
u(j,n+1)=(Q(2,j,n+1))/d(j,n+1);
et(j,n+1)=(Q(3,j,n+1))/d(j,n+1);
m(j,n+1)=d(j,n+1)*u(j,n+1);
p(j,n+1)=(et(j,n+1)-(0.5*u(j,n+1)*u(j,n+1)))*d(j,n+1)*(y-1);
F(:,j,n+1)=S(j)*[d(j,n+1)*u(j,n+1);d(j,n+1)*u(j,n+1)*u(j,n+1)+p(j,n+1);u(j,n+1)*(d(j,n+1)*et(j,n+1)+p(j,n+1))];
Fp(:,j,n+1)=F(:,j,n+1);
Fm(:,j,n+1)=[0 0 0];
H(:,j,n+1)=[0;dvt(j)*p(j,n+1);0];
end
end

x= 0:dx:10;
for n=9000
    plot(x,p(:,n))
    axis([0 10 500 5000])
    hold on;
end