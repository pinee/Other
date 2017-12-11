%1D compressible euler equation applied over staggered rectangular grid
%Supersonic inflow and outflow in a wind tunnel
%Notations used in the code 
%u:velocity; d:density;p:pressure; et:total specific energy; m:momentum; 
%E:total energy; Q:variable matrix; F:Flux matrix


%Discretization of domain
JM=51;
NM=1500;
dx=0.2;
dt=0.0001;

%Known constants
y=1.4;
R=1716;

%Setting up of Initial conditions for mid-points of all control volumes
%Set all the properties in the domain equal to the given initial values
for j=1:JM
u(j,1)=1676.55;
d(j,1)=0.002241;
et(j,1)=3636204;
p(j,1)=(et(j,1)-0.5*u(j,1)*u(j,1))*d(j,1)*(y-1);
m(j,1)=d(j,1)*u(j,1);
E(j,1)=et(j,1)*d(j,1);
end

%Value of x, cross sectional area of the nozzle and the corresponding
%derivative at midpoints of every control volume
for j=1:JM
x(j)=(j-1)*dx;
%cross sectional area of wind tunnel at x(j)
S(j)=1.398+0.347*tanh(0.8*x(j)-4);
syms z;
f=@(z) (1.398+0.347*tanh(0.8*z-4));
%derivative of S(j)
dvt(j)= eval( (subs(diff(f,z,1),z,x(j))) );
end
%Flux splitting is used to find the flux at cell interfaces
%Setting up of initial values for variable matrix and flux matrices at
%every control volume
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
%Time advancing by first order forward difference approximation 
for n=1:NM 
 %Inlet boundary condition   
Q(:,1,n+1)=(Q(:,1,n));
%1st order Steger and Warming flux vector explicit scheme
for j=2:JM-1
    %A backward difference approximation is used for the positive terms and
    %a forward difference approximation is used for the negative terms
    Q(:,j,n+1)=Q(:,j,n)-0.5*(dt/dx)*(1/S(j))*(Fp(:,j,n)-Fp(:,j-1,n))+((1/S(j))*dt*(H(:,j,n)));
end
%Properties at the exit plane using simple extrapolation
Q(:,JM,n+1)=(Q(:,JM-1,n+1));
%Extraction of values of primary variables from the variable matrix
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
%Plotting of graph 
x= 0:dx:10;
for n=1:100:1500
    plot(x,p(:,n))
    axis([0 10 500 2500])
    hold on;
end
figure(2);
n=1:100:1500;
plot(n,p(JM,n)-p(1,n));

