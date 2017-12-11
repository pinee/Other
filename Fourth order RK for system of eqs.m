%Modified 4th order Runga Kutte method for system of equations
dx=0.1;
dt=0.1;
IM=101;
NM=121;
y=1.4;
%Initial conditions
for i=1:51
    u(i,1)=0;
    d(i,1)=1;
    p(i,1)=1;
end
for i=51:IM
    u(i,1)=0;
    p(i,1)=0.1;
    d(i,1)=0.125;
end
for i=1:IM
    et(i,1)=(0.5*u(i,1)*u(i,1))+p(i,1)/(d(i,1)*(y-1));
    %3D array with each column as vector Q at ith position and nth time, i
    %is the 2nd dimension and n is the 3rd dimension
    Q(:,i,1)=[d(i,1);d(i,1)*u(i,1);d(i,1)*et(i,1)];
    E(:,i,1)=[d(i,1)*u(i,1);d(i,1)*u(i,1)*u(i,1)+p(i,1);u(i,1)*(d(i,1)*et(i,1)+p(i,1))];
end
for n=1:120
Q2(:,1,n)=Q(:,2,n);
for i=2:IM-1
K1(:,i,n)=(E(:,i+1,n)-E(:,i-1,n))/(2*dx);
Q2(:,i,n)=Q(:,i,n)-(dt/4)*(K1(:,i,n));
end
Q2(:,IM,n)=Q(:,IM-1,n);
for i=1:IM
E2(:,i,n)=myfunc(Q2(:,i,n));
end

Q3(:,1,n)=Q3(:,2,n);
for i=2:IM-1
K2(:,i,n)=(E2(:,i+1,n)-E2(:,i-1,n))/(2*dx);
Q3(:,i,n)=Q(:,i,n)-(dt/3)*(K2(:,i,n));
end
Q3(:,IM,n)=Q3(:,IM-1,n);
for i=1:IM
E3(:,i,n)=myfunc(Q3(:,i,n));
end

Q4(:,1,n)=Q4(:,2,n);
for i=2:IM-1
K3(:,i,n)=(E3(:,i+1,n)-E3(:,i-1,n))/(2*dx);
Q4(:,i,n)=Q(:,i,n)-(dt/2)*(K3(:,i,n));
end
Q4(:,IM,n)=Q4(:,IM-1,n);
for i=1:IM
E4(:,i,n)=myfunc(Q4(:,i,n));
end

Q(:,1,n+1)=Q(:,2,n+1);
for i=2:IM-1
K4(:,i,n)=(E4(:,i+1,n)-E4(:,i-1,n))/(2*dx);
Q(:,i,n+1)=Q(:,i,n)-dt*(K4(:,i,n));
end
Q(:,IM,n+1)=Q(:,IM-1,n+1);
for i=1:IM
d(i,n+1)=(Q(1,i,n+1));
u(i,n+1)=(Q(2,i,n+1))/d(i,n+1);
et(i,n+1)=(Q(3,i,n+1))/d(i,n+1);
p(i,n+1)=(et(i,n+1)-(0.5*u(i,n+1)*u(i,n+1)))*d(i,n+1)*(y-1);
E(:,i,n+1)=[d(i,n+1)*u(i,n+1);d(i,n+1)*u(i,n+1)*u(i,n+1)+p(i,n+1);u(i,n+1)*(d(i,n+1)*et(i,n+1)+p(i,n+1))];
end
%j=i+0.5
for j=1:IM-1
U(j,n)=(u(j,n)+u(j+1,n))/2;
DE(j,n)=(d(j,n)+d(j+1,n))/2;
P(j,n)=(p(j,n)+p(j+1,n))/2;
A(j,n)=sqrt(y*p(j+1,n)/d(j+1,n));
end



%here F is the jacobian matrix as A is used for speed of sound at cell
%interfaces
for j=1:IM-1
F=[0,1,0;0.5*(y-3)*U(j,n)*U(j,n),-(y-3)*U(j,n),(y-1);-(U(j,n)*A(j,n)*A(j,n)/(y-1))+(0.5*y-1)*U(j,n)*U(j,n)*U(j,n),(A(j,n)*A(j,n)/(y-1))+(1.5-y)*U(j,n)*U(j,n),y*U(j,n)];
[X,D]=eig(F);
W=inv(X);
X(:,:,j)=X;
al(:,j)=[U(j,n);U(j,n)+A(j,n);U(j,n)-A(j,n)];
del(:,j)=W*(Q(:,j+1,n)-Q(:,j,n));
end


for j=2:IM-1
for k=1:3
G(k,j)=minmod(del(k,j-1),del(k,j));
end
end
 
 for j=1:IM-1
 for k=1:3
if del(k,j)==0
    beta(k,j)=0;
else
    beta(k,j)=sigma(al(k,j))*((G(i+1)-G(i))/del(k,j));
end
 end
end
 
for j=2:IM-2
for k=1:3
phi(k,j)=sigma(al(k,j))*(G(j+1)+G(j))-si(al(k,j)+beta(k,j))*del(k,j);
end
R(:,j)=X(:,:,j)*phi(:,j);
end

for i=3:IM-2
Q(:,i,n+1)=Q(:,i,n+1)-(dt/dx)*0.5*(R(:,i)-R(:,i-1));
end
Q(:,1,n+1)=Q(:,3,n+1);
Q(:,2,n+1)=Q(:,3,n+1);
Q(:,IM,n+1)=Q(:,IM-2,n+1);
Q(:,IM-1,n+1)=Q(:,IM-2,n+1);

for i=1:IM
d(i,n+1)=(Q(1,i,n+1));
u(i,n+1)=(Q(2,i,n+1))/d(i,n+1);
et(i,n+1)=(Q(3,i,n+1))/d(i,n+1);
p(i,n+1)=(et(i,n+1)-(0.5*u(i,n+1)*u(i,n+1)))*d(i,n+1)*(y-1);
end
end
x=0:dx:10;
for n=1:120
    plot(x,d(:,n));
    hold on;
end
