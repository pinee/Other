
clear all; clc;
%IV is the vector which stores the different initial values of f'
IV(1)=0;
IV(2)=1;
%O is the vector which stores the f' value at the boundary after applying
%2nd order RK
%RK 2nd order method is coded in the rk function
O(1)=rk(0);
O(2)=rk(1);
%temp matrix stores the deviation of boundary value from the desired value1
temp=O-1;
%Here we have assumed two initial values 0 and 1 for f'. Then we
%overshooted the soln my RK method with this assumed initial value. The
%initial value was changed using bisection method for 10 iterations in
%order to get the desired given boundary value of f'
t1=1;
t2=2;
%Number of iterations for initial value
n=10;
%Code for bisection method
for j=3:n
    IV(j)=(IV(t1)+IV(t2))/2;
    O(j)=rk(IV(j));
    temp(j)=O(j)-1;
if (temp(j)*temp(t1)<0)
    t1=t1;
else
    t1=t2;
end
t2=j;
end

fprintf('The appropriate initial value for fdot is')
display(IV(10));
%Final plots
h=0.1;
N=150;
eta=0:h:15;
y10=IV(n);
y20=0;
y30=0;
Y(1,:)=[y10 y20 y30];
%Runga-Kutta scheme
%Here we take a very large finite value instead of infinity
for i=1:N
k1=h*f(Y(i,:));
k2=h*f(Y(i,:)+k1/2);
Y(i+1,:)=Y(i,:)+k2;
end
figure('Name','eta Vs f');
plot(Y(:,3),eta');
xlabel('f')
ylabel('eta')
pause;
figure('Name','eta Vs fdot') ;
plot(Y(:,2),eta');
xlabel('fdot')
ylabel('eta')
pause;
figure ('Name','eta Vs f"');
plot(Y(:,1),eta');
xlabel('f"')
ylabel('eta')
