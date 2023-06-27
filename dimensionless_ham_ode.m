clc
clear all


[t,x] = ode45(@derivative,[0,100],[pi/5,0,0,0]);

plot(x(:,1),x(:,2))


function dxdt = derivative(t,X)


dxdt = [X(3)-1;2*X(4)-1;-2*X(1);-X(2)];


end