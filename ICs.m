clc
clear all

p2 = (0.127+0.023)/2;

E0 = [-2.12132034355964];

p1 = 1/2*(2*p2+sqrt(4*p2^2+8*(3+E0-p2^2)));
x0 = [0,0,p1,p2];
options = ['RelTol',1e-8,'AbsTol',1e-10];


[t,Y] = ode45(@derivative,[0:0.1:100],x0,options);

function dxdt = derivative(t,X)
mu = 1;
lam = 1;

dxdt1 = (X(3)*1/lam-X(4)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);

dxdt2 = (lam*(mu+1)*X(4)-X(3)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);

A = -(X(3)*X(4)*sin(X(1)-X(2)))/((mu+(sin(X(1)-X(2)))^2));
B = (sin(2*(X(1)-X(2)))*(1/2*1/lam*(X(3))^2+1/2*lam*(mu+1)*(X(4))^2-X(3)*X(4)*cos(X(1)-X(2))))/(mu+(sin(X(1)-X(2)))^2)^2;

dxdt3 = A+B-(mu+1)*sin(X(1));

dxdt4 = -A-B-1/lam*sin(X(2));

dxdt = [dxdt1;dxdt2;dxdt3;dxdt4];



end