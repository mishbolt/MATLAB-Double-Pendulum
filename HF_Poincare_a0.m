
clear all; close all; clc

qold1 = pi/5; pold1 = 0;
qold2 = -qold1; pold2 = 0;

% Hamiltonian value for P sect 
E0 = -2*cos(qold1) - cos(qold1 + qold2);

% solve for set of ICs (alpha,beta(alpha),0,0)
% with the same energy
N=15;
qold1 = linspace(-qold1,qold1,N);
qold2 = acos(-(E0 + 2*cos(qold1))) - qold1;

% Set the options so ode45 detects the events
options = odeset('Events',@pendevents,...
    'RelTol',1e-9,'AbsTol',1e-12);

% Set the loop for various initial condtions
for i=1:N
% Set initial parameters
m1 = 1;
m2 = 1;
L1 = 1;
L2 = 1;
g = 9.8;
tf = 2000;
x0 = [qold1(i);qold2(i);0;0];

% p2 = (0.02557+0.02948)/2;
% E0 = -2.99954308546917;
% p1 = 1/2*(2*p2+sqrt(4*p2^2+8*(3+E0-p2^2)));
% 
% N =0 ;
% x0 = [0;0;p1;p2];
[T,x,te,ye,ie] = ode45(@derivative,[0,tf],x0,options);

plot(ye(:,2),ye(:,4),'r.')
xlabel('\beta'); ylabel('p_\beta')
title(['Poincare Section for E_0 = ',num2str(E0)]);
hold on

end

% plot boundary curve; alphadot = 0 is boundary
% since Psect: alphadot > 0, f1 = 0 = H - E0 defines 
% implicit function for beta, p_beta (1d curve)
maxB = acos(-2-E0);

f1 = @(x,y) -E0 - 2 - cos(x) + y.^2/2;
fimplicit(f1,'k-')


% Function for ode45
function dxdt = derivative(t,X)
dxdt1 = 2*(X(3)-(1+cos(X(2)))*X(4))/(3-cos(2*X(2)));
dxdt2 = 2*(-(1+cos(X(2)))*X(3)+(3+2*cos(X(2)))*X(4))/(3-cos(2*X(2)));
dxdt3 = -2*sin(X(1))-sin(X(1)+X(2));

A = -sin(X(1)+X(2))-(2*sin(X(2))*(X(3)-X(4))*X(4))/(3-cos(2*X(2)));

B = (2*sin(2*X(2)))*((X(3))^2-2*(1+cos(X(2)))*X(3)*X(4)+(X(4))^2*(3+2*cos(X(2))))/((3-cos(2*X(2)))^2);

dxdt4 = A+B;

dxdt = [dxdt1;dxdt2;dxdt3;dxdt4];

end


% Function for events
function [value,isterminal,direction] = pendevents(t,y)

value = y(1);
isterminal = 0;
direction = 1;

end