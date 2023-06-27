close all
clc
clear
mu = 1;
lam = 1;

q1(1) = 0; p1(1) = 0.5;
q2(1) = 0; p2(1) = 0;

E0 = H([q1(1),q2(1),p1(1),p2(1)]);

N = 20;

p1 = linspace(p1(1),0,N);
p2 = (p1+sqrt(p1.^2-2*lam*(mu+1)*...
    (1/2*lam^(-1)*p1.^2-mu^2-mu-mu*lam^(-1)-E0*mu)))/(lam*(mu+1));

den = mu+(sin(q1-q2)).^2;

omega1 = 1/lam.*p1./den-p2.*cos(q1-q2)./den;
omega2 = lam.*(mu+1).*p2./den-p1.*cos(q1-q2)./den;

options = odeset('Events',@pendevents,'RelTol',1e-6,'AbsTol',1e-6);
% Set the loop for various initial condtions

for i=1:N
% Set initial parameters

tf = 2000;

% omega2 = 0.25;
% omega1 = (-omega2+sqrt(omega2.^2-2*lam*(mu+1)*...
%     (1/2*lam^(-1)*omega2.^2-mu-1-lam^(-1)-E0)))/(lam*(mu+1));

%x0 = [0;0;omega1;omega2];

x0 = [q1(1);q2(1);omega1(i);omega2(i)];

[T,x,te,ye,ie] = ode45(@derivative,[0,tf],x0,options);
 

plot(ye(:,2),ye(:,4),'r.')
xlabel('\theta_2'); ylabel('p_\theta_2')

title(['Poincare Section, $\theta_1 = 0$, $\dot{\theta_1}>0$, $E_0 \approx$', num2str(round(E0,2))],'Interpreter','latex')
hold on

end

%Boundary Curve

%boundary = @(x,y) E0-(-1-mu+1/2.*lam.*y.^2-cos(x)./lam);
%fimplicit(boundary,'k-')



% Function for ode45

%Algebraically computed
function dxdt = derivative(t,X)
mu = 1;
lam = 1;

dxdt1 = X(3);
dxdt2 = X(4);
A = (1/lam*sin(X(2))*cos(X(1)-X(2))-(X(3))^2*sin(2*(X(1)-X(2)))-(X(4))^2*...
    sin(X(1)-X(2))*lam-lam*(mu+1)*sin(X(1)))/(lam^2*(mu+1)-(cos(X(1)-X(2)))^2);

dxdt3 = A;
dxdt4 = (-1/lam*sin(X(2))-A*cos(X(1)-X(2))+(X(3))^2*sin(X(1)-X(2)))/lam;


dxdt = [dxdt1;dxdt2;dxdt3;dxdt4];

end


function ham = H(X)
mu = 1;
lam = 1;

ham = (-1).*(1+mu).*cos(X(1))+lam.^(-1).*((X(3).^2+lam.^2.*(1+mu).*X(4).^2+( ...
  -2).*lam.*X(3).*X(4).*cos(X(1)+(-1).*X(2))).*(1+2.*mu+(-1).*cos(2.*(X(1)+( ...
  -1).*X(2)))).^(-1)+(-1).*cos(X(2)));



end


% Function for events
function [value,isterminal,direction]=pendevents(t,y)

value = y(1);
isterminal = 0;
direction = 1;


end