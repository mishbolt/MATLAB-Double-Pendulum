
close all; 
clear; clc

lam = 1; mu = 1;

X0 = [0.001,0,0,0];
tf = 300;

options = odeset('RelTol',1e-10,'AbsTol',1e-13);

[t,x] = ode45(@derivs,[0 tf],X0,options);

figure(1)
plot(t,x(:,1),'k.-')

wp = sqrt(2+sqrt(2)); wn = sqrt(2-sqrt(2));
Xp = x(:,1) + 1/sqrt(2)*x(:,2);
Xn = x(:,1) - 1/sqrt(2)*x(:,2);
Xpdot = x(:,3) + 1/sqrt(2)*x(:,4);
Xndot = x(:,3) - 1/sqrt(2)*x(:,4);

Hp = 1/2*(Xpdot.^2 + wp^2*Xp.^2);
Hn = 1/2*(Xndot.^2 + wn^2*Xn.^2);

a = 5e-5;

figure(2)
plot(t,Hp,'k.-')
ylim([-a a])
xlabel('time'); ylabel('energy')
title('H_+ as a function of time')

figure(3)
plot(t,Hn,'k.-')
ylim([-a a])
xlabel('time'); ylabel('energy')
title('H_- as a function of time')


% ordering of variables: X = [th1,th2,th1dot,th2dot]

function Ydot = derivs(t,X)

Ydot = zeros(4,1);

lam = 1; mu = 1;

Ydot(1) = X(3);

Ydot(2) = X(4);

Ydot(3) = (-1).*lam.^(-1).*(1+mu+(-1).*cos (X(1)+(-1).*X(2)).^2).^(-1).*(sin(X(1)) ...
  +mu.*sin(X(1))+X(4).^2.*sin(X(1)+(-1).*X(2))+lam.*X(3).^2.*cos(X(1)+(-1).*X(2)) ...
  .*sin(X(1)+(-1).*X(2))+(-1).*cos (X(1)+(-1).*X(2)).*sin(X(2)));

Ydot(4) = (-1).*((-1)+(-1).*mu+cos (X(1)+(-1).*X(2)).^2).^(-1).*(cos (X(1)+(-1).* ...
  X(2)).*sin(X(1))+mu.*cos (X(1)+(-1).*X(2)).*sin(X(1))+lam.*X(3).^2.*sin(X(1)+(-1) ...
  .*X(2))+lam.*mu.*X(3).^2.*sin(X(1)+(-1).*X(2))+X(4).^2.*cos (X(1)+(-1).*X(2)).* ...
  sin(X(1)+(-1).*X(2))+(-1).*sin(X(2))+(-1).*mu.*sin(X(2)));

end

