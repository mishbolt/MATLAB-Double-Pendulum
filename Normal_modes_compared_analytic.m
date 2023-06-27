
m1 = 1;
m2 = 200;
L1 = 5;
L2 = L1;
L=L1;
g = 9.8;

tf = 30;

time = 0:0.01:tf;

x0 = [pi/10;0;0;0];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[T,y] = ode45(@derivative,[0,tf],x0,options);


%Eigenfrquencies

omega1 = sqrt(2+sqrt(2));
omega2 = sqrt(2-sqrt(2));

%------Normalization Constants and Mode Amplitudes------%

%Generalized eigenvectors
Phi_1 = [1;-sqrt(2)];
Phi_2 = [1;sqrt(2)];

      
%Normalization Constants     
N1 = 1/(sqrt(2-sqrt(2)));
N2 = 1/(sqrt(2+sqrt(2)));

%Mode Amplitudes
A1 = 1/2*N1*x0(1)*(2-sqrt(2));
A2 = 1/2*N2*x0(1)*(2+sqrt(2));


%Solutions
theta_1 = 1/2*x0(1)*cos(omega1*time)*Phi_1(1)+1/2*x0(1)*cos(omega2*time)*Phi_2(1);

%WHY MINUS???
theta_2 = 1/2*x0(1)*cos(omega1*time)*Phi_1(2)+1/2*x0(1)*cos(omega2*time)*Phi_2(2);
  
figure(1)
plot(y(:,1),y(:,2),'-k')
xlabel('\theta_1')
ylabel('\theta_2')
hold on
plot (theta_1,theta_2,'-r*')
legend('ode45','Normal modes')
title('Phase Space')

o = y(:,1);
o = o(1:10:end);

figure(2)
plot (time,theta_1,'-r*')
hold on
plot(T(1:10:end),o,'-ko')
xlabel('time')
ylabel('\theta_1')
legend('ode45','Normal modes')

o = y(:,2);
o = o(1:10:end);

figure(3)
plot(time,theta_2,'-r*')
hold on
plot(T(1:10:end),o,'-ko')
legend('Normal modes','ode45')
xlabel('time')
ylabel('\theta_2')


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