
m1 = 10;
m2 = 1;
L1 = 1;
L2 = L1;
L=L1;
g = 9.8;

tf = 100;

time = 0:0.01:tf;

x0 = [pi/20;0;0;0];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[T,y] = ode45(@derivative,[0,tf],x0,options);


%Eigenfrquencies

omega1 = sqrt(2+sqrt(2));
omega2 = sqrt(2-sqrt(2));

%------Normalization Constants and Mode Amplitudes------%

%Generalized eigenvectors
Phi_1 = [1;-sqrt(2)];
Phi_2 = [1;sqrt(2)];

%t-matrix
T_matrix = 1/2*[2,1;
                1,1];
      
%Normalization Constants     
N1 = sqrt(1./(Phi_1'*T_matrix*Phi_1));
N2 = sqrt(1./(Phi_2'*T_matrix*Phi_2));

%Mode Amplitudes
A1 = N1*Phi_1'*T_matrix*[x0(1);x0(2)];
A2 = N2*Phi_2'*T_matrix*[x0(1);x0(2)];


%Solutions
theta_1 = A1*N1*cos(omega1*time)*Phi_1(1)+A2*N2*cos(omega2*time)*Phi_2(1);

%WHY MINUS???
theta_2 = A1*N1*cos(omega1*time)*Phi_1(2)+A2*N2*cos(omega2*time)*Phi_2(2);
  
figure(1)
plot(y(:,1),y(:,2),'-k')
xlabel('\theta_1')
ylabel('\theta_2')
hold on
%plot (theta_1,theta_2,'-r*')
legend('ode45','Normal modes')
title('Phase Space')

figure(2)
plot(T,y(:,1),'-ko')
hold on
plot (time,theta_1,'-r*')
xlabel('time')
ylabel('\theta_1')
legend('ode45','Normal modes')

figure(3)
plot(T,y(:,2),'-ko')
hold on
plot(time,theta_2,'-r*')
legend('ode45','Normal modes')
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