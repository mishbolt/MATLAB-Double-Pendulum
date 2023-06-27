
m1 = 10;
m2 = 1;
L1 = 1;
L2 = L1;
L=L1;
g = 9.8;

tf = 10;

time = 0:0.01:tf;

x0 = [pi/180;0;0;0];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[T,y] = ode45(@derivative,[0,tf],x0,options);


%Eigenfrquencies

omega1 = sqrt(g/(m1*L)*(m1+m2+sqrt((m1+m2)*m2)));
omega2 = sqrt(g/(m1*L)*(m1+m2-sqrt((m1+m2)*m2)));

%------Normalization Constants and Mode Amplitudes------%

%Generalized eigenvectors
Phi_1 = [1/sqrt(m1+m2);1/sqrt(m2)];
Phi_2 = [1/sqrt(m1+m2);-1/sqrt(m2)];

%t-matrix
T_matrix = 1/2*[(m1+m2)*L^2,m2*L^2;
                 m2*L^2  ,  m2*L^2];
      
%Normalization Constants     
N1 = sqrt(1./(Phi_1'*T_matrix*Phi_1));
N2 = sqrt(1./(Phi_2'*T_matrix*Phi_2));

%Mode Amplitudes
A1 = N1*Phi_1'*T_matrix*[x0(1);x0(2)];
A2 = N2*Phi_2'*T_matrix*[x0(1);x0(2)];


%Solutions
theta_1 = A1*N1*cos(omega1*time)*Phi_1(1)+A2*N2*cos(omega2*time)*Phi_2(1);

%WHY MINUS???
theta_2 = -A1*N1*cos(omega1*time)*Phi_1(2)-A2*N2*cos(omega2*time)*Phi_2(2);
  
figure(1)
plot(y(:,1),y(:,2),'-k')
xlabel('\theta_1')
ylabel('\theta_2')
hold on
plot (theta_1,theta_2,'-r*')
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
m1 = 10;
m2 = 1;
L1 = 1;
L2 = L1;
g = 9.8;
A = (X(3)*X(4)*sin(X(1)-X(2)))/(L1*L2*(m1+m2*(sin(X(1)-X(2)))^2));

B = (((X(3))^2*L2^2*m2+(X(4))^2*(m1+m2)*L1^2-m2*L1*L2*X(3)*X(4)*cos(X(1)-X(2))))...
    *sin(2*(X(1)-X(2)))/(2*L1^2*L2^2*(m1+m2*(sin(X(1)-X(2)))^2)^2);

dxdt1 = (X(3)*L2-X(4)*L1*cos(X(1)-X(2)))/(L1^2*L2*(m1+m2*(sin(X(1)-X(2)))^2));
dxdt2 = (X(4)*(m1+m2)*L1-X(3)*m2*L2*cos(X(1)-X(2)))/(L1*L2^2*m2*(m1+m2*(sin(X(1)-X(2)))^2));
dxdt3 =  B-A-g*L1*(m1+m2)*sin(X(1));
dxdt4 =  A-B-g*m2*L2*sin(X(2));

dxdt = [dxdt1;dxdt2;dxdt3;dxdt4];
end