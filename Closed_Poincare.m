clc
clear



x0 = [-pi/180;pi/180;0;0];

k = 0:100;

% theta_1_0(1) = -pi/180;
% theta_2_0(1) = -pi/180;
% 
% 
% E0 = 
% theta_1_0 = linspace(-pi/180,pi/180);


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


for i = 1:20
    
    
    
    
end
%Mode Amplitudes
A1 = N1*Phi_1'*T_matrix*[x0(1);x0(2)];
A2 = N2*Phi_2'*T_matrix*[x0(1);x0(2)];

C1 = 1/2*(x0(1)-sqrt(2)/2*x0(2));
C2 = 1/2*(x0(1)+sqrt(2)/2*x0(2));


time = real(1/(omega2-omega1)*acos((C1^1+C2^2)/(2*C1*C2))+2*pi*k/(omega2-omega1));
t = linspace(0,10);
theta_1 = C1*cos(omega1*time)+C2*cos(omega2*time);
theta_2 = C1*(-sqrt(2))*cos(omega1*time)+C2*(sqrt(2))*cos(omega2*time);

theta_2_d = -omega1*C1*(-sqrt(2))*sin(omega1*time)-C2*omega2*(sqrt(2))*sin(omega2*time);

theta_1_test = @(t) C1*cos(omega1*t)+C2*cos(omega2*t);
%Solutions

%plot(t,theta_1_test,'r.')

tt = fsolve(theta_1_test,9);
function ham = E(X)
mu = 1;
lam = 1;

ham = 1/2*(mu+1)*(X(1))^2+1/2*lam^(-1)*(X(2))^2;



end




