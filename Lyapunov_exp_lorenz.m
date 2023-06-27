
clear all

%Integrate the system first to get points closer to the attractor
%or use ICs below

[TT,XX] = ode45(@v,[0,10],[1;2;3]);

x0 = [XX(end,1);XX(end,1);XX(end,3);1;0;0;0;1;0;0;0;1];

% ICs by Sandri

%x0 = [19;20;50;1;0;0;0;1;0;0;0;1];

% Total integration time, time-step, #of iterations
T0 = 800;

dt = 0.1;
N = T0/dt;

% Initial region of initial conditions
w = eye(3);

% Initial sums for LCEs
sum1 = 0;
sum2 = 0;
sum3 = 0;

% Set the loop for LCEs calculation
for i = 1:N
    %Integrate the system for small dt
[jac,X] = J((i-1)*dt,i*dt,x0);

% Image of the vectors under the transformation
w = jac*w;
 
% %Reorthogonalize by GS-process
w1 = w(:,1);
e1 = w1/norm(w1);
w2 = w(:,2)-(e1'*w(:,2))*e1;
e2 = w2/norm(w2);
w3 = w(:,3)-(e1'*w(:,3))*e1-(e2'*w(:,3))*e2;
e3 = w3/norm(w3);
w = [e1,e2,e3];

%Sum up the separation 
w1 = log(norm(w1));
w2 = log(norm(w2));
w3 = log(norm(w3));

sum1(i+1) = sum1(i)+w1;
sum2(i+1) = sum2(i)+w2;
sum3(i+1) = sum3(i)+w3;



% Advance Ics for the next iteration
x0 = [X(:,end);w(1,:)';w(2,:)';w(3,:)'];%1;0;0;0;1;0;0;0;1];

    
end
 
% Compute the LCEs
exp1 = sum1/(T0);
exp2 = sum2/(T0);
exp3 = sum3/(T0);
t = linspace(1,N,N+1);
plot(t,exp1,t,exp2,t,exp3)
a = [exp1(end);exp2(end);exp3(end)]
% Function for calculation of the Jacobian(by Chaos book)/Monodromy
%/Flow matrix
function [jacobian,sol] = J(t0,t,X)

X(:,1) = X;

dt = 0.001;

N = (t-t0)/dt; 

t(1) = t0;

%Runge-kutta method
for i = 1:N
 
   k1 = diffjac(t(i),X(:,i));
   k2 = diffjac(t(i),X(:,i)+dt*k1/2);
   k3 = diffjac(t(i),X(:,i)+dt*k2/2);
   k4 = diffjac(t(i),X(:,i)+dt*k3);
   
   X(:,i+1) = X(:,i)+1/6*dt*(k1+2*k2+2*k3+k4);
   t(i+1) = t(i)+dt;
   


end

% The flow (monodromy) matrix (ChaosBook calls this jacobian) is the last 9 components of X
jac = X(:,end);
jac  = jac(4:12);
jacobian = zeros(3,3);

for i = 1:3
    
    jacobian(i,1:3) = jac((3*i-2):3*i);
    
end
% Define a variable for the trajectory f^(t-t0)(x0)
sol = X(1:3,:);


end

% Differental equation for Lorenz system (3 eqs) and the monodromy
%(3x3 => 9 eqs)
function jacobian = diffjac(t,X)

sig = 16;
rho = 45.92;
bet = 4;

dxdt1 = sig*(X(2)-X(1));
dxdt2 = X(1)*(rho-X(3))-X(2);
dxdt3 = X(1)*X(2)-bet*X(3);

% Actual Jacobian D(dx/dt)
A = [-sig,sig,0;rho-X(3),-1,-X(1);X(2),X(1),-bet];

% Diff.eq for the Monodromy matrix
jacobian = A*[X(4),X(5),X(6);
              X(7),X(8),X(9);
              X(10),X(11),X(12)];
              
jacobian = [dxdt1;dxdt2;dxdt3;jacobian(1,:)';jacobian(2,:)';jacobian(3,:)'];
          


end

function dxdt = v(t,X);

sig = 16;
rho = 45.92;
bet = 4;

dxdt1 = sig*(X(2)-X(1));
dxdt2 = X(1)*(rho-X(3))-X(2);
dxdt3 = X(1)*X(2)-bet*X(3);

dxdt = [dxdt1;dxdt2;dxdt3];

end