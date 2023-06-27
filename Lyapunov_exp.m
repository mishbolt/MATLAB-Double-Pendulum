
clear all

%Integrate the system first to get points closer to the attractor
%or use ICs below

%[TT,XX] = ode45(@v,[0,10],[1;2;3]);

%x0 = [XX(end,1);XX(end,1);XX(end,3);1;0;0;0;1;0;0;0;1];

% ICs


I = [1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;1];
x0 = [1.3;-1.3;0;0;I];

% Total integration time, time-step, #of iterations
T0 = 800;

dt = 0.1;
N = T0/dt;

% Initial region of initial conditions
w = eye(4);

% Initial sums for LCEs
sum1(1) = 0;
sum2(1) = 0;
sum3(1) = 0;
sum4(1) = 0;

options = odeset('RelTol',10^(-6),'AbsTol',10^(-8));
% Set the loop for LCEs calculation
for i = 1:N
    %Integrate the system for small dt

[~,X] = ode45(@diffjac,[(i-1)*dt,i*dt],x0, options);
jac = X(end,:)';
jac  = jac(5:20);
jacobian = zeros(4,4);

for j = 1:4
    
    jacobian(j,:) = jac((4*j-3):4*j);
    
end
% Image of the vectors under the transformation
     
w = jacobian*w;



% %Reorthogonalize by GS-process
w1 = w(:,1);
e1 = w1/norm(w1);
w2 = w(:,2)-(e1'*w(:,2))*e1;
e2 = w2/norm(w2);
w3 = w(:,3)-(e1'*w(:,3))*e1-(e2'*w(:,3))*e2;
e3 = w3/norm(w3);
w4 = w(:,4)-(e1'*w(:,4))*e1-(e2'*w(:,4))*e2-(e3'*w(:,4))*e3;
e4 = w4/norm(w4);
w = [e1,e2,e3,e4];

%Sum up the separation 
w1 = log(norm(w1));
w2 = log(norm(w2));
w3 = log(norm(w3));
w4 = log(norm(w4));

sum1(i+1) = (sum1(i)+w1);%/(dt*(i));
sum2(i+1) = (sum2(i)+w2);%/(dt*(i));
sum3(i+1)= (sum3(i)+w3);%/(dt*(i));
sum4(i+1) = (sum4(i)+w4);%/(dt*(i));



% Advance Ics for the next iteration
x0 = [X(end,1:4)';w(:,1);w(:,2);w(:,3);w(:,4);];

    

end
% Compute the LCEs
exp1 = sum1/T0;
exp2 = sum2/T0;
exp3 = sum3/T0;
exp4 = sum4/T0;
t = 0:dt:T0;
plot(t,exp1,t,exp2,t,exp3,t,exp4)
a = [exp1(end);exp2(end);exp3(end);exp4(end)]
%[~,A] = double_pend(x0);
%trace_A = trace(A)
% Function for calculation of the Jacobian(by Chaos book)/Monodromy
%/Flow matrix
% function [jacobian,sol] = J(t0,t,X)
% 
% X(:,1) = X;
% 
% dt = 0.001;
% 
% N = (t-t0)/dt; 
% 
% t(1) = t0;
% 
% %Runge-kutta method
% for i = 1:N
%  
%    k1 = diffjac(t(i),X(:,i));
%    k2 = diffjac(t(i),X(:,i)+dt*k1/2);
%    k3 = diffjac(t(i),X(:,i)+dt*k2/2);
%    k4 = diffjac(t(i),X(:,i)+dt*k3);
%    
%    X(:,i+1) = X(:,i)+1/6*dt*(k1+2*k2+2*k3+k4);
%    t(i+1) = t(i)+dt;
%    
% 
% 
% end
% 
% % The matrix is the last 9 components of X
% jac = X(:,end);
% jac  = jac(5:20);
% jacobian = zeros(4,4);
% 
% for i = 1:4
%     
%     jacobian(i,1:4) = jac((4*i-3):4*i);
%     
% end
% % Define a variable for the trajectory f^(t-t0)(x0)
% sol = X(1:4,:);
% 
% 
% end

% Differental equation for Lorenz system (4 eqs) and the monodromy
%(3X(3) => 16 eqs)
function [jacobian] = diffjac(t,X)

mu = 1;
lam = 1;


dxdt1 = (X(3)*1/lam-X(4)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);
dxdt2 = (lam*(mu+1)*X(4)-X(3)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);
dxdt3 = -(X(3)*X(4)*sin(X(1)-X(2)))/((mu+(sin(X(1)-X(2)))^2))+...
    (sin(2*(X(1)-X(2)))*(1/2*1/lam*(X(3))^2+1/2*lam*(mu+1)*(X(4))^2-X(3)*X(4)*cos(X(1)-X(2))))/...
    (mu+(sin(X(1)-X(2)))^2)^2-(mu+1)*sin(X(1));
dxdt4 = (X(3)*X(4)*sin(X(1)-X(2)))/((mu+(sin(X(1)-X(2)))^2))...
    -(sin(2*(X(1)-X(2)))*(1/2*1/lam*(X(3))^2+1/2*lam*(mu+1)*...
    (X(4))^2-X(3)*X(4)*cos(X(1)-X(2))))/(mu+(sin(X(1)-X(2)))^2)^2-1/lam*sin(X(2));


% Actual Jacobian D(dx/dt)
% Y1 = (der(t,[X(1)+h,X(2),X(3),X(4)])-der(t,[X(1)-h,X(2),X(3),X(4)]))/(2*h);
% Y2 = (der(t,[X(1),X(2)+h,X(3),X(4)])-der(t,[X(1),X(2)-h,X(3),X(4)]))/(2*h);
% Y3 = (der(t,[X(1),X(2),X(3)+h,X(4)])-der(t,[X(1),X(2),X(3)-h,X(4)]))/(2*h);
% Y4 = (der(t,[X(1),X(2),X(3),X(4)+h])-der(t,[X(1),X(2),X(3),X(4)-h]))/(2*h);
% 
% for i =1:4
%     
%     A(i,1) = Y1(i);
%     A(i,2) = Y2(i);
%     A(i,3) = Y3(i);
%     A(i,4) = Y4(i);
%     
% end
A = zeros(4,4);

A(1,1) = lam.^(-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*(lam.*(1+mu).* ...
  X(4)+(-2).*X(3).*cos(X(1)+(-1).*X(2))+lam.*X(4).*cos(X(1)+(-1).*X(2)).^2).*sin( ...
  X(1)+(-1).*X(2));
A(2,1) = (1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*((1+mu).*X(3)+(-2).*lam.*(1+ ...
  mu).*X(4).*cos(X(1)+(-1).*X(2))+X(3).*cos(X(1)+(-1).*X(2)).^2).*sin(X(1)+(-1).* ...
  X(2));
A(3,1) = lam.^(-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-3).*((-4).*cos(X(1)+( ...
  -1).*X(2)).*sin(X(1)+(-1).*X(2)).*((-1/4).*lam.*(1+mu).*((-1)+(-2).*mu+ ...
  cos(2.*(X(1)+(-1).*X(2)))).^2.*sin(X(1))+((-1).*lam.*(1+mu).*X(4)+X(3).*cos( ...
  X(1)+(-1).*X(2))).*(X(3)+(-1).*lam.*X(4).*cos(X(1)+(-1).*X(2))).*sin(X(1)+(-1).* ...
  X(2)))+(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).*(cos(X(1)+(-1).*X(2)).*((-1).* ...
  lam.*(1+mu).*X(4)+X(3).*cos(X(1)+(-1).*X(2))).*(X(3)+(-1).*lam.*X(4).*cos(X(1)+( ...
  -1).*X(2)))+(-1/4).*lam.*(1+mu).*cos(X(1)).*((-1)+(-2).*mu+cos(2.*(X(1)+ ...
  (-1).*X(2)))).^2+lam.*X(4).*((-1).*lam.*(1+mu).*X(4)+X(3).*cos(X(1)+(-1).* ...
  X(2))).*sin(X(1)+(-1).*X(2)).^2+X(3).*((-1).*X(3)+lam.*X(4).*cos(X(1)+(-1).*X(2))) ...
  .*sin(X(1)+(-1).*X(2)).^2+lam.*(1+mu).*((-1)+(-2).*mu+cos(2.*(X(1)+(-1) ...
  .*X(2)))).*sin(X(1)).*sin(2.*(X(1)+(-1).*X(2)))));

A(4,1) = lam.^(-1).*(1+2.*mu+(-1).*cos(2.*(X(1)+(-1).*X(2)))).^(-3).*(6.*X(3).^2+ ...
  6.*lam.^2.*X(4).^2+6.*lam.^2.*mu.*X(4).^2+lam.*((-11)+4.*mu+8.*mu.^2) ...
  .*X(3).*X(4).*cos(X(1)+(-1).*X(2))+(-4).*(1+2.*mu).*(X(3).^2+lam.^2.*(1+mu) ...
  .*X(4).^2).*cos(2.*(X(1)+(-1).*X(2)))+(21/2).*lam.*X(3).*X(4).*cos(3.*(X(1)+( ...
  -1).*X(2)))+12.*lam.*mu.*X(3).*X(4).*cos(3.*(X(1)+(-1).*X(2)))+(-2).*X(3).^2.* ...
  cos(4.*(X(1)+(-1).*X(2)))+(-2).*lam.^2.*X(4).^2.*cos(4.*(X(1)+(-1).*X(2)))+( ...
  -2).*lam.^2.*mu.*X(4).^2.*cos(4.*(X(1)+(-1).*X(2)))+(1/2).*lam.*X(3).*X(4).* ...
 cos(5.*(X(1)+(-1).*X(2))));

A(1,2) = (-1).*lam.^(-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*(lam.*(1+ ...
  mu).*X(4)+(-2).*X(3).*cos(X(1)+(-1).*X(2))+lam.*X(4).*cos(X(1)+(-1).*X(2)).^2).* ...
  sin(X(1)+(-1).*X(2));

A(2,2) = (-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*((1+mu).*X(3)+(-2).* ...
  lam.*(1+mu).*X(4).*cos(X(1)+(-1).*X(2))+X(3).*cos(X(1)+(-1).*X(2)).^2).*sin( ...
  X(1)+(-1).*X(2));

A(3,2) = lam.^(-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-3).*(4.*cos(X(1)+(-1) ...
  .*X(2)).*sin(X(1)+(-1).*X(2)).*((-1/4).*lam.*(1+mu).*((-1)+(-2).*mu+cos( ...
  2.*(X(1)+(-1).*X(2)))).^2.*sin(X(1))+((-1).*lam.*(1+mu).*X(4)+X(3).*cos(X(1)+( ...
  -1).*X(2))).*(X(3)+(-1).*lam.*X(4).*cos(X(1)+(-1).*X(2))).*sin(X(1)+(-1).*X(2))) ...
  +(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).*((-1).*(X(3).^2+lam.^2.*(1+mu).* ...
  X(4).^2).*cos(X(1)+(-1).*X(2)).^2+lam.*X(3).*X(4).*cos(X(1)+(-1).*X(2)).^3+lam.* ...
  X(3).*X(4).*cos(X(1)+(-1).*X(2)).*(mu+cos(2.*(X(1)+(-1).*X(2))))+(X(3).^2+ ...
  lam.^2.*(1+mu).*X(4).^2).*sin(X(1)+(-1).*X(2)).^2+lam.*(1+mu).*(1+2.*mu+ ...
  (-1).*cos(2.*(X(1)+(-1).*X(2)))).*sin(X(1)).*sin(2.*(X(1)+(-1).*X(2)))));

A(4,2) = lam.^(-1).*((-1)+(-2).*mu+cos(2.*(X(1)+(-1).*X(2)))).^(-3).*((-4).* ...
  sin(2.*(X(1)+(-1).*X(2))).*(2.*lam.*X(3).*X(4).*(3+2.*mu+cos(2.*(X(1)+(-1).* ...
  X(2)))).*sin(X(1)+(-1).*X(2))+(-2).*(X(3).^2+lam.^2.*(1+mu).*X(4).^2).*sin( ...
  2.*(X(1)+(-1).*X(2)))+(-1).*((-1)+(-2).*mu+cos(2.*(X(1)+(-1).*X(2)))).^2.* ...
  sin(X(2)))+((-1)+(-2).*mu+cos(2.*(X(1)+(-1).*X(2)))).*(4.*(X(3).^2+ ...
  lam.^2.*(1+mu).*X(4).^2).*cos(2.*(X(1)+(-1).*X(2)))+(-2).*lam.*X(3).*X(4).* ...
  cos(X(1)+(-1).*X(2)).*(3+2.*mu+cos(2.*(X(1)+(-1).*X(2))))+(-1).*((-1)+(-2) ...
  .*mu+cos(2.*(X(1)+(-1).*X(2)))).^2.*cos(X(2))+8.*lam.*X(3).*X(4).*cos(X(1)+( ...
  -1).*X(2)).*sin(X(1)+(-1).*X(2)).^2+(-4).*((-1)+(-2).*mu+cos(2.*(X(1)+(-1) ...
  .*X(2)))).*sin(2.*(X(1)+(-1).*X(2))).*sin(X(2))));

A(1,3) = (lam+lam.*mu+(-1).*lam.*cos(X(1)+(-1).*X(2)).^2).^(-1);

A(2,3) = cos(X(1)+(-1).*X(2)).*((-1)+(-1).*mu+cos(X(1)+(-1).*X(2)).^2).^(-1);

A(3,3) = (-1).*lam.^(-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*(lam.*(1+ ...
  mu).*X(4)+(-2).*X(3).*cos(X(1)+(-1).*X(2))+lam.*X(4).*cos(X(1)+(-1).*X(2)).^2).* ...
  sin(X(1)+(-1).*X(2));

A(4,3) = 2.*lam.^(-1).*((-1)+(-2).*mu+cos(2.*(X(1)+(-1).*X(2)))).^(-2).*((-4).* ...
  X(3).*cos(X(1)+(-1).*X(2))+lam.*X(4).*(3+2.*mu+cos(2.*(X(1)+(-1).*X(2))))).* ...
  sin(X(1)+(-1).*X(2));

A(1,4) = cos(X(1)+(-1).*X(2)).*((-1)+(-1).*mu+cos(X(1)+(-1).*X(2)).^2).^(-1);
A(2,4) = lam.*(1+mu).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-1);

A(3,4) = (-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*((1+mu).*X(3)+(-2).* ...
  lam.*(1+mu).*X(4).*cos(X(1)+(-1).*X(2))+X(3).*cos(X(1)+(-1).*X(2)).^2).*sin( ...
  X(1)+(-1).*X(2));

A(4,4) = 2.*((-1)+(-2).*mu+cos(2.*(X(1)+(-1).*X(2)))).^(-2).*((-4).*lam.*(1+mu) ...
  .*X(4).*cos(X(1)+(-1).*X(2))+X(3).*(3+2.*mu+cos(2.*(X(1)+(-1).*X(2))))).*sin( ...
  X(1)+(-1).*X(2));

% Diff.eq for the Monodromy matrix
jacobian = A*[X(5),X(6),X(7),X(8);
              X(9),X(10),X(11),X(12);
              X(13),X(14),X(15),X(16);
              X(17),X(18),X(19),X(20)];
              
jacobian = [dxdt1;dxdt2;dxdt3;dxdt4;jacobian(1,:)';jacobian(2,:)';jacobian(3,:)';jacobian(4,:)'];
       


end

function dxdt = der(t,X);

mu = 1;
lam = 1;

dxdt1 = (X(3)*1/lam-X(4)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);
dxdt2 = (lam*(mu+1)*X(4)-X(3)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);
dxdt3 = -(X(3)*X(4)*sin(X(1)-X(2)))/((mu+(sin(X(1)-X(2)))^2))+...
    (sin(2*(X(1)-X(2)))*(1/2*1/lam*(X(3))^2+1/2*lam*(mu+1)*(X(4))^2-X(3)*X(4)*cos(X(1)-X(2))))/...
    (mu+(sin(X(1)-X(2)))^2)^2-(mu+1)*sin(X(1));
dxdt4 = (X(3)*X(4)*sin(X(1)-X(2)))/((mu+(sin(X(1)-X(2)))^2))...
    -(sin(2*(X(1)-X(2)))*(1/2*1/lam*(X(3))^2+1/2*lam*(mu+1)*...
    (X(4))^2-X(3)*X(4)*cos(X(1)-X(2))))/(mu+(sin(X(1)-X(2)))^2)^2-1/lam*sin(X(2));

dxdt = [dxdt1;dxdt2;dxdt3;dxdt4];

end