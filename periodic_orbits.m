
clc
clear all

p2 = (0.02557+0.02948)/2;
E0 = -2.99954308546917;
p1 = 1/2*(2*p2+sqrt(4*p2^2+8*(3+E0-p2^2)));

N =1 ;
x0 = [7.08175230799863e-18;4.85314520672181e-18;0.685147016443462;0.430000000000000];
T0 = [8.12091254489115];
err = 1;
func = (f(T0,x0))';


while err>=1e-4 % WHILE abs(f(x)-x)>tol
    
    velocity = -der(func(:,end));
    jacobian = J(T0,x0);
    A = [eye(4)-jacobian,-der(func(:,end));
        der(func(:,end))',               0];
    
    
    Y = -[x0-(func(:,end))  ;0];
    
    
    XT = A\Y;
    
  [Eigv, Eig] = eig(J(T0,x0));
    
   x0 = x0+XT(1:4);
   T0 = T0+XT(5);
   
   
    func = (f(T0,x0))';
    
    err = abs(func(:,end)-x0);
    err = norm(err);
   
    N = N+1;
    
end

perorb = f(T0,x0);

det(J(T0,x0))


function x = f(t,x0)

options = [];

[~,x] = ode45(@derivative,[0,t],x0,options);


end

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


function A = stab(t,X)

h = 1e-9;
X0 = X;

[~,X] = ode45(@derivative,[0,t],X0);

X = X(:,end);

Y1 = (derivative(t,[X(1)+h,X(2),X(3),X(4)])-derivative(t,[X(1),X(2),X(3),X(4)]))/h;
Y2 = (derivative(t,[X(1),X(2)+h,X(3),X(4)])-derivative(t,[X(1),X(2),X(3),X(4)]))/h;
Y3 = (derivative(t,[X(1),X(2),X(3)+h,X(4)])-derivative(t,[X(1),X(2),X(3),X(4)]))/h;
Y4 = (derivative(t,[X(1),X(2),X(3),X(4)+h])-derivative(t,[X(1),X(2),X(3),X(4)]))/h;

A(:,1) = Y1;
A(:,2) = Y2;
A(:,3) = Y3;
A(:,4) = Y4;


end



function dxdt = der(X)
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


function jacobian = J(t,X)




jac0(:,1) = [1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;1];
options = odeset('RelTol',1e-9,'AbsTol',1e-8);


N = 1e3;

dt = t/N;

t0 = 0;

for i = 1:N
    if t0==0
        h = 1e-9;
        Y1 = (derivative(t,[X(1)+h,X(2),X(3),X(4)])-derivative(t,[X(1),X(2),X(3),X(4)]))/h;
        Y2 = (derivative(t,[X(1),X(2)+h,X(3),X(4)])-derivative(t,[X(1),X(2),X(3),X(4)]))/h;
        Y3 = (derivative(t,[X(1),X(2),X(3)+h,X(4)])-derivative(t,[X(1),X(2),X(3),X(4)]))/h;
        Y4 = (derivative(t,[X(1),X(2),X(3),X(4)+h])-derivative(t,[X(1),X(2),X(3),X(4)]))/h;
        
        for i =1:4
            
            A(i,1) = Y1(i);
            A(i,2) = Y2(i);
            A(i,3) = Y3(i);
            A(i,4) = Y4(i);
            
        end
    else
        A = stab(t0,X);
    end
    [~,jac] = ode45(@diffjac,[t0,t0+dt],jac0,options,A);
    
    t0 = t0+dt;
    jac0 = jac(end,:)';
end

jac = jac(end,:)';
jacobian = zeros(4,4);

for i = 1:4
    
    jacobian(i,1:4) = jac((4*i-3):4*i);
    
end


end
function jacobian = diffjac(t,X,A)


jacobian = [A(1,1)...
    *X(1)+A(1,2)*X(5)+A(1,3)*X(9)+A(1,4)*X(13);
    A(1,1)*X(2)+A(1,2)*X(6)+A(1,3)*X(10)+A(1,4)*X(14);
    A(1,1)*X(3)+A(1,2)*X(7)+A(1,3)*X(11)+A(1,4)*X(15);
    A(1,1)*X(4)+A(1,2)*X(8)+A(1,3)*X(12)+A(1,4)*X(16);
    A(2,1)*X(1)+A(2,2)*X(5)+A(2,3)*X(9)+A(2,4)*X(13);
    A(2,1)*X(2)+A(2,2)*X(6)+A(2,3)*X(10)+A(2,4)*X(14);
    A(2,1)*X(3)+A(2,2)*X(7)+A(2,3)*X(11)+A(2,4)*X(15);
    A(2,1)*X(4)+A(2,2)*X(8)+A(2,3)*X(12)+A(2,4)*X(16);
    A(3,1)*X(1)+A(3,2)*X(5)+A(3,3)*X(9)+A(3,4)*X(13);
    A(3,1)*X(2)+A(3,2)*X(6)+A(3,3)*X(10)+A(3,4)*X(14);
    A(3,1)*X(3)+A(3,2)*X(7)+A(3,3)*X(11)+A(3,4)*X(15);
    A(3,1)*X(4)+A(3,2)*X(8)+A(3,3)*X(12)+A(3,4)*X(16);
    A(4,1)*X(1)+A(4,2)*X(5)+A(4,3)*X(9)+A(4,4)*X(13);
    A(4,1)*X(2)+A(4,2)*X(6)+A(4,3)*X(10)+A(4,4)*X(14);
    A(4,1)*X(3)+A(4,2)*X(7)+A(4,3)*X(11)+A(4,4)*X(15);
    A(4,1)*X(4)+A(4,2)*X(8)+A(4,3)*X(12)+A(4,4)*X(16);];




end

% function dxdt = newder(X)
% mu = 1;
% lam = 1;
% 
% dxdt1 = (X(3)*1/lam-X(4)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);
% dxdt2 = (lam*(mu+1)*X(4)-X(3)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);
% dxdt3 = -(X(3)*X(4)*sin(X(1)-X(2)))/((mu+(sin(X(1)-X(2)))^2))+...
%     (sin(2*(X(1)-X(2)))*(1/2*1/lam*(X(3))^2+1/2*lam*(mu+1)*(X(4))^2-X(3)*X(4)*cos(X(1)-X(2))))/...
%     (mu+(sin(X(1)-X(2)))^2)^2-(mu+1)*sin(X(1));
% dxdt4 = (X(3)*X(4)*sin(X(1)-X(2)))/((mu+(sin(X(1)-X(2)))^2))...
%     -(sin(2*(X(1)-X(2)))*(1/2*1/lam*(X(3))^2+1/2*lam*(mu+1)*...
%     (X(4))^2-X(3)*X(4)*cos(X(1)-X(2))))/(mu+(sin(X(1)-X(2)))^2)^2-1/lam*sin(X(2));
% 
% dxdt = [-dxdt3;-dxdt4;dxdt1;dxdt2];
% 
% end



