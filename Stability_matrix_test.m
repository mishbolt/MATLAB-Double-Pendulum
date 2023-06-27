clc
clear all

p2 = (0.02557+0.02948)/2;
E0 = -2.99954308546917;
p1 = 1/2*(2*p2+sqrt(4*p2^2+8*(3+E0-p2^2)));

x0 = [0;0;p1;p2];
T0 = 8.18537366184949;

stab(x0)
jac = J(T0,x0);
det(jac)
[EIG,EIGV] = eig(jac);
A = stab(x0);
eig(A)
omega = [zeros(2,2),eye(2);
         -eye(2),zeros(2,2)];
iszero1 = A'*omega+omega*A;
isomega = jac'*omega*jac;

velocity = jac*der(x0);

function A = stab(X)

h = 1e-9;

Y1 = (der([X(1)+h,X(2),X(3),X(4)])-der([X(1),X(2),X(3),X(4)]))/h;
Y2 = (der([X(1),X(2)+h,X(3),X(4)])-der([X(1),X(2),X(3),X(4)]))/h;
Y3 = (der([X(1),X(2),X(3)+h,X(4)])-der([X(1),X(2),X(3),X(4)]))/h;
Y4 = (der([X(1),X(2),X(3),X(4)+h])-der([X(1),X(2),X(3),X(4)]))/h;

for i =1:4
    
        A(i,1) = Y1(i);
        A(i,2) = Y2(i);
        A(i,3) = Y3(i);
        A(i,4) = Y4(i);
    
end


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

A = stab(X);


jac0 = [1;0;0;0;0;1;0;0;0;0;1;0;0;0;0;1];
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

[~,jac] = ode45(@diffjac,[0,20],jac0,options,A);


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