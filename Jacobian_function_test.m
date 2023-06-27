 
clear all

p2 = (0.02557+0.02948)/2;
E0 = -2.99954308546917;
p1 = 1/2*(2*p2+sqrt(4*p2^2+8*(3+E0-p2^2)));

x0 =[-1.954144432525574e-04;
    -2.766601233358286e-04;
     3.883865839844799e-02;
     2.746384754317892e-02];

T0 = 8.18537366184949;

JAC = J(T0,x0)
det(JAC)


function jacobian = J(t,X)

h0 = 1e-7;
options = [];
jac0 = X;

[~,f0] = ode45(@diffjac,[0,t],jac0,options,[0,0,0,0]);

[~,f1] = ode45(@diffjac,[0,t],jac0,options,[h0,0,0,0]);
[~,f2] = ode45(@diffjac,[0,t],jac0,options,[0,h0,0,0]);
[~,f3] = ode45(@diffjac,[0,t],jac0,options,[0,0,h0,0]);
[~,f4] = ode45(@diffjac,[0,t],jac0,options,[0,0,0,h0]);

f0 = f0(end,:)';
f1 = f1(end,:)';
f2 = f2(end,:)';
f3 = f3(end,:)';
f4 = f4(end,:)';

jacobian(:,1) = (f1-f0)/h0;
jacobian(:,2) = (f2-f0)/h0;
jacobian(:,3) = (f3-f0)/h0;
jacobian(:,4) = (f4-f0)/h0;


end


function dxdt = diffjac(t,X,h)
mu = 1;
lam = 1;

dxdt1 = ((X(3)+h(3))*1/lam-(X(4)+h(4))*cos((X(1)+h(1))-(X(2)+h(2))))/(mu+(sin((X(1)+h(1))-(X(2)+h(2))))^2);

dxdt2 = (lam*(mu+1)*(X(4)+h(4))-(X(3)+h(3))*cos((X(1)+h(1))-(X(2)+h(2))))/(mu+(sin((X(1)+h(1))-(X(2)+h(2))))^2);

A = -((X(3)+h(3))*(X(4)+h(4))*sin((X(1)+h(1))-(X(2)+h(2))))/((mu+(sin((X(1)+h(1))-(X(2)+h(2))))^2));
B = (sin(2*((X(1)+h(1))-(X(2)+h(2))))*(1/2*1/lam*((X(3)+h(3)))^2+1/2*lam*(mu+1)*((X(4)+h(4)))^2-(X(3)+h(3))*(X(4)+h(4))*cos((X(1)+h(1))-(X(2)+h(2)))))/(mu+(sin((X(1)+h(1))-(X(2)+h(2))))^2)^2;

dxdt3 = A+B-(mu+1)*sin((X(1)+h(1)));

dxdt4 = -A-B-1/lam*sin((X(2)+h(2)));

dxdt = [dxdt1;dxdt2;dxdt3;dxdt4];



end