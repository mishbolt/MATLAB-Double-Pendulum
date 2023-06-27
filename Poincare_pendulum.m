
clear all

q1(1) = -pi/180; p1(1) = 0;
q2(1) = 0; p2(1) = 0;

E0 = H([q1(1),q2(1),0,0]);


% solve for set of ICs with the same energy
N=15;
q1 = linspace(-pi/180,pi/180,N);

for k = 1:N
q2(k) = fsolve(@(q2) H2(q1(k),q2)-E0,pi/180);
end

options = odeset('Events',@pendevents,'RelTol',1e-6,'AbsTol',1e-10);
% Set the loop for various initial condtions
for i=1:N
% Set initial parameters
m1 = 1;
m2 = 1;
L1 = 1;
L2 = 1;
g = 9.8;
tf = 1000;
x0 = [q1(i);q2(i);0;0];

[T,x,te,ye,ie] = ode45(@derivative,[0,tf],x0,options);

plot(ye(:,2),ye(:,4),'r.')
xlabel('\beta'); ylabel('p_\beta')
hold on

end

% Function for ode45
function dxdt = derivative(t,X)
m1 = 1;
m2 = 1;
L1 = 1;
L2 = 1;
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
function ham = H(X)

m1 = 1;
m2 = 1;
L1 = 1;
L2 = 1;
g = 9.8;

T = (((X(3))^2*L2^2*m2+(X(4))^2*(m1+m2)*L1^2-m2*L1*L2*X(3)*X(4)*cos(X(1)-X(2))))/...
    (2*L1^2*L2^2*(m1+m2*(sin(X(1)-X(2)))^2));

U = -g*L1*(m1+m2)*cos(X(1))-g*m2*L2*cos(X(2));


ham = T+U;



end

function ham = H2(q1,q2)

m1 = 1;
m2 = 1;
L1 = 1;
L2 = 1;
g = 9.8;
p1 = 0;
p2 = 0;

T = (((p1)^2*L2^2*m2+(p2)^2*(m1+m2)*L1^2-m2*L1*L2*p1*p2*cos(q1-q2)))./...
    (2*L1^2*L2^2*(m1+m2*(sin(q1-q2)).^2));

U = -g*L1*(m1+m2)*cos(q1)-g*m2*L2*cos(q2);


ham = T+U;



end

% Function for events
function [value,isterminal,direction]=pendevents(t,y)

value = y(1);
isterminal = 0;
direction = 1;


end