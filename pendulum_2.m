
clear all

% Set initial parameters
m1 = 1;
m2 = 1;
L1 = 1;
L2 = 1;
g = 9.8;
qold1 = pi/100; pold1 = 0;qold2 =pi/100 ;pold2 = 0;
tf = 30;
N = 200000;
h = tf/N;
x0 = [qold1;qold2;pold1;pold2];
x = zeros(4,N+1);
x(:,1) = x0;
error = 1;
tol = 1e-12;
t = 0:h:tf;
options = [];

Ham(1) = H(t,x0);


for k=1:N
    
    % Verlet+ Newton for p1
    pguess1 = pold1;
    error = 1;
    while error>tol
        
        [~,~,f,~] = deriv(t,[qold1,qold2,pguess1,pold2]);
        
        
        
        ptemp1 = pguess1-(pold1 + h/2*f-pguess1)/(h/2*derivp1(t,[qold1,qold2,pguess1,pold2])-1);
        
        error = abs((ptemp1-pguess1)/ptemp1);
        
        pguess1 = ptemp1;
        
        
    end
    qguess1 = qold1;
    error = 1;
    
    % Verlet+ Newton for q1
    while error>tol
        [f1,~,~,~] = deriv(t,[qguess1,qold2,ptemp1,pold2]);
        [f2,~,~,~] = deriv(t,[qold1,qold2,ptemp1,pold2]);
        
        qnew1 = qguess1-(qold1+h/2*(f1+f2)-qguess1)/(h/2*derivt1(t,[qguess1,qold2,ptemp1,pold2])-1);
        
        error = abs((qnew1-qguess1)/qnew1);
        
        qguess1 = qnew1;
        
    end
    [~,~,f,~] = deriv(t,[qnew1,qold2,ptemp1,pold2]);
    
    pnew1 = ptemp1+h/2*f;
    
    x(1,k+1) = qnew1;
    x(3,k+1) = pnew1;
    
    
    %**************************************************
    
    % Verlet+ Newton for p2
    
    pguess2 = pold2;
    error = 1;
    while error>tol
        
        [~,~,~,f] = deriv(t,[qold1,qold2,pold1,pguess2]);
        
        
        
        ptemp2 = pguess2-(pold2 + h/2*f-pguess2)/(h/2*derivp2(t,[qold1,qold2,pold1,pguess2])-1);
        
        error = abs((ptemp2-pguess2)/ptemp2);
        
        pguess2 = ptemp2;
        
        
    end
    qguess2 = qold2;
    error = 1;
    
    % Verlet+ Newton for q2
    while error>tol
        [~,f1,~,~] = deriv(t,[qold1,qguess2,pold1,ptemp2]);
        [~,f2,~,~] = deriv(t,[qold1,qold2,pold1,ptemp2]);
        
        qnew2 = qguess2-(qold2+h/2*(f1+f2)-qguess2)/(h/2*derivt2(t,[qold1,qguess2,pold1,ptemp2])-1);
        
        error = abs((qnew2-qguess2)/qnew2);
        
        qguess2 = qnew2;
        
    end
    [~,~,~,f] = deriv(t,[qold1,qnew2,pold1,ptemp2]);
    
    pnew2 = ptemp2+h/2*f;
    
    x(2,k+1) = qnew2;
    x(4,k+1) = pnew2;
    
    
    % Calculate Hamiltonian for new values
    Ham(k+1) = H(t,[qnew1,qnew2,pnew1,pnew2]);
    
    
    pold1 = pnew1;
    qold1 = qnew1;
    pold2 = pnew2;
    qold2 = qnew2;
    
    
    
end



% ode45
[T,y] = ode45(@derivative,[0,tf],x0);


% Plot results+compare with ode45
figure
plot(t,x(1,:),'-ro')
hold on
plot(T,y(:,1),'-ko')

xlabel('time')
ylabel('\theta_1')
legend('Implicit Verlet','ode45')

figure
plot(x(2,:),x(4,:),'-r')
hold on
plot(y(:,2),y(:,4),'-k')
xlabel('\theta_2')
ylabel('p_2')
legend('Implicit Verlet','ode45')


figure
plot(x(1,:),x(2,:),'-r')
hold on
plot(y(:,1),y(:,2),'-k')
xlabel('\theta_1')
ylabel('\theta_2')
legend('Implicit Verlet','ode45')

 
figure
plot(x(1,:),x(3,:),'-r')
hold on
plot(y(:,1),y(:,3),'-k')
xlabel('\theta_1')
ylabel('p_1')
legend('Implicit Verlet','ode45')


figure
plot(x(3,:),x(4,:),'-r')
hold on
plot(y(:,3),y(:,4),'-k')
xlabel('p_1')
ylabel('p_2')
legend('Implicit Verlet','ode45')

% Plot total energy aka Hamiltonian
figure
plot (t,Ham,'-r')

xlabel('time')
ylabel('total energy H = T+V')
legend('Implicit Verlet')

% Function for the diff eqs
function [dxdt1,dxdt2,dxdt3,dxdt4] = deriv(t,X)
m1 = 1;
m2 = 0.6;
L1 = 1;
L2 = 0.5;
g = 9.8;
A = (X(3)*X(4)*sin(X(1)-X(2)))/(L1*L2*(m1+m2*(sin(X(1)-X(2)))^2));

B = (((X(3))^2*L2^2*m2+(X(4))^2*(m1+m2)*L1^2-m2*L1*L2*X(3)*X(4)*cos(X(1)-X(2))))...
    *sin(2*(X(1)-X(2)))/(2*L1^2*L2^2*(m1+m2*(sin(X(1)-X(2)))^2)^2);

dxdt1 = (X(3)*L2-X(4)*L1*cos(X(1)-X(2)))/(L1^2*L2*(m1+m2*(sin(X(1)-X(2)))^2));
dxdt2 = (X(4)*(m1+m2)*L1-X(3)*m2*L2*cos(X(1)-X(2)))/(L1*L2^2*m2*(m1+m2*(sin(X(1)-X(2)))^2));
dxdt3 =  B-A-g*L1*(m1+m2)*sin(X(1));
dxdt4 =  A-B-g*m2*L2*sin(X(2));

end



% Derivatives for the Newton's Method at each iteration
function dxdt = derivp1(t,X)


h = 1e-8;
[~,~,f1,~] = deriv(t,[X(1),X(2),X(3)+h,X(4)]);
[~,~,f2,~] = deriv(t,[X(1),X(2),X(3)-h,X(4)]);




dxdt = (f1-f2)/(2*h);

end

function dxdt = derivp2(t,X)


h = 1e-8;

[~,~,~,f1] = deriv(t,[X(1),X(2),X(3),X(4)+h]);
[~,~,~,f2] = deriv(t,[X(1),X(2),X(3),X(4)-h]);


dxdt = (f1-f2)/(2*h);


end

function dxdt = derivt1(t,X)


h = 1e-8;
[f1,~,~,~] = deriv(t,[X(1)+h,X(2),X(3),X(4)]);
[f2,~,~,~] = deriv(t,[X(1)-h,X(2),X(3),X(4)]);

dxdt = (f1-f2)/(2*h);


end

function dxdt = derivt2(t,X)

h = 1e-8;
[~,f1,~,~] = deriv(t,[X(1),X(2)+h,X(3),X(4)]);
[~,f2,~,~] = deriv(t,[X(1),X(2)-h,X(3),X(4)]);

dxdt = (f1-f2)/(2*h);


end


% Function for ode45
function dxdt = derivative(t,X)
m1 = 1;
m2 = 0.6;
L1 = 1;
L2 = 0.5;
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


% Hamiltonian
function ham = H(t,X)

m1 = 1;
m2 = 0.6;
L1 = 1;
L2 = 0.5;
g = 9.8;

T = (((X(3))^2*L2^2*m2+(X(4))^2*(m1+m2)*L1^2-m2*L1*L2*X(3)*X(4)*cos(X(1)-X(2))))/...
    (2*L1^2*L2^2*(m1+m2*(sin(X(1)-X(2)))^2));

U = -g*L1*(m1+m2)*cos(X(1))-g*m2*L2*cos(X(2));


ham = T+U;



end





