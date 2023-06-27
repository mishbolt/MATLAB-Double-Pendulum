clear all;
clc;

V=13;
f=2/5;
R=3;
L=6;
C=1/5;
%skip=250;

% Initial Conditions
tf=25;
dt=.05;
N=tf/dt;

% Arrays
myt=0:dt:tf;
myx=zeros(N,2);
myx(1,:)=[0,8];

% Heun's Method
for n=1:N
    xtmp=myx(n,:)+dt*deriv(myt(n),myx(n,:))';
    myx(n+1,:)=myx(n,:)+.5*(deriv(myt(n),myx(n,:))'+deriv(myt(n+1),xtmp)')*dt;
end

% ODE45 Solution
[t,x]=ode45(@deriv,[0,tf],myx(1,:));

% Plot
figure(1);
%plot(myt(1:skip:end),myx(1:skip:end,1),'k-o')
plot(myt,myx(:,1),'-rd',t,x(:,1),'-ko')
grid on;
legend('ODE45','Heun''s method');
xlabel('time in seconds');
ylabel('current in amps');
title('Current vs. Time');

figure(2);
plot(myt,myx(:,2),'-rd',t,x(:,2),'-ko');
grid on;
legend('ODE45','Heun''s method');
xlabel('time in seconds');
ylabel('potential difference in volts');
title('Potential Difference Acorss Inductor vs, Time');

% Subfunction
function dxdt=deriv(t,x)
V=13;
f=2/5;
R=3;
L=6;
C=1/5;
w=2*pi*f;
dxdt=[x(2);-R/L*x(2)-1/(L*C)*x(1)-1/L*V*w*sin(w*t)]; 
end