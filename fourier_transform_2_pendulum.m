
m1 = 10;
m2 = 1;
L1 = 1;
L2 = L1;
g = 9.8;

Fs = 1/10;                  % Sampling frequency                    
T = 1/Fs;                   % Sampling period       
L = 1000;                   % Length of signal
t = 0:0.01:(L-1)*T;         % Time vector


omega1 = sqrt(g/(m1*L1)*(m1+m2+sqrt((m1+m2)*m2)));
omega2= sqrt(g/(m1*L1)*(m1+m2-sqrt((m1+m2)*m2)));


x0 = [pi/180;0;0;0];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);


[T,y] = ode45(@derivative,[0,t(end)],x0,options);

Y = fft(y);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

omega1_fft = max(P1)*(2*pi);
P1new = setdiff(P1,[max(P1)]);
omega2_fft = max(P1new)*(2*pi);

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