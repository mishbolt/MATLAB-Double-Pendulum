%-------------------------------------------------------------------------
%Programer:  Grant Fuchs
%
%Purpose: This script numerically solves a differential equation with 
%a Forward Euler Method and a MATLAB built in function. It also generates
%a plot of the solutions.
%
%Variables:
%dt, tf    -timestep(s),end time(s)
%x0        -initial condidtions
%t,x       -solutions from ode45
%g,R      -gravitational acceleration(m/s^2),Radius of the Earth(m)
%-------------------------------------------------------------------------



dt = 0.002;
tf = 20000;
 
x0 = [-6370 1];
X(1,:) = x0;
%variable initialization




[t, x] = ode45(@spring,[0 tf],x0);
%ode 45 fuction call

figure(1); 
plot(t,x(:,1),'kx-')
xlim([0 20000]);
xlabel('Time (s)');
grid on 
ylabel('Displacement (rad)'); 
legend('ode45','location', 'best') 
title('Displacement as a Function of Time')
%plotting for fig 1

figure(2);
plot(x(:,1),x(:,2),'kx-')
xlim([0 20000]); 
xlabel('Time (s)'); 
grid on 
ylabel('Velocity (m/s)'); 
legend('ode45','location', 'best') 
title('Velocity as a Function of Time')
%plotting fig 2








%THIS WHOLE THING IS THE DENSITY MODEL
function mass = integr(r)
dens1 = @(x) (4*pi.*x.^2).*(13.0885-8.8381*(x./6371).^2)*10^12;
dens2 = @(x) (4*pi.*x.^2).*(12.5815-1.2638.*(x./6371)-3.6426.*(x./6371).^2-5.5281.*(x./6371).^3)*10^12;
dens3 = @(x) (4*pi.*x.^2).*(7.9565-6.4761.*(x./6371)+5.5283.*(x./6371).^2-3.0807.*(x./6371).^3)*10^12;
dens4 = @(x) (4*pi.*x.^2).*(5.3197-1.4836.*(x./6371))*10^12;
dens5 = @(x) (4*pi.*x.^2).*(11.2494-8.0298.*(x./6371))*10^12;
dens6 = @(x) (4*pi.*x.^2).*(7.1089 -3.8045.*(x./6371))*10^12;
dens7 = @(x) (4*pi.*x.^2).*(2.6910 +0.6924.*(x./6371))*10^12;
dens8 =  2.900*10^12;
dens9 =  2.600*10^12;
dens10 =  1.020*10^12;
mass1max = integral(dens1,0,1221.5);
mass2max = mass1max + integral(dens2,1221.5,3480);
mass3max = mass2max + integral(dens3,3480,5701);
mass4max = mass3max + integral(dens4,5701,5771);
mass5max = mass4max + integral(dens5,5771,5971);
mass6max = mass5max + integral(dens6,5971,6151);
mass7max = mass6max + integral(dens7,6151,6346.6);
mass8max = mass7max + dens8*(4/3*pi*6356^3-4/3*pi*6346.6^3);
mass9max = mass8max + dens9*(4/3*pi*6368^3-4/3*pi*6356^3);
maxmass = 5.972*10^24;

if (0 <= abs(r)) && (abs(r) <= 1221.5)
    dens = dens1;
    mass = integral(dens,0,abs(r));
    
end
if (1221.5 < abs(r)) && (abs(r) <=3480)
    dens = dens2;
    mass2 = integral(dens,1221.5,abs(r));
    mass = mass1max + mass2;
end 
if (3480 < abs(r)) && (abs(r) <= 5701)
    
    dens = dens3;
    mass3 = integral(dens,3480,abs(r));
    mass = mass2max + mass3;
end
if (5701 < abs(r)) && (abs(r) <=5771)
    
    dens = dens4;
    mass4 = integral(dens,5701,abs(r));
    mass = mass3max + mass4;
end 
if (5771 < abs(r)) && (abs(r) <= 5971)
    
    dens = dens5;
    mass5 = integral(dens,5771,abs(r));
    mass = mass4max + mass5;
    
end 
 if (5971 < abs(r)) && (abs(r) <=6151)
     
    dens = dens6;
    mass6 = integral(dens,5971,abs(r));
    mass = mass5max+ mass6;
    
end 
if (6151 < abs(r)) && (abs(r) <=6346.6)
    
    dens = dens7;
    mass7 = integral(dens,6151,abs(r));
    mass = mass6max+ mass7;
end 
if (6346.6 < abs(r)) && (abs(r) <= 6356)
    
    dens = dens8;
    mass8 = dens*(4*1/3*pi*abs(r)^3-4*1/3*pi*6346^3);
    mass = mass7max+ mass8;
end 
if (6356 < abs(r)) && (abs(r) <=6386)
    
    dens = dens9;
    mass9 = dens*(4*1/3*pi*abs(r)^3-4*1/3*pi*6356^3);
    mass = mass8max+ mass9;
end 
if (6386 < abs(r)) && (abs(r) <= 6371)
    
    dens = dens10;
    mass10 = dens*(4*1/3*pi*abs(r)^3-4*1/3*pi*6386^3);
    mass = mass9max+ mass10;
end 
  if abs(r) > 6371
      mass = maxmass;
  end
end
%THIS WHOLE THING IS THE DENSITY MODEL

function dxdt = spring(t,x) 
    G = 6.674*10^-20;
    dxdt = [x(2); -(G*integr(x(1)))*abs(x(1))/(x(1)^3)]; 
end
%just the derivative subfunction