
m = 4*1.66*10^(-27);
T = 1000;
k = 1.38*10^(-23);

P =@(v) exp(-m*v.^2/(2*k*T))*(2*pi*m*k*T)^(-3/2)*m^3*4*pi.*v.^2;

probability = integral(P,2.38*10^3,inf)