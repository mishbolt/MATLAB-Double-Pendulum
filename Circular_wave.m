syms t
A = 2;
k = 2*pi;
z0 = 0;
w = 1;

f_z = A*cos(k*z0-w*t);
f_y = -A*sin(k*z0-w*t);
x = t;

ax = gca;
ax.XAxisLocation = 'origin';
fplot3(x,f_y,f_z,[0,10])
