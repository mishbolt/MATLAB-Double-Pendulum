q1(1) = -pi/180; p1(1) = 0;
q2(1) = pi/180; p2(1) = 0;

E0 = H([q1(1),q2(1),p1(1),p2(1)]);

mu = 1;
lam = 1;
 %solve for set of ICs with the same energy
N=20;
 q1 = linspace(q1(1),q1(1)+2*abs(q1(1)),N);
% 
% 
 q2 = acos((-E0-(mu+1)*cos(q1))*lam);

%p1 = linspace(p1(1),0,N);
%p2 = (p1+sqrt(p1.^2-2*lam*(mu+1)*...
   %(1/2*lam^(-1)*p1.^2-mu^2-mu-mu*lam^(-1)-E0*mu)))/(lam*(mu+1));

i = 1;
a = -1;
while a<0

[T,Res]=lyapunov(4,@double_pend,@ode45,0,0.1,600,[q1(i),q2(i),p1(i),p2(i)],1);
a = max(Res(end,:));
i = i+1;
end

function ham = H(X)
mu = 1;
lam = 1;

ham = (-1).*(1+mu).*cos(X(1))+lam.^(-1).*((X(3).^2+lam.^2.*(1+mu).*X(4).^2+( ...
  -2).*lam.*X(3).*X(4).*cos(X(1)+(-1).*X(2))).*(1+2.*mu+(-1).*cos(2.*(X(1)+( ...
  -1).*X(2)))).^(-1)+(-1).*cos(X(2)));



end