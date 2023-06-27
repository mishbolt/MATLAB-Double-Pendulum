
mu = 1;
lam = 1;

q1 = pi/180; p1(1) = 0;
q2 = pi/90; p2(1) = 0;


 %solve for set of ICs with the same energy
N=40;
 
options = odeset('Events',@pendevents,'RelTol',1e-6,'AbsTol',1e-8);
% Set the loop for various initial condtions
step = pi/100;
%h = figure;
%M(N) = struct('cdata',[],'colormap',[]);
%h.Visible = 'off';
for i=1:N
% Set initial parameters

q1 = q1-step;
q2 = 0;
E0 = H([q1,q2,p1(1),p2(1)]);
q2 = real(acos((-E0-(mu+1)*cos(q1))*lam));
tf = 5000;


x0 = [q1;q2;p1(1);p2(1)];
[T,x,te,ye,ie] = ode45(@derivative,[0,tf],x0,options);

figure
plot(ye(:,2),ye(:,4),'r.')
ylim([-1,1])
xlim([-1,1])

xlabel('\theta_2'); ylabel('p_\theta_2')

title(['Poincare Section, $\theta_1 = 0$, $\dot{\theta_1}>0$, $E_0 \approx$', num2str(round(E0,2))],'Interpreter','latex')
text(0,mean(ye(:,4)),['$E_0 \approx$', num2str(round(E0,2))],'Interpreter','latex')


M(i) = getframe;

end
%h.Visible = 'on';
%movie(M,5);
%Boundary Curve

% boundary = @(x,y) E0-(-1-mu+1/2.*lam.*y.^2-cos(x)./lam);
% fimplicit(boundary,'k-')



% Function for ode45

%Algebraically computed
function dxdt = derivative(t,X)
mu = 1;
lam = 1;

dxdt1 = (X(3)*1/lam-X(4)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);

dxdt2 = (lam*(mu+1)*X(4)-X(3)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);

A = -(X(3)*X(4)*sin(X(1)-X(2)))/((mu+(sin(X(1)-X(2)))^2));
B = (sin(2*(X(1)-X(2)))*(1/2*1/lam*(X(3))^2+1/2*lam*(mu+1)*(X(4))^2-X(3)*X(4)*cos(X(1)-X(2))))/(mu+(sin(X(1)-X(2)))^2)^2;

dxdt3 = A+B-(mu+1)*sin(X(1));

dxdt4 = -A-B-1/lam*sin(X(2));

dxdt = [dxdt1;dxdt2;dxdt3;dxdt4];



end


function ham = H(X)
mu = 1;
lam = 1;

ham = (-1).*(1+mu).*cos(X(1))+lam.^(-1).*((X(3).^2+lam.^2.*(1+mu).*X(4).^2+( ...
  -2).*lam.*X(3).*X(4).*cos(X(1)+(-1).*X(2))).*(1+2.*mu+(-1).*cos(2.*(X(1)+( ...
  -1).*X(2)))).^(-1)+(-1).*cos(X(2)));



end



% function ham = H2(q1,q2)
% 
% mu = 1;
% lam = 1;
% p1 = 0;
% p2 = 0;
% 
% ham = (-1).*(1+mu).*cos(q1)+lam.^(-1).*((p1.^2+lam.^2.*(1+mu).*p2.^2+( ...
%   -2).*lam.*p1.*p2.*cos(q1+(-1).*q2)).*(1+2.*mu+(-1).*cos(2.*(q1+( ...
%   -1).*q2))).^(-1)+(-1).*cos(q2));
% 
% 
% 
% end

% Function for events
function [value,isterminal,direction]=pendevents(t,y)

value = y(1);
isterminal = 0;
direction = 1;


end