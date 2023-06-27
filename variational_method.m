close all;
clear; %clc;
mu = 1;
lam = 1;
N=100;

T = 10;
s = linspace(0,2*pi,N+1);
s = s(1:N);

dN = T/N;

%Cool ones [0,0,1,1] T = 24; 


%Approx periodic orbit location:
%works->

E_0 = -2.12;
q1 = 0;
q2 = acos((-E_0-(mu+1)*cos(q1))*lam);

x0 = [0;q2;-0.5;0];


%x0 = [0,0,1.39853100185856,0.0750000000000000];
[time,xy] = ode45(@derivative, 0:dN:T, x0);

figure(1)
plot(xy(:,1),xy(:,3),'-r.')
xlabel('\theta_1')
ylabel('p_1')
title([['x_0 = [',num2str((x0(1)),3)],[',',num2str(x0(2))],[',',num2str(x0(3),3)],[',',num2str(x0(4),3)],[']'],[',T_{initial} = ',num2str(T)]])

%delete last entry
%xy = xy(1:N,:);


% store loop points in 2Nx1 vector
% XY is 2*Nx1 (loop as column vector for matrix product derivative)
% XY = [x1;y1;x2;y2;...xN;yN]
XY=zeros(4*N,1);
for k=1:N
    XY(4*k-3:4*k)=xy(k,:)';
end

% construct finite difference matrix, d/ds
% this is for approximating loop tangent vector (v-twiddle)
D1=diag(8*ones(4*N-4,1),4);
D2=diag(-8*ones(4*N-4,1),-4);
D3=diag(-1*ones(4*N-8,1),8);
D4=diag(ones(4*N-8,1),-8);
D=D1+D2+D3+D4;

% handle periodic boundary conditions (top right/bottom left)
M1=eye(8)+diag(-8*ones(4,1),4);
M2=-1*M1';
D(1:8,4*N-7:4*N)=M1;
D(4*N-7:4*N,1:8)=M2;
% h is delta_s = delta_t (because i know T apriori)
% I don't think delta_t/delta_s will be the same in your case
h=2*pi/N;

D=(1/(12*h))*D;

% define scale factor
% YOU WILL NEED TO THINK ABOUT HOW TO DEAL WITH LAMBDA!
lam = dN/(s(2)-s(1));

% perturbation for forward diff approximation of dv/dx (shear matrix, A)
del=1e-7;


% set initial Euler step to something smaller than 1
% increase Euler step as n_RHS gets small
n_RHS=1; dtau=.2; number = 0;
while(n_RHS>1e-13)
    number = number+1;
    vel=zeros(4*N,1); t=0;
    
    %The diagonal stability matrix
    A = zeros(4*N,4*N);
    
    x2 = XY(5:8,1);
    for i = 1:N
        
        DD = stab(XY((4*i-3):4*i));
        
        A(4*i-3:4*i,4*i-3:4*i) = DD;
        
        %(actual)Velocity vector
        vel(4*i-3:4*i) = derivative(time,XY(4*i-3:4*i));
        %gradient of the Hamiltonian
        X2(4*i-3:4*i) = x2;
        gH(4*i-3:4*i) = gradham(XY(4*i-3:4*i)); 
        
    end
   
    
    
    % setup discretized equation for Newton step
    % I AM MYSTIFIED BY THIS...NOT HOW Ahat IS DEFINED IN Eq 18!!
    Ahat=D-lam*A;
    
    % v_L is loop tangent: v-twiddle
    v_L=D*XY;
    
    % increase Euler step as norm of RHS gets small (don't increase past 1)
    if n_RHS < 0.75
        dtau=min(1,2*dtau);
    end
    
    RHS=dtau*[lam*vel-v_L;0];
    
    n_RHS=norm(lam*vel-v_L);
    
    M = [Ahat, -vel;
        gH,0];
    
    dN= M\RHS;
    
    XY=XY+dN(1:4*N);
    lam = lam+dN(end);
    
    %eignls = sort(abs(eig(M)),'descend')
    xy=zeros(N,4);
    for m=1:N
        xy(m,1:4)=XY(4*m-3:4*m,1)';
    end
    %fft XY (fft shift, to order)
    
    %Fs = N/(time(2)-time(1));
%     
%     for l = 1:4
%     fXY = fft(xy(:,l));
%     f = Fs*(1:N)/N;
%     fXY = fftshift(fXY);
%     fXY(end) = 0;
%     xy(:,l) = ifft(fXY);
%     end
    
%     for q=1:N
%         XY(4*q-3:4*q,1) = xy(k,:)';
%     end
    % inspect how loop changes on each Newton step
    % construct Nx2 phase space loop from column vector XY
%     xy=zeros(N,2);
%     
%     for m=1:N
%         xy(m,1:4)=XY(4*m-3:4*m,1)';
%     end
%     
    figure(2); plot(xy(:,1),xy(:,3),'g.-');
    xlabel('\theta_1')
    ylabel('p_1')
    
%     xlabel('X'); ylabel('Y');
for i =1:N-4
E(i) = H(XY(4*i-3:4*i));
end
end
hold on
E_err = abs(max(E)-min(E))/min(E);
% construct Nx2 phase space loop from column vector XY
% xy=zeros(N,2);
% for m=1:N
%     xy(m,1:2)=XY(2*m-1:2*m,1)';
% end
% 
% figure(2); plot(xy(:,1),xy(:,2),'g.-');
% xlabel('X'); ylabel('Y');
% hold on
% % plot starting pt and ending pt
% plot(xy(1,1),xy(1,2),'rx',xy(end,1),xy(end,2),'b*')
% 
% figure(3); plot3(xy(:,1),xy(:,2),tsto,'g.-'); grid on
% xlabel('X'); ylabel('Y'); zlabel('Time')
T = 2*pi*lam;
[t,perorb] = ode45(@derivative,[0,T],XY(1:4));
title(['T = ',num2str(T)])
figure(3)
plot(perorb(:,1),perorb(:,3))
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

function A = stab(X)
mu = 1;
lam = 1;
A = zeros(4,4);

A(1,1) = lam.^(-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*(lam.*(1+mu).* ...
    X(4)+(-2).*X(3).*cos(X(1)+(-1).*X(2))+lam.*X(4).*cos(X(1)+(-1).*X(2)).^2).*sin( ...
    X(1)+(-1).*X(2));
A(2,1) = (1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*((1+mu).*X(3)+(-2).*lam.*(1+ ...
    mu).*X(4).*cos(X(1)+(-1).*X(2))+X(3).*cos(X(1)+(-1).*X(2)).^2).*sin(X(1)+(-1).* ...
    X(2));
A(3,1) = lam.^(-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-3).*((-4).*cos(X(1)+( ...
    -1).*X(2)).*sin(X(1)+(-1).*X(2)).*((-1/4).*lam.*(1+mu).*((-1)+(-2).*mu+ ...
    cos(2.*(X(1)+(-1).*X(2)))).^2.*sin(X(1))+((-1).*lam.*(1+mu).*X(4)+X(3).*cos( ...
    X(1)+(-1).*X(2))).*(X(3)+(-1).*lam.*X(4).*cos(X(1)+(-1).*X(2))).*sin(X(1)+(-1).* ...
    X(2)))+(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).*(cos(X(1)+(-1).*X(2)).*((-1).* ...
    lam.*(1+mu).*X(4)+X(3).*cos(X(1)+(-1).*X(2))).*(X(3)+(-1).*lam.*X(4).*cos(X(1)+( ...
    -1).*X(2)))+(-1/4).*lam.*(1+mu).*cos(X(1)).*((-1)+(-2).*mu+cos(2.*(X(1)+ ...
    (-1).*X(2)))).^2+lam.*X(4).*((-1).*lam.*(1+mu).*X(4)+X(3).*cos(X(1)+(-1).* ...
    X(2))).*sin(X(1)+(-1).*X(2)).^2+X(3).*((-1).*X(3)+lam.*X(4).*cos(X(1)+(-1).*X(2))) ...
    .*sin(X(1)+(-1).*X(2)).^2+lam.*(1+mu).*((-1)+(-2).*mu+cos(2.*(X(1)+(-1) ...
    .*X(2)))).*sin(X(1)).*sin(2.*(X(1)+(-1).*X(2)))));

A(4,1) = lam.^(-1).*(1+2.*mu+(-1).*cos(2.*(X(1)+(-1).*X(2)))).^(-3).*(6.*X(3).^2+ ...
    6.*lam.^2.*X(4).^2+6.*lam.^2.*mu.*X(4).^2+lam.*((-11)+4.*mu+8.*mu.^2) ...
    .*X(3).*X(4).*cos(X(1)+(-1).*X(2))+(-4).*(1+2.*mu).*(X(3).^2+lam.^2.*(1+mu) ...
    .*X(4).^2).*cos(2.*(X(1)+(-1).*X(2)))+(21/2).*lam.*X(3).*X(4).*cos(3.*(X(1)+( ...
    -1).*X(2)))+12.*lam.*mu.*X(3).*X(4).*cos(3.*(X(1)+(-1).*X(2)))+(-2).*X(3).^2.* ...
    cos(4.*(X(1)+(-1).*X(2)))+(-2).*lam.^2.*X(4).^2.*cos(4.*(X(1)+(-1).*X(2)))+( ...
    -2).*lam.^2.*mu.*X(4).^2.*cos(4.*(X(1)+(-1).*X(2)))+(1/2).*lam.*X(3).*X(4).* ...
    cos(5.*(X(1)+(-1).*X(2))));

A(1,2) = (-1).*lam.^(-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*(lam.*(1+ ...
    mu).*X(4)+(-2).*X(3).*cos(X(1)+(-1).*X(2))+lam.*X(4).*cos(X(1)+(-1).*X(2)).^2).* ...
    sin(X(1)+(-1).*X(2));

A(2,2) = (-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*((1+mu).*X(3)+(-2).* ...
    lam.*(1+mu).*X(4).*cos(X(1)+(-1).*X(2))+X(3).*cos(X(1)+(-1).*X(2)).^2).*sin( ...
    X(1)+(-1).*X(2));

A(3,2) = lam.^(-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-3).*(4.*cos(X(1)+(-1) ...
    .*X(2)).*sin(X(1)+(-1).*X(2)).*((-1/4).*lam.*(1+mu).*((-1)+(-2).*mu+cos( ...
    2.*(X(1)+(-1).*X(2)))).^2.*sin(X(1))+((-1).*lam.*(1+mu).*X(4)+X(3).*cos(X(1)+( ...
    -1).*X(2))).*(X(3)+(-1).*lam.*X(4).*cos(X(1)+(-1).*X(2))).*sin(X(1)+(-1).*X(2))) ...
    +(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).*((-1).*(X(3).^2+lam.^2.*(1+mu).* ...
    X(4).^2).*cos(X(1)+(-1).*X(2)).^2+lam.*X(3).*X(4).*cos(X(1)+(-1).*X(2)).^3+lam.* ...
    X(3).*X(4).*cos(X(1)+(-1).*X(2)).*(mu+cos(2.*(X(1)+(-1).*X(2))))+(X(3).^2+ ...
    lam.^2.*(1+mu).*X(4).^2).*sin(X(1)+(-1).*X(2)).^2+lam.*(1+mu).*(1+2.*mu+ ...
    (-1).*cos(2.*(X(1)+(-1).*X(2)))).*sin(X(1)).*sin(2.*(X(1)+(-1).*X(2)))));

A(4,2) = lam.^(-1).*((-1)+(-2).*mu+cos(2.*(X(1)+(-1).*X(2)))).^(-3).*((-4).* ...
    sin(2.*(X(1)+(-1).*X(2))).*(2.*lam.*X(3).*X(4).*(3+2.*mu+cos(2.*(X(1)+(-1).* ...
    X(2)))).*sin(X(1)+(-1).*X(2))+(-2).*(X(3).^2+lam.^2.*(1+mu).*X(4).^2).*sin( ...
    2.*(X(1)+(-1).*X(2)))+(-1).*((-1)+(-2).*mu+cos(2.*(X(1)+(-1).*X(2)))).^2.* ...
    sin(X(2)))+((-1)+(-2).*mu+cos(2.*(X(1)+(-1).*X(2)))).*(4.*(X(3).^2+ ...
    lam.^2.*(1+mu).*X(4).^2).*cos(2.*(X(1)+(-1).*X(2)))+(-2).*lam.*X(3).*X(4).* ...
    cos(X(1)+(-1).*X(2)).*(3+2.*mu+cos(2.*(X(1)+(-1).*X(2))))+(-1).*((-1)+(-2) ...
    .*mu+cos(2.*(X(1)+(-1).*X(2)))).^2.*cos(X(2))+8.*lam.*X(3).*X(4).*cos(X(1)+( ...
    -1).*X(2)).*sin(X(1)+(-1).*X(2)).^2+(-4).*((-1)+(-2).*mu+cos(2.*(X(1)+(-1) ...
    .*X(2)))).*sin(2.*(X(1)+(-1).*X(2))).*sin(X(2))));

A(1,3) = (lam+lam.*mu+(-1).*lam.*cos(X(1)+(-1).*X(2)).^2).^(-1);

A(2,3) = cos(X(1)+(-1).*X(2)).*((-1)+(-1).*mu+cos(X(1)+(-1).*X(2)).^2).^(-1);

A(3,3) = (-1).*lam.^(-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*(lam.*(1+ ...
    mu).*X(4)+(-2).*X(3).*cos(X(1)+(-1).*X(2))+lam.*X(4).*cos(X(1)+(-1).*X(2)).^2).* ...
    sin(X(1)+(-1).*X(2));

A(4,3) = 2.*lam.^(-1).*((-1)+(-2).*mu+cos(2.*(X(1)+(-1).*X(2)))).^(-2).*((-4).* ...
    X(3).*cos(X(1)+(-1).*X(2))+lam.*X(4).*(3+2.*mu+cos(2.*(X(1)+(-1).*X(2))))).* ...
    sin(X(1)+(-1).*X(2));

A(1,4) = cos(X(1)+(-1).*X(2)).*((-1)+(-1).*mu+cos(X(1)+(-1).*X(2)).^2).^(-1);
A(2,4) = lam.*(1+mu).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-1);

A(3,4) = (-1).*(1+mu+(-1).*cos(X(1)+(-1).*X(2)).^2).^(-2).*((1+mu).*X(3)+(-2).* ...
    lam.*(1+mu).*X(4).*cos(X(1)+(-1).*X(2))+X(3).*cos(X(1)+(-1).*X(2)).^2).*sin( ...
    X(1)+(-1).*X(2));

A(4,4) = 2.*((-1)+(-2).*mu+cos(2.*(X(1)+(-1).*X(2)))).^(-2).*((-4).*lam.*(1+mu) ...
    .*X(4).*cos(X(1)+(-1).*X(2))+X(3).*(3+2.*mu+cos(2.*(X(1)+(-1).*X(2))))).*sin( ...
    X(1)+(-1).*X(2));

end
function dxdt = gradham(X)
mu = 1;
lam = 1;

dxdt1 = (X(3)*1/lam-X(4)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);
dxdt2 = (lam*(mu+1)*X(4)-X(3)*cos(X(1)-X(2)))/(mu+(sin(X(1)-X(2)))^2);
dxdt3 = -(X(3)*X(4)*sin(X(1)-X(2)))/((mu+(sin(X(1)-X(2)))^2))+...
    (sin(2*(X(1)-X(2)))*(1/2*1/lam*(X(3))^2+1/2*lam*(mu+1)*(X(4))^2-X(3)*X(4)*cos(X(1)-X(2))))/...
    (mu+(sin(X(1)-X(2)))^2)^2-(mu+1)*sin(X(1));
dxdt4 = (X(3)*X(4)*sin(X(1)-X(2)))/((mu+(sin(X(1)-X(2)))^2))...
    -(sin(2*(X(1)-X(2)))*(1/2*1/lam*(X(3))^2+1/2*lam*(mu+1)*...
    (X(4))^2-X(3)*X(4)*cos(X(1)-X(2))))/(mu+(sin(X(1)-X(2)))^2)^2-1/lam*sin(X(2));

dxdt = [-dxdt3;-dxdt4;dxdt1;dxdt2];

end
function ham = H(X)
mu = 1;
lam = 1;

ham = (-1).*(1+mu).*cos(X(1))+lam.^(-1).*((X(3).^2+lam.^2.*(1+mu).*X(4).^2+( ...
  -2).*lam.*X(3).*X(4).*cos(X(1)+(-1).*X(2))).*(1+2.*mu+(-1).*cos(2.*(X(1)+( ...
  -1).*X(2)))).^(-1)+(-1).*cos(X(2)));



end