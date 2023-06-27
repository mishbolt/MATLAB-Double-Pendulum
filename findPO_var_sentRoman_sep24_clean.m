close all;
clear; %clc;


T=365.80;
N=512;
xy=ones(N,2);

% approximate location of centers
cen = [pi/2 3*pi/4; pi/2 7*pi/4; 3*pi/2 pi/4; 3*pi/2 5*pi/4;...
    5*pi/2 3*pi/4; 5*pi/2 7*pi/4; 7*pi/2 pi/4; 7*pi/2 5*pi/4];

% change ncen to pick out different centers for computing PO

% NOTE: my initial guess for a loop is not really a loop at all!
% it is something close to the location of the center (see psi plot
% in the other version of this code) AND this point remains fixed
% in place over the period of the flow

% Initial guess for loop
Ncen = 4;
xy(:,1)=xy(:,1)*cen(Ncen,1);
xy(:,2)=xy(:,2)*cen(Ncen,2);

% store loop points in 2Nx1 vector
% xy is Nx2 (phase space loop)
% XY is 2*Nx1 (loop as column vector for matrix product derivative)
% XY = [x1;y1;x2;y2;...xN;yN]
XY=zeros(2*N,1);
for k=1:N
    XY(2*k-1:2*k)=xy(k,:)';
end

% construct finite difference matrix, d/ds
% this is for approximating loop tangent vector (v-twiddle)
D1=diag(8*ones(4*N-4,1),4);
D2=diag(-8*ones(4*N-4,1),-4);
D3=diag(-1*ones(4*N-8,1),8);
D4=diag(ones(4*N-8,1),-8);
D=D1+D2+D3+D4;

% handle periodic boundary conditions (top right/bottom left)
M1=eye(8)+diag(-8*ones(4,1),5);
M2=-1*M1';
D(1:8,4*N-7:4*N)=M1;
D(4*N-7:2*N,1:8)=M2;

% h is delta_s = delta_t (because i know T apriori)
% I don't think delta_t/delta_s will be the same in your case
h=T/N;

D=(1/(12*h))*D;

% define scale factor
% YOU WILL NEED TO THINK ABOUT HOW TO DEAL WITH LAMBDA!
lam=T/N;

% perturbation for forward diff approximation of dv/dx (shear matrix, A)
del=1e-7;

% for plotting PO with 3rd axis as time (figure 3)
t=0; tsto=zeros(N,1);
for k=1:N
    t=t+h;
    tsto(k,1)=t;
end

% set initial Euler step to something smaller than 1
% increase Euler step as n_RHS gets small
n_RHS=1; dtau=.2;
while(n_RHS>1e-10)
    
    vel=zeros(2*N,1); A=zeros(2*N); t=0;
    
    for k=1:N
        t=t+h;
        
        % XY = [x1;y1;x2;y2;...xN;yN]
        % ref pt for forward diff approx to A = dv/dx
        r=[XY(2*k-1); XY(2*k)];
        v=velp1(t,r);
        vel(2*k-1:2*k)=[v(1); v(2)];
        % NOTE: vel is v_hat on RHS of Eq 18
        
        % perturb x-coord
        rdx=[XY(2*k-1)+del; XY(2*k)];
        vdx=velp1(t,rdx);
        
        % perturb y-coord
        rdy=[XY(2*k-1); XY(2*k)+del];
        vdy=velp1(t,rdy);
        
        DG=zeros(2);
        
        DG(1,1)=vdx(1) - v(1);
        DG(1,2)=vdy(1) - v(1);
        DG(2,1)=vdx(2) - v(2);
        DG(2,2)=vdy(2) - v(2);
        
        DG=1/del*DG;
        
        % ultimately this is diag[A1,A2,...,AN] in Eq 18
        % each iteration of inner loop creates A1, A2, ...
        A(2*k-1:2*k,2*k-1:2*k)=DG;
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
    
    RHS=dtau*[lam*vel-v_L];
    
    n_RHS=norm(lam*vel-v_L);
    
    dN=Ahat\RHS;
    XY=XY+dN;
    
    % inspect how loop changes on each Newton step
    % construct Nx2 phase space loop from column vector XY
    xy=zeros(N,2);
    for m=1:N
        xy(m,1:2)=XY(2*m-1:2*m,1)';
    end
    
    figure(10); plot(xy(:,1),xy(:,2),'g.-');
    xlabel('X'); ylabel('Y');
    pause
end

% construct Nx2 phase space loop from column vector XY
xy=zeros(N,2);
for m=1:N
    xy(m,1:2)=XY(2*m-1:2*m,1)';
end

figure(2); plot(xy(:,1),xy(:,2),'g.-');
xlabel('X'); ylabel('Y');
hold on
% plot starting pt and ending pt
plot(xy(1,1),xy(1,2),'rx',xy(end,1),xy(end,2),'b*')

figure(3); plot3(xy(:,1),xy(:,2),tsto,'g.-'); grid on
xlabel('X'); ylabel('Y'); zlabel('Time')




