
close all; 
clear; clc;

Nx=128; Lx=4*pi;
Ny=64; Ly=2*pi;

T=365.80;

x=zeros(Nx,1); y=zeros(Ny,1);
psi=zeros(Ny,Nx);
for k=1:Ny
    for l=1:Nx
        x(l)=(Lx*l)/Nx;
        y(k)=(Ly*k)/Ny;
        psi(k,l)=psip1(0,[x(l),y(k)]);      
    end
end

figure(1);
contour(x,y,psi)
xlabel('X'); ylabel('Y');
title('Stream function contours with initial guess of center PO')
hold on

% locations of the 8 saddles
sad=[pi/2 pi/4; pi/2 5*pi/4; 3*pi/2 3*pi/4; 3*pi/2 7*pi/4;...
    5*pi/2 pi/4; 5*pi/2 5*pi/4; 7*pi/2 3*pi/4; 7*pi/2 7*pi/4];

N=512; 
nsad=5;
xy=ones(N,2);

% % initial guess for saddle
% xy(:,1)=xy(:,1)*sad(nsad,1)+.1;
% xy(:,2)=xy(:,2)*sad(nsad,2)+.1;

center=[3*pi/2; 5*pi/4];
xy(:,1)=xy(:,1)*center(1);
xy(:,2)=xy(:,2)*center(2);

figure(1);
plot(center(1),center(2),'rx')

XY=zeros(2*N,1);
for k=1:N
    XY(2*k-1:2*k)=xy(k,:)';
end

D1=diag(8*ones(2*N-2,1),2);
D2=diag(-8*ones(2*N-2,1),-2);
D3=diag(-1*ones(2*N-4,1),4);
D4=diag(ones(2*N-4,1),-4);
D=D1+D2+D3+D4;

M1=eye(4)+diag(-8*ones(2,1),2);
M2=-1*M1';
D(1:4,2*N-3:2*N)=M1;
D(2*N-3:2*N,1:4)=M2;

h=T/N; D=(1/(12*h))*D;

lam=T/N; 
del=1e-7;

t=0; tsto=zeros(N,1);
for k=1:N
    t=t+h;
    tsto(k,1)=t;
end

n_RHS=1; dtau=1; 
while(n_RHS>1e-10)
    vel=zeros(2*N,1); A=zeros(2*N); t=0;
    for k=1:N
        t=t+h;
        
        r=[XY(2*k-1); XY(2*k)];
        v=velp1(t,r);
        vel(2*k-1:2*k)=[v(1); v(2)];
       
        rdx=[XY(2*k-1)+del; XY(2*k)];
        vdx=velp1(t,rdx);
        
        rdy=[XY(2*k-1); XY(2*k)+del];
        vdy=velp1(t,rdy);
        
        DG=zeros(2);
        
        DG(1,1)=vdx(1) - v(1);
        DG(1,2)=vdy(1) - v(1);
        DG(2,1)=vdx(2) - v(2);
        DG(2,2)=vdy(2) - v(2);
        
        DG=1/del*DG;
        
        A(2*k-1:2*k,2*k-1:2*k)=DG;
    end

    Ahat=D-lam*A;
    
    v_L=D*XY;
    RHS=dtau*[lam*vel-v_L];
    
    n_RHS=norm(lam*vel-v_L);
    
    dN=Ahat\RHS;
    XY=XY+dN;
end

xy=zeros(N,2);
for m=1:N
    xy(m,1:2)=XY(2*m-1:2*m,1)';
end

figure(2); plot(xy(:,1),xy(:,2),'b.-');

figure(3); plot3(xy(:,1),xy(:,2),tsto,'b.-'); grid on



