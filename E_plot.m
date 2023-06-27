% mu = 1;
% lam = 1;
% q1 = linspace(-pi,pi);
% q2 = linspace(-pi,pi);
% 
% [qq1,qq2] = meshgrid(q1,q2);
% 
% H  = (-1).*(1+mu).*cos(qq1)+lam.^(-1).*((0.^2+lam.^2.*(1+mu).*0.^2+( ...
%   -2).*lam.*0.*0.*cos(qq1+(-1).*qq2)).*(1+2.*mu+(-1).*cos(2.*(qq1+( ...
%   -1).*qq2))).^(-1)+(-1).*cos(qq2));
% 
% figure
% surfc(q1,q2,H,'LineStyle','none')
% grid on
% figure
% contour(q1,q2,H)

% lam = linspace(1/4,4);
% W = linspace(1/4,3);
% 
% [lam1,W1] = meshgrid(lam,W);
% 
% mu = (1+lam1).^2.*W1./(lam1-W1-lam1.^2.*W1+lam1.*W1.^2);
% 
% surfc(lam,W,mu)


lam1 = 1;
W1 = 3/5;

% [lam1,W1] = meshgrid(lam,W);

mu = (1+lam1).^2.*W1.^2./(lam1-W1.^2-lam1.^2.*W1.^2+lam1.*W1.^4);
%  plot(lam1,mu)
% surfc(lam,W,mu)

% mu1 = linspace(1/3,3,100);
% lam1 = linspace(1/3,3,100);

% [mu,lam] = meshgrid(mu1,lam1);
% w1 = real(((mu.*lam+lam+mu+1)+sqrt((mu.*lam+lam+mu+1).^2-4.*mu.^2.*lam.*(1+mu.^(-1))))./(2.*mu.*lam));
% w2 = real(((mu.*lam+lam+mu+1)-sqrt((mu.*lam+lam+mu+1).^2-4.*mu.^2.*lam.*(1+mu.^(-1))))./(2.*mu.*lam));
% W = real(sqrt(w1/w2));
% 
% contour(mu1,lam1,W)






