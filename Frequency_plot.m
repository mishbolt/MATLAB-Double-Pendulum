mu = linspace(0.1,2,1000);
lam = 1;

omega1 = sqrt(1/2*(1+1./mu+1./lam+1./mu.*1./lam)+1/2*sqrt((1+1./mu+1./lam+1./mu.*1./lam).^2-4*1./lam.*(1+1./mu)));


omega2 = sqrt(1/2*(1+1./mu+1./lam+1./mu.*1./lam)-1/2*sqrt((1+1./mu+1./lam+1./mu.*1./lam).^2-4*1./lam.*(1+1./mu)));

figure
plot(mu,omega1,mu,omega2)

%figure
%plot(omega1,lam,omega2,lam)