a_max = sqrt(b^2/(4*Pdry^2)+9*etaT*d./((t_fallAVG-t_fallSD)*2*rhop*g))-b/(2*Pdry);
a_min = sqrt(b^2/(4*Pdry^2)+9*etaT*d./((t_fallAVG+t_fallSD)*2*rhop*g))-b/(2*Pdry);
a_avg = (a_min+a_max)/2;
a_SD = a_max-a_avg;

q_max= 6*pi*D/V*(9*etaT^3/(2*rhop*g))^(1/2).*1./(1+b./(a_max*Pdry)).^(3/2).*...
    (d./(t_fallAVG-t_fallSD)).^(1/2).*(d./(t_fallAVG-t_fallSD)+d./(t_riseAVG-t_riseSD));
q_min= 6*pi*D/V*(9*etaT^3/(2*rhop*g))^(1/2).*1./(1+b./(a_min*Pdry)).^(3/2).*...
    (d./(t_fallAVG+t_fallSD)).^(1/2).*(d./(t_fallAVG+t_fallSD)+d./(t_riseAVG+t_riseSD));
q_AVG = (q_min+q_max)/2;
q_total = mean(q_AVG);

q_SD = q_max-q_AVG;

[q_AVG_asc,index] = sort(q_AVG);



integer_e =(q_max./min(q_max));
int_e_r = round(integer_e);

X = sum(q_max.^2)/(sum(int_e_r.*q_max));

integer_e_new = q_max/X;
charge_max = max(integer_e_new-int_e_r)*X;

integer_e =(q_min./min(q_min));
int_e_r = round(integer_e);

X = sum(q_min.^2)/(sum(int_e_r.*q_min));

integer_e_new = q_min/X;
charge_min = max(integer_e_new-int_e_r)*X;

charge_AVG = (charge_max+charge_min)/2;

charge_SD = charge_max-charge_AVG;




for k=1:15
err(k) = q_SD(index(k));
end

figure()
errorbar(1:15,q_AVG_asc,err,'ro')
grid on
title('Data points for the charge,q, with uncertainties')
ylabel('q (C)')
xlabel('number of the drop')
hold on
plot(linspace(1,15),e*((linspace(1,15))-1)+q_AVG_asc(1))



