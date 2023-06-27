error = 100;
n = 2;
x = 1/2*(log(13)-log(5));
f = @(x) x/(x^2+4);
a = 1;
b = 3;
while error>=0.1
   
    value = numint(f,n,a,b,'mid');
    error = abs((x-value)/x)*100;
    n = n+2;
    
end
n = n-2
m = n/2+1