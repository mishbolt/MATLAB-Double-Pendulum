
[t,x] = ode45(@f,[0,10],5);

function func = f(t,X)

func = 5*X;

end