function estimate = numint(f,n,a,b,method)
% Thanks to Mary Begay for providing this template.
% Input Notes:
% f: @(x) FUNCTION	=> for example @(x) x^2
% a: lower limit	=> number
% b: upper limit 	=> number
% n: Interval       => number
% method:           => string, 'trp' or 'mid' or 'sim'

% Checking method input is a string
rule = ischar(method);
if rule == 0
    disp('Error: Method input must be a string')
    return
end
    
    %% Trapezoid sum Theorem 4.5
if method == 'trp'
    
    h = (b-a)/n;
    
    % Set initial values
    sum = 0;
    yold = feval(f,a);
    
    %Set the loop according to Theorem 4.5
    for j = 0:1:n-1
        xj = a+j*h; %x_j
        xjp1 = xj+h; %x_j+1
        ynew = feval(f,xjp1); % y value at x_j+1
        sum = sum+1/2*(yold+ynew); %area of a trapezoid
        yold = ynew;
    end
    % Output the result
    estimate = h*sum;
    
    %%
    %% Midpoint sum Theorem 4.6
elseif method == 'mid'
    
    % Set the condition that n has to be even
    if mod(n,2)~=0
        fprintf('Error: Number of intervals n has to be even')
        return
    end
    
    % Calculate h
    h = (b-a)/(n+2);
    
    sum = 0;
    
    % Set the loop according to theorem 4.6
    for j = 0:1:n/2
        xj = a+h*(2*j+1); % find x_2j
       
        sum = sum+f(xj); % calculate partial sum
    end
    estimate = 2*h*sum; % calculate the area=estimate
    
    %% Simpson's rule, Algorithm 4.1
    %%
elseif method == 'sim'
    
     % Set the condition that n has to be even
    if mod(n,2)~=0
        fprintf('Error: Number of intervals n has to be even')
        return
    end
    
    % Calculate h
    h = (b-a)/n;
    
    % Set initial conditions
    XI0 = f(a)+f(b);
    XI1 = 0;
    XI2 = 0;
    
    % Set the loop according to Algorithm 4.1
    for i = 1:n-1
        X = a+i*h; %calculate jth x
        
        if mod(i,2) == 0 % if i is even increment x2
            XI2 = XI2+f(X);
        else
            XI1 = XI1+f(X); % else increment x1
        end
        
        % Calculate final area under the parabola
        estimate = h/3*(XI0+2*XI2+4*XI1);
        
    end
elseif true
    fprintf('Error: method ');
    fprintf(method)
    fprintf(' has not been programmed yet.\n');
    estimate = 0;
    %%
end
end

