
for k=1:850
    [m,n]=size(X(k).h0_t20_n);
    H0_BD(1:m,k:k+1) = X(k).h0_t20_n;
end