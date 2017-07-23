x = linspace(-10,10,201);
randoms = randn(2*n+1);
epsilon = .01;
xi = 0;
for i=1:2*n+1
    xi = xi+rhohat(epsilon*k)*randoms(i)*exp(i*k*x);
end;
    