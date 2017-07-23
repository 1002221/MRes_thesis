function xi = xi(epsilon,x)
n=10;
k=linspace(-4/epsilon,4/epsilon,2*n+1);
randoms = randn(2*n+1);
xi=sum(rhohat(epsilon*k).*randoms.*exp(1i*k*x));
