N=10;
soln=0;
for n=max(1,(2*N+1)/2-100):min(2*N+1,(2*N+1)/2+100)
    soln = soln+ rhohat(n)*exp(1i*(n-N-1)*1);
end