function a=sumfourier_2D(N,lower,upper,x,y,u,j,m,T)
temp=0;
for n=lower:upper
    for l=lower:upper
        temp = temp+ u(2^m*T+1,n,l,j)*...
            exp(1i*(n-N-1)*x)*exp(1i*(l-N-1)*y);
    end
end
a=temp;
end