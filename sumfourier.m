function a=sumfourier(N,lower,upper,x,u)
temp=0;
for n=lower:upper
        temp= temp+ u(n)*exp(1i*(n-N-1).*x);
end  
a=temp;
end