function a=approx_integral(u,j,d1,d2,d3,w,eps,m,N,k,n,source)
temp=0;
for i=1:k
    temp1=sum_nonlinearterm(u,j,d1,d2,d3,i,w,eps);
    temp = temp+ exp((n-(N+1))^2*((i-1)-(k-1))/(2^m))*...
        (source(n,i)-temp1);       
end        
a=(1/(2^m))*temp;
end
