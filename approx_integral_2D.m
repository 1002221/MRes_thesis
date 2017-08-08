function a=approx_integral_2D(u,j,d1,d2,d3,w1,w2,eps,m,N,k,n,l,source)
temp=0;
for i=1:k
    temp1=sum_nonlinearterm_2D(u,j,d1,d2,d3,i,w1,w2,eps);
    temp = temp+ exp(((n-(N+1))^2+(l-(N+1))^2)*(i-1-k-1)/(2^m))*...
        (source(n,l,i)-temp1);     
end        
a=(1/(2^m))*temp;
end