function a=approx_soln(u,q,j,d1,d2,d3,eps,m,N,k,source)
for n=N+1:2*N+1
    w=q(:,n-N);%find triples that sum to n. First column contains those
    %that sum to 0, second to 1, etc...
    w=w(w~=0);
    u(k,n,j)=approx_integral(u,j,d1,d2,d3,w,eps,m,N,k,n,source);
    u(k,2*(N+1)-n,j)=conj(u(k,n,j));          
end
a=u;
end