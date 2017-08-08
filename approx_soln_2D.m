function a=approx_soln_2D(u,q1,q2,j,d1,d2,d3,eps,m,N,k,source)
for l=1:2*N+1
    for n=1:2*N+1
        w1=q1(:,n);%find triples that sum to n
        w1=w1(w1~=0);
        w2=q2(:,n);%find triples that sum to n
        w2=w2(w2~=0);
        u(k,n,l,j)=approx_integral_2D(u,j,d1,d2,d3,w1,w2,eps,m,N,k,n,l,source);
    end
end
a=u;
end