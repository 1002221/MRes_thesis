function a=sum_nonlinearterm_2D(u,j,d1,d2,d3,i,w1,w2,eps)
temp1=0;
for o1=1:length(w1)
    for o2=1:length(w2)
        temp1=temp1+sum(eps.*u(i,d1(w1(o1)),d1(w2(o2)),j-1)...
        .*u(i,d2(w1(o1)),d2(w2(o2)),j-1).*u(i,d3(w1(o1)),d3(w2(o2)),j-1));   
    end
end
a=temp1;