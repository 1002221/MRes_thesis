function a=sourceterm(N,m,randoms_init,r,T)
source = zeros(2*N+1,2^m+1);
for n=N+2:2*N+1
     for i=1:2^m*T+1
        x1=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,1));
        x2=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,2));
        %terms where n>0
        source(n,i)=1/(2*sqrt(pi))*(x1-1i*x2);
        %terms where n<0
        source(2*(N+1)-n,i)=conj(source(n,i));
     end
end
%take care of the term where n=0, which has a different normalising
%constant
for i=1:2^m*T+1
    x1=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,N+1,1));
    source(N+1,i)=1/(sqrt(2*pi))*(x1); 
end
a=source;