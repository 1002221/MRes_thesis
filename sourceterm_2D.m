function a=sourceterm_2D(N,m,T,r,randoms_init)
source = zeros(2*N+1,2*N+1,2^m+1); %source term (i.e. the 'f')
for n=N+2:2*N+1
    for l=N+2:2*N+1
         for i=1:2^m*T+1
            xi1=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,l,1));
            xi2=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,l,2));
            xi3=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,l,3));
            xi4=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,l,4));
            %terms where n>0, l>0
            source(n,l,i)=1/(4*pi)*(xi1-1i*xi2-1i*xi3-xi4); 
            %terms where n<0, l<0
            source(2*(N+1)-n,2*(N+1)-l,i)=1/(4*pi)*(xi1+1i*xi2+1i*xi3-xi4);
            %terms where n>0, l<0
            source(n,2*(N+1)-l,i)=1/(4*pi)*(xi1+1i*xi2-1i*xi3+xi4);
            %terms where n<0, l>0
            source(2*(N+1)-n,l,i)=1/(4*pi)*(xi1-1i*xi2+1i*xi3+xi4);
         end
    end
end
for n=N+2:2*N+1
     for i=1:2^m*T+1
        xi1=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,1));
        xi3=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,3));
        %terms where n>0, l=0
        source(n,N+1,i)=1/(sqrt(2)*pi)*(xi1-1i*xi3); 
        %terms where n<0, l=0
        source(2*(N+1)-n,N+1,i)=1/(sqrt(2)*pi)*(xi1+1i*xi3);
     end
end
for l=N+2:2*N+1
     for i=1:2^m*T+1
        xi1=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,1));
        xi2=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,2));
        %terms where n=0, l>0
        source(N+1,l,i)=1/(sqrt(2)*pi)*(xi1-1i*xi2); 
        %terms where n=0, l<0
        source(N+1,2*(N+1)-l,i)=1/(sqrt(2)*pi)*(xi1+1i*xi2);
     end
end
%take care of term where n,l=0
for i=1:2^m*T+1
        xi1=timewhitenoise(r,(i-1)/(2^m),randoms_init(:,n,1));
        source(N+1,N+1,i)=1/(4*pi)*(xi1); 
end
a=source;