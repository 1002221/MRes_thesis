function answer=test3(r,J,m,T,N,randoms_init)
%USAGE
%J:= number of iterations
%m:= we subdivide the unit interval into 2^m equally spaced points
%T:= units of time
%N:= subscripts of the Fourier coefficients go from -N to N
%epsilon:= parameter with which we smooth the white noise (=1/r)
%randoms_init:= random variables supplied to the function: for example,
%randn(2^10,30,4)
u = zeros(2^m*T+1,2*N+1,2*N+1,J);
eps = 1; %parameter by which we multiply u^3. Allows us to get rid of the 
%non-linear term if we wish
source = zeros(2*N+1,2*N+1,2^m+1); %source term (i.e. the 'f')
%set up space domain
x=linspace(0,2*pi,2^2);
y=linspace(0,2*pi,2^2);
for n=N+2:2*N+1
    for l=N+2:2*N+1
         for i=1:2^m*T+1
            xi1=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,1));
            xi2=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,2));
            xi3=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,3));
            xi4=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,4));
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
        xi1=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,1));
        xi3=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,3));
        %terms where n>0, l=0
        source(n,N+1,i)=1/(sqrt(2)*pi)*(xi1-1i*xi3); 
        %terms where n<0, l=0
        source(2*(N+1)-n,N+1,i)=1/(sqrt(2)*pi)*(xi1+1i*xi3);
     end
end
for l=N+2:2*N+1
     for i=1:2^m*T+1
        xi1=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,1));
        xi2=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,2));
        %terms where n=0, l>0
        source(N+1,l,i)=1/(sqrt(2)*pi)*(xi1-1i*xi2); 
        %terms where n=0, l<0
        source(N+1,2*(N+1)-l,i)=1/(sqrt(2)*pi)*(xi1+1i*xi2);
     end
end
%take care of term where n,l=0
for i=1:2^m*T+1
        xi1=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,1));
        source(N+1,N+1,i)=1/(4*pi)*(xi1); 
end
cutoff=10;
lower=max(1,N+1-cutoff); %compute the truncated upper and lower bounds for the 
%sum outside the loop for efficiency's sake
upper=min(2*N+1,N+1+cutoff);
%the next step will be used when finding triples of values that sum to a
%particular value
[d1, d2, d3] = ndgrid(lower:upper, lower:upper, lower:upper);
soln=zeros(J,size(x,2),size(y,2));
h = waitbar(0,'Please wait...');
for j=2:J
    for k=1:2^m*T+1
        for l=1:2*N+1
            for n=1:2*N+1       
                q = find(d1+d2+d3-(N+1)*3==n-(N+1));
                q2= find(d1+d2+d3-(N+1)*3==l-(N+1));
                temp = 0;
                for i=1:k
                    temp1=0;
                    for o=1:size(q)
                        for o2=1:size(q2)
                            temp1=temp1+sum(eps...
                                .*u(i,d1(q(o)),d1(q2(o2)),j-1)...
                                .*u(i,d2(q(o)),d2(q2(o2)),j-1)...
                                .*u(i,d3(q(o)),d3(q2(o2)),j-1));
                        end
                    end
                    temp = temp+ exp(((n-(N+1))^2+(l-(N+1))^2)...
                        *(i-1-k-1)/(2^m*T))*...
                        (source(n,l,i)-temp1);
                end        
                u(k,n,l,j)=(1/(2^m))*temp;
            end 
        end
    end
    %find the solution at each iteration by summing over the Fourier 
    %transforms
%     for n=lower:upper
%         soln(j) = soln(j)+ u(2^m*T+1,n,l,j)*exp(1i*(n-N-1)*x)*...
%         exp(1i*(l-N-1)*y);
%     end
    waitbar((j-1) / (J-1),h)
end
close(h)
%find solution at (T, x,y) but summing over the Fourier transforms.
for n=lower:upper
    for l=lower:upper
        soln(J,:,:) = soln(J,:,:)+ u(2^m*T+1,n,l,J)*exp(1i*(n-N-1).*x)...
        .*exp(1i*(l-N-1).*y);
    end
end
answer=soln(J,:,:);
% plot(linspace(1,J,J),soln);

%NOTE: by setting 
%source = zeros(2*N+1,2*N+1,1);
%source(N+1,N+1)=1; 
%we eliminate the
%dependence on space. Thus, the problem becomes u' = 1-u^3; u(0)=0. This
%can be solved via separation of variables, and we get that u(:,N+1,N+1,J) 
%should equal to 0.823041 (for T=1, that is). Use the code 
%'exactsolution.m' to find this.