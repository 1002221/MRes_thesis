function answer=test3(r,J,m,T,N,x,y,randoms_init)
%J=1; %number of iterations
%m=1; %we subdivide the unit interval into 2^m equally spaced points
%T = 1; %units of time
%N=1; %subscripts of the Fourier coefficients go from -N to N
%epsilon=.01;
u = zeros(2^m*T+1,2*N+1,2*N+1,J);
%randoms=randn(2*N+1,2);
%randoms= randn(2*N+1,1); %standard normal random variables, one for each n
eps = 1; %parameter by which we multiply u^3. Allows us to get rid of the 
%non-linear term if we wish
source = zeros(2*N+1,2*N+1,2^m+1); %source term (i.e. the 'f')
% source(N+1,N+1)=1;

for n=N+2:2*N+1
    for l=N+2:2*N+1
         for i=1:2^m*T+1
            xi1=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,1));
            xi2=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,2));
            xi3=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,2));
            xi4=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,2));
            source(n,:)=1/(2*pi)*(xi1pos-1i*xi2pos); 
            source(2*(N+1)-n,:)=1/(2*pi)*(xi1neg+1i*xi2neg);
         end
    end
end

source(N+1,N+1)=randoms(1,1,1)*randoms(1,1,2);
lower=max(1,N+1-7); %compute the truncated upper and lower bounds for the 
%sum outside the loop for efficiency's sake
upper=min(2*N+1,N+1+7);
[d1, d2, d3] = ndgrid(lower:upper, lower:upper, lower:upper);
soln=zeros(J,1);
for j=2:J
    for k=1:2^m*T+1
        for l=1:2*N+1
            for n=1:2*N+1       
                q = find(d1+d2+d3-(N+1)*3==n-(N+1));
                q2= find(d1+d2+d3-(N+1)*3==l-(N+1));
                temp = 0;
                for i=1:k
                    %truncate at -100 and 100 (take reindexing into account)
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
                        (source(n,l)-temp1);
                end        
                u(k,n,l,j)=(1/(2^m))*temp;
            end 
        end
    end
        %find the solution at each iteration by summing over the Fourier transforms
%     for n=lower:upper
%         soln(j) = soln(j)+ u(2^m*T+1,n,l,j)*exp(1i*(n-N-1)*x)*exp(1i*(l-N-1)*y);
%     end
end
u;
%find solution at (T, x) but summing over the Fourier transforms. Truncate
%at -100, 100 (take re-indexing into account)
for n=lower:upper
    for l=lower:upper
        soln(J) = soln(J)+ u(2^m*T+1,n,l,J)*exp(1i*(n-N-1)*x)*exp(1i*(l-N-1)*y);
    end
end
answer=soln(J);
% plot(linspace(1,J,J),soln);



%NOTE: by setting 
%source = zeros(2*N+1,2*N+1,1);
%source(N+1,N+1)=1; 
%we eliminate the
%dependence on space. Thus, the problem becomes u' = 1-u^3; u(0)=0. This
%can be solved via separation of variables, and we get that u(:,N+1,N+1,J) 
%should equal to 0.823041 (for T=1, that is). Use the code 
%'exactsolution.m' to find this.