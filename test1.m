function soln=test1(epsilon,J,m,T,N,x0,randoms)
%J=1; %number of iterations
%m=1; %we subdivide the unit interval into 2^m equally spaced points
%T = 1; %units of time
%N=1; %subscripts of the Fourier coefficients go from -N to N
%epsilon=.01;
u = zeros(2^m*T+1,2*N+1,J);
%randoms=randn(2*N+1,2);
%randoms= randn(2*N+1,1); %standard normal random variables, one for each n
eps = 1; %parameter by which we multiply u^3. Allows us to get rid of the 
%non-linear term if we wish
source = zeros(2*N+1,1); %source term (i.e. the 'f')
%source(N+1)=1;
%randoms=randn(2*N+1,1); %standard normal RVs that will get multiplied by
%the Fourier transforms of the source term
for n=N+2:2*N+1
    source(n)=rhohat(epsilon*(n-N-1))*.5*(randoms(n-(N+1),1)-1i*...
        randoms(n-(N+1),2)); 
    %compute the 
    %Fourier transform, but take re-indexing into account
end
for n=1:N
    source(n)=rhohat(epsilon*(n-N-1))*.5*(randoms(N+1-n,1)+...
        1i*randoms(N+1-n,2)); %compute the 
    %Fourier transform, but take re-indexing into account
end
source(N+1)=randoms(N+1,1);
lower=max(1,N+1-100);
upper=min(2*N+1,N+1+100);
for j=1:J 
    for k=1:2^m*T+1
        for n=1:2*N+1
            temp = 0;
            for i=1:k
                temp1 = 0;
                %truncate at -100 and 100 (take reindexing into account)
                for x=lower:upper 
                    for y=lower:upper
                        for z=lower:upper
                            if (x+y+z-(N+1)*3==n-(N+1)) %reindexing
                            temp1 = temp1 + eps*u(i,x,j)*u(i,y,j)*...
                                u(i,z,j);
                            if (temp1 ~= 0)
                                temp1;
                            end
                            end
                        end
                    end
                end
                temp = temp+ exp((n-(N+1))^2*(i-1-k-1)/(2^m*T))*...
                    (source(n)-temp1);
                
            end
            u(k,n,j)=(1/(2^m*T+1))*temp;
        end
    end
end
%u
%find solution at (T, x) but summing over the Fourier transforms. Truncate
%at -100, 100 (take re-indexing into account)
soln=0;
for n=lower:upper
    soln = soln+ u(2^m*T+1,n,J)*exp(1i*(n-N-1)*x0);
end
%soln=u(2^m*T+1,:,J);
%setting source = zeros(2*N+1,1); %source term (i.e. the 'f') and 
%source(N+1)=1; we eliminate the
%dependence on space. Thus, the problem becomes u' = 1-u^3; u(0)=0. This
%can be solved via separation of variables, and we get that u(1) should
%equal to 0.823041

