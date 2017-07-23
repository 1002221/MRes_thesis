function a=test2(epsilon,J,m,T,N,randoms_init)
%USAGE
%J:= number of iterations
%m:= we subdivide the unit interval into 2^m equally spaced points
%T:= units of time
%N:= subscripts of the Fourier coefficients go from -N to N
%epsilon:= parameter with which we smooth the white noise
%randoms:= random variables supplied to the function: should be at least
%(N+1) times 1 in size.

x=linspace(0,2*pi,2^5);
u = zeros(2^m*T+1,2*N+1,J); %this will contain the Fourier coefficients. 
%The first argument is the position in time, the second is the (reindexed) 
%index of the Fourier coefficient, and the last is the iteration number.
%Entry N+1 corresponds to the zeroth Fourier coefficient.
eps = 1; %parameter by which we multiply u^3. Allows us to get rid of the 
%non-linear term if we wish, by setting it to 0
source = zeros(2*N+1,1); %preallocate matrix for the source term 
%(i.e. the 'f'). First, we populate its terms corresponding to a positive
%index
% source(N+1)=1; 
randoms=randoms_init(1:N+1,:);
for n=N+2:2*N+1
    source(n)=rhohat(epsilon*(n-N-1))*.5*(randoms(n-(N+1)+1,1)-1i*...
        randoms(n-(N+1)+1,2)); 
    
end
%Next, we populate its terms corresponding to a negative index. As
%rhohat(-x) = rhohat(x), then because of the re-indexing, 
%rhohat(epsilon*(n-N-1)) = rhohat(epsilon*(N+1-n))
for n=1:N
    source(n)=rhohat(epsilon*(n-N-1))*.5*(randoms(N+1-n+1,1)+...
        1i*randoms(N+1-n+1,2)); 
end
%Finally, the term corresponding to 0
source(N+1)=randoms(1,1); 
%compute the truncated upper and lower bounds for the 
%sum outside the loop for efficiency's sake. There doesn't seem to be any
%increase in accuracy with taking more than 7 fourier coefficients
lower=max(1,N+1-7); 
upper=min(2*N+1,N+1+7);
[d1, d2, d3] = ndgrid(lower:upper, lower:upper, lower:upper);
soln=zeros(J,size(x,2));
for j=2:J
    for k=1:2^m*T+1
        for n=1:2*N+1
            %we want the sum of the indices of the u's at the previous
            %iteration to equal the current index. This is what we ensure
            %here
            q = find(d1+d2+d3-(N+1)*3==n-(N+1));
            temp = 0;
            for i=1:k
                temp1=0;
                for o=1:size(q)
                    temp1=temp1+sum(eps.*u(i,d1(q(o)),j-1)...
                        .*u(i,d2(q(o)),j-1).*u(i,d3(q(o)),j-1));
                end
                temp = temp+ exp((n-(N+1))^2*((i-1)-(k-1))/(2^m))*...
                    (source(n)-temp1);
            end        
            u(k,n,j)=(1/(2^m))*temp; %approximate the integral by 
            %dividing by the difference between successive points in the
            %partition
        end    
    end
%     comment out to find the solution at each iteration by summing over 
%     the Fourier coefficients
    for n=lower:upper
        soln(j,:) = soln(j,:)+ u(2^m*T+1,n,j)*exp(1i*(n-N-1).*x);
    end
end
surf(x,linspace(2,J,J-1),real(soln(2:J,:)));
xlabel('space')
ylabel('iteration')
zlabel('solution')
% for n=lower:upper
%     soln(J,:) = soln(J,:)+ u(2^m*T+1,n,J)*exp(1i*(n-N-1).*x);
% end
a=soln(2,:);
%7 iterations seem to be enough

%NOTE: by setting 
%source = zeros(2*N+1,1);
%source(N+1)=1; 
%we eliminate the
%dependence on space. Thus, the problem becomes u' = 1-u^3; u(0)=0. This
%can be solved via separation of variables, and we get that u(:,N+1,J) 
%should equal to 0.823041 (for T=1, that is). Use the code 
%'exactsolution.m' to find this.
