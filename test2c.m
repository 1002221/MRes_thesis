function answer=test2c(epsilon,J,m,T,N,x,randoms_init)
%USAGE
%J:= number of iterations
%m:= we subdivide the unit interval into 2^m equally spaced points
%T:= units of time
%N:= subscripts of the Fourier coefficients go from -N to N
%epsilon:= parameter with which we smooth the white noise
%randoms:= random variables supplied to the function: should be a 201x1
%matrix or bigger
J=min(J,6); % There doesn't seem to be any
%increase in accuracy with taking more than 6 iterations
u = zeros(2^m*T+1,2*N+1,J); %this will contain the Fourier coefficients. 
%The first argument is the position in time, the second is the (reindexed) 
%index of the Fourier coefficient, and the last is the iteration number
eps = 1; %parameter by which we multiply u^3. Allows us to get rid of the 
%non-linear term if we wish, by setting it to 0
source = zeros(2*N+1,1); %preallocate matrix for the source term 
%(i.e. the 'f'). First, we populate its terms corresponding to a positive
%index
source(N+1)=1; 
% randoms=randoms_init(1:N+1,:);
% for n=N+2:2*N+1
%     source(n)=rhohat(epsilon*(n-N-1))*.5*(randoms(n-(N+1)+1,1)-1i*...
%         randoms(n-(N+1)+1,2)); 
%     
% end
% %Next, we populate its terms corresponding to a negative index. As
% %rhohat(-x) = rhohat(x), then because of the re-indexing, 
% %rhohat(epsilon*(n-N-1)) = rhohat(epsilon*(N+1-n))
% for n=1:N
%     source(n)=rhohat(epsilon*(n-N-1))*.5*(randoms(N+1-n+1,1)+...
%         1i*randoms(N+1-n+1,2)); 
% end
% %Finally, the term corresponding to 0
% source(N+1)=randoms(1,1); 
%compute the truncated upper and lower bounds for the 
%sum outside the loop for efficiency's sake. There doesn't seem to be any
%increase in accuracy with taking more than 7 fourier coefficients
lower=max(1,N+1-7); 
upper=min(2*N+1,N+1+7);
[d1, d2, d3] = ndgrid(lower:upper, lower:upper, lower:upper);
soln=zeros(J,1);
for j=2:J
    temp = 0;
    for k=2:2^m*T+1
        for n=1:2*N+1
            %we want the sum of the indices of the u's at the previous
            %iteration to equal the current index. This is what we ensure
            %here
            q = find(d1+d2+d3-(N+1)*3==n-(N+1));
           
            temp1=0;
            for o=1:size(q)
                temp1=temp1+sum(eps.*u(k-1,d1(q(o)),j-1)...
                    .*u(k-1,d2(q(o)),j-1).*u(k-1,d3(q(o)),j-1));
            end
            temp = temp+ exp((n-(N+1))^2*(k-1-T)/(2^m))*...
                (source(n)-temp1);     
            u(k,n,j)=(1/(2^m+1))*temp; %approximate the integral by 
            %dividing by the number of points in the partition
        end    
    end
    %find the solution at each iteration by summing over the Fourier transforms
%     for n=lower:upper
%         soln(j) = soln(j)+ u(2^m*T+1,n,j)*exp(1i*(n-N-1)*x);
%     end
end
%plot(linspace(1,J,J),soln);
for n=lower:upper
    soln(J) = soln(J)+ u(2^m*T+1,n,J)*exp(1i*(n-N-1)*x);
end
answer=soln(J);
u(2^m*T+1,:,J)
%7 iterations seem to be enough

% %find solution at (T, x) but summing over the Fourier transforms. Truncate
% %at -100, 100 (take re-indexing into account)
% soln=0;
% for n=lower:upper
%     soln = soln+ u(2^m*T+1,n,J)*exp(1i*(n-N-1)*x);
% end
%NOTE: by setting 
%source = zeros(2*N+1,1);
%source(N+1)=1; 
%we eliminate the
%dependence on space. Thus, the problem becomes u' = 1-u^3; u(0)=0. This
%can be solved via separation of variables, and we get that u(1) should
%equal to 0.823041 (for T=1, that is). Wolfram code: 
%1/6(2sqrt(3)tan^(-1)((1+2x)/sqrt(3))-2log(1-x)+log(1+x+x^2))-
%pi/(6*sqrt(3))=1

%this is hard to get my head round. as epsilon gets small, the effect of
%the RV is cancelled out. Ok, so the fact that we get convergence as
%epsilon shrinks isn't at all surprising. What we should do, perhaps, is to
%keep N big enough - but like this, the implementation time will go right
%through the roof. With epsilon = .01, for example, it would require N=100,
%which requires...actually, does it require N=100? Let's do the testing N
%thing, but with a different epsilon.

%why am I evaluating that integral for each k? So confusing.

%ok, maybe we do need to store all the integrals at the intermediate times
%- otherwise, how do we calculate u epsilon (s) at the previous iteration?
