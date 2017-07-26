function a=soln_1D(r,J,m,T,N,randoms_init)
%USAGE
%J:= number of iterations
%m:= we subdivide the unit interval into 2^m equally spaced points
%T:= units of time
%N:= subscripts of the Fourier coefficients go from -N to N
%epsilon:= parameter with which we smooth the white noise (=1/r)
%randoms_init:= random variables supplied to the function: for example,
%randn(2^10,30,4)
x=linspace(0,2*pi,2^2); %set up space domain
u = zeros(2^m*T+1,2*N+1,J); %this will contain the Fourier coefficients. 
%The first argument is the position in time, the second is the (reindexed) 
%index of the Fourier coefficient, and the last is the iteration number.
%Entry N+1 corresponds to the zeroth Fourier coefficient.
eps = 1; %parameter by which we multiply u^3. Allows us to get rid of the 
%non-linear term if we wish, by setting it to 0

source = zeros(2*N+1,2^m+1);
for n=N+2:2*N+1
     for i=1:2^m*T+1
        x1=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,1));
        x2=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,n,2));
        %terms where n>0
        source(n,i)=1/(2*pi)*(x1-1i*x2);
        %terms where n<0
        source(2*(N+1)-n,i)=1/(2*pi)*(x1+1i*x2);
     end
end
%take care of the term where n=0
for i=1:2^m*T+1
    x1=timewhitenoise1D(r,(i-1)/(2^m),randoms_init(:,N+1,1));
    source(N+1,i)=1/(2*pi)*(x1); 
end
%compute the truncated upper and lower bounds for the 
%sum outside the loop for efficiency's sake. There doesn't seem to be any
%increase in accuracy with taking more than 7 fourier coefficients
cutoff=10;
lower=max(1,N+1-cutoff); 
upper=min(2*N+1,N+1+cutoff);
%the next step will be used when finding triples of values that sum to a
%particular value
[d1, d2, d3] = ndgrid(lower:upper, lower:upper, lower:upper);
soln=zeros(J,size(x,2));
h = waitbar(0,'Please wait...');
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
                    (source(n,i)-temp1);
            end        
            u(k,n,j)=(1/(2^m))*temp; %approximate the integral by 
            %dividing by the difference between successive points in the
            %partition
        end    
    end
%     for n=lower:upper
%         soln(j,:) = soln(j,:)+ u(2^m*T+1,n,j)*exp(1i*(n-N-1).*x);
%     end
    waitbar((j-1) / (J-1),h)
end
close(h)
% *** comment out next part of the code to see how the approximation 
% improves with each iteration ***
% surf(x,linspace(2,J,J-1),real(soln(2:J,:)));
% xlabel('space')
% ylabel('iteration')
% zlabel('solution')
for n=lower:upper
    soln(J,:) = soln(J,:)+ u(2^m*T+1,n,J)*exp(1i*(n-N-1).*x);
end
a=soln(J,:);

%NOTE: by setting 
%source = zeros(2*N+1,1);
%source(N+1)=1; 
%we eliminate the
%dependence on space. Thus, the problem becomes u' = 1-u^3; u(0)=0. This
%can be solved via separation of variables, and we get that u(:,N+1,J) 
%should equal to 0.823041 (for T=1, that is). Use the code 
%'exactsolution.m' to find this.
