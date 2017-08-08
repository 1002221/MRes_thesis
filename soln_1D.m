function a=soln_1D(r,J,m,T,N,randoms_init)
%USAGE
%J:= number of iterations
%m:= we subdivide the unit interval into 2^m equally spaced points
%T:= units of time
%N:= subscripts of the Fourier coefficients go from -N to N. Because of
%Matlab's indexing rules, 1 corresponds to -N, and N+1 to 0.
%epsilon:= parameter with which we smooth the white noise (=1/r)
%randoms_init:= random variables supplied to the function: for example,
%randn(2^10,30,4)

x=linspace(0,2*pi,2^5); %set up space domain

u = zeros(2^m*T+1,2*N+1,J); %this will contain the Fourier coefficients. 
eps = 1; %set to 0 to take out non-linear term from code

source=sourceterm(N,m,randoms_init,r,T);
%compute the truncated upper and lower bounds for N in the sum

cutoff=10;
lower=max(1,N+1-cutoff); 
upper=min(2*N+1,N+1+cutoff);

soln=zeros(J,size(x,2)); %this will contain the solutions at each iteration

[d1, d2, d3] = ndgrid(lower:upper, lower:upper, lower:upper);
q=sum_to_n(N,d1,d2,d3);%finding triples that sum to each n <= N

h = waitbar(0,'Please wait...');
for j=2:J %carry out the numerical scheme described in the project
%     waitbar((j-1) / (J-1),h)
    for k=1:2^m*T+1
        u=approx_soln(u,q,j,d1,d2,d3,eps,m,N,k,source);
        waitbar(((j-1-1)*(2^m*T+1)+k) / ((J-1)*(2^m*T+1)),h)
    end
    soln(j,:)=sumfourier(N,lower,upper,x,u(2^m*T+1,:,j));
end
close(h)

surf(x,linspace(2,J,J-1),real(soln(2:J,:)));%plot soln at each iteration
xlabel('space')
ylabel('iteration')
zlabel('solution')

a=soln(J,:);%return solution at last iteration

%NOTE: by setting 
%source = zeros(2*N+1,1);
%source(N+1)=1; 
%we eliminate the
%dependence on space. Thus, the problem becomes u' = 1-u^3; u(0)=0. This
%can be solved via separation of variables, and we get that u(:,N+1,J) 
%should equal to 0.823041 (for T=1, that is). Use the code 
%'exactsolution.m' to find this.