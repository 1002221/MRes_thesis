function a=soln_2D(r,J,m,T,N,randoms_init)
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
%set up space domain
x=.5;
y=.5;
source=sourceterm_2D(N,m,T,r,randoms_init);
cutoff=10;
lower=max(1,N+1-cutoff); %compute the truncated upper and lower bounds for the 
%sum outside the loop for efficiency's sake
upper=min(2*N+1,N+1+cutoff);
%the next step will be used when finding triples of values that sum to a
%particular value
[d1, d2, d3] = ndgrid(lower:upper, lower:upper, lower:upper);
q1=sum_to_n_2D(N,d1,d2,d3);%finding triples that sum to each n <= N
q2=sum_to_n_2D(N,d1,d2,d3);%finding triples that sum to each n <= N

soln=zeros(J);

h = waitbar(0,'Please wait...');
for j=2:J
    for k=1:2^m*T+1
        u=approx_soln_2D(u,q1,q2,j,d1,d2,d3,eps,m,N,k,source);%implement 
        %numerical scheme described in the report - see called functions
        %and report for details
        waitbar(((j-1-1)*(2^m*T+1)+k) / ((J-1)*(2^m*T+1)),h)
    end
    %find the solution at each iteration by summing over the Fourier 
    %transforms
    soln(j)=sumfourier_2D(N,lower,upper,x,y,u,j,m,T);
end
 close(h)
u(2^m*T+1,:,:,J);
a=soln(J);
% plot(linspace(2,J,J-1),soln(2:J));

%NOTE: by setting 
%source = zeros(2*N+1,2*N+1,1);
%source(N+1,N+1)=1; 
%we eliminate the
%dependence on space. Thus, the problem becomes u' = 1-u^3; u(0)=0. This
%can be solved via separation of variables, and we get that u(:,N+1,N+1,J) 
%should equal to 0.823041 (for T=1, that is). Use the code 
%'exactsolution.m' to find this.

%try taking 3 iterations, but re-run it with smaller and smaller epsilons.
%shouldn't converge

%also check, that in the 1D case, reducing epsilon further doesn't afffect
%the limiting solution