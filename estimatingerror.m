epsilons = linspace(1,9,9); %increasingly small values of epsilon
solns=zeros(size(epsilons,2),2^2); %preallocate matrix where we will store
%the solutions (u_epsilon)(t,x) for a given (t,x) and the increasingly
%small values of epsilon
randoms=randn(2^10,30,4); %fix a matrix of random variables to be used for the 
%whole exercise
N=7;
for i=1:size(epsilons,2)
    i %check progress as you run the code
    solns(i,:)=soln_1D(epsilons(i),10,5,1,N,randoms);
end
surf(linspace(0,2*pi,2^2),epsilons,real(solns)); %plot (u_epsilon)(t,x) as a function of epsilon
%set(gca,'xdir','reverse') %reverse direction of epsilons so they go from
%large to small

