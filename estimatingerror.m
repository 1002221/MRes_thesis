epsilons = logspace(0,-2,5); %increasingly small values of epsilon
solns=zeros(size(epsilons,2),1); %preallocate matrix where we will store
%the solutions (u_epsilon)(t,x) for a given (t,x) and the increasingly
%small values of epsilon
randoms=randn(8,2); %fix a matrix of random variables to be used for the 
%whole exercise
N=7;
for i=1:size(epsilons,2)
    i %check progress as you run the code
    solns(i)=test2(epsilons(i),6,5,1,N,1,randoms);
end
semilogx(epsilons,solns); %plot (u_epsilon)(t,x) as a function of epsilon
set(gca,'xdir','reverse') %reverse direction of epsilons so they go from
%large to small

