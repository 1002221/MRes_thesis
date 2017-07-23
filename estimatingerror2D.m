epsilons = logspace(0,-2,4);
solns=zeros(size(epsilons,2),1);
N=6;
randoms=randn(N+1,N+1,2);
for i=1:size(epsilons,2)
    i
    solns(i)=test3(epsilons(i),7,5,1,N,1,1,randoms);
end
semilogx(epsilons,solns);
set(gca,'xdir','reverse')
%always converges to something at 10^(-2)