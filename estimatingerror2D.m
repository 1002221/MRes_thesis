epsilons = linspace(5,10,6);
solns=zeros(size(epsilons,2),1);
% N=6;
% randoms=randn(N+1,N+1,2);
for i=1:size(epsilons,2)
    i
    solns(i)=soln_2D(epsilons(i),3,5,1,5,randoms_init);
end
figure
plot(epsilons,solns);
% set(gca,'xdir','reverse')
%always converges to something at 10^(-2)