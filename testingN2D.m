maxN=10;
solns=zeros(maxN,1);
randoms_init=randn(8,8,2);
Ns=linspace(1,maxN,maxN);
%randoms=randn(N+1,2);
for i=1:maxN
    i
    solns(i)=test3(.1,8,5,1,Ns(i),1,1,randoms_init);
    solns(i)
end
plot(Ns,solns);
%set(gca,'xdir','reverse')
%always converges to something at 10^(-2)

%something to think about: maybe you should be able to start with
%an arbitrarily large random matrix, and then only select the part
%you're interested in, so you can compare up to what point increasing
%N actually helps

%beyond N=7 there's no advantage