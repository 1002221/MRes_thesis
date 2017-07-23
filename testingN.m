maxN=2;
solns=zeros(maxN,1);
randoms_init=randn(101,2);
Ns=linspace(1,maxN,maxN);
%randoms=randn(N+1,2);
for i=1:maxN
    i
    solns(i)=test2(.01,6,5,1,Ns(i),1,randoms_init);
    solns(i)
end
plot(Ns,solns);
%set(gca,'xdir','reverse')
%always converges to something at 10^(-2)

%something to think about: maybe you should be able to start with
%an arbitrarily large random matrix, and then only select the part
%you're interested in, so you can compare up to what point increasing
%N actually helps

%beyond N=7 there's no advantage, for epsilon=.1 and J=6.

%but, even with a smaller epsilon, it appears that there may still be no
%advantage at taking a massive N. Got to think about this more carefully.
-0.622851055387264
 -0.138700777141528 - 0.140868150546134i
 
 %so the various modes actually stay the same. Is this what I'd expect?
 %Wouldn't I expect them to change? I mean, if you're summing to 1 and you
 %take more numberes, then surely you have more things that can sum to 1?
 
 %in fact, they do change. They get more accurate as you go on - and the
 %more N's you put in, the more accurate it gets - but after a certain
 %point, the 'marginal' modes are so small they don't actually make a
 %difference.
 
 %but if that's the only place where we actually use epsilon, then what
 %difference does it make? we might as well just multiply by 1, if we're
 %gonna use it like this. Doesn't make any difference. It only makes a
 %difference if N gets massively big - but even then, it doesn't matter.