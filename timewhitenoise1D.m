function xi=timewhitenoise1D(m,t,randoms_init)
randoms=randoms_init(1:2^m+1,1);
a=linspace(0,1,2^m+1);
b=0;
for i=1:size(a,2)
b=b+i*(a(i)<=t && a(i+1)>t);
end
xi=2^(m/2)*randoms(b);
