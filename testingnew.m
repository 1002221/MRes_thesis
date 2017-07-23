m=4;
temp=0;
t=.6;
for i=1:2^m
    temp=temp+2^m*randn(m)*(i/(2^m)<=t)*(t<=(i+1)/(2^m)<=t);
end