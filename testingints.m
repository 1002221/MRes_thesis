%let's start by trying to approximate simpler integrals that we can
%compute. Like, say, x^2 between 0 and 2. Should get 8/3.

%If we follow the same method we have here, it should be: sum the
%integrand, then divide by however many points you've used. Let's try this.
m=10;
T=2;

ints=zeros(2^m*T+1,1);
for k=1:2^m*T+1
    temp=0;
    for i=1:k
    temp = temp+ ((i-1)/(2^m))^2;
    end        
    ints(k)=(1/(2^m+1))*temp; %approximate the integral by 
    %dividing by the number of points in the partition
end
ints
%ok, so this is correctly approximating 1/3, which is the integral between
%0 and 1.
   
%ok, this works.