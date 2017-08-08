function a=sum_to_n_2D(N,d1,d2,d3)
q=zeros(100,N+1);%preallocate q - remove extra zeros later
for n=1:2*N+1
            f=find(d1+d2+d3-(N+1)*3==n-(N+1));
            q(1:length(f),n) = f;%q's columns don't go from -N to N, but
            %from 0 to N. As Matlab starts from 1,  1 to N+1 -
            %this is why we subtract N.
end
a=q;
end