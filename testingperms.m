N=1;
n=3;
no=0;
for x=max(1,(2*N+1)/2-100):min(2*N+1,(2*N+1)/2+100) 
    for y=max(1,(2*N+1)/2-100):min(2*N+1,(2*N+1)/2+100)
        for z=max(1,(2*N+1)/2-100):min(2*N+1,(2*N+1)/2+100)
            if (x+y+z-(N+1)*3==n-(N+1)) %reindexing
            x
            y
            z
            end
        end
    end
end
no
[d1, d2, d3] = ndgrid(lower:upper, lower:upper, lower:upper);
d = d1+d2+d3;
i = find(d-(N+1)*3==n-(N+1));
A=[d1(i),d2(i),d3(i)];
size0=size(A,1);
for p=1:size0
A(p,1)
A(p,2)
A(p,3)
end