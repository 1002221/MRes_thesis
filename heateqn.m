N = 20;
x = linspace(0, 1, 200);
L = 1;
a = [];
f = zeros(1, length(x));
for k = 1:N
    if(k==1)
 a(k) = 1;
    else
        a(k)=0;
    end
end

%here the initial condition is 1

t = linspace(0, 0.2, 20);
[X, T] = meshgrid(x,t);
U = 0*X;
figure;
plot3([0 1], [0 0], [-1 -1], 'LineWidth', 2);
hold on;
for j = 1:length(t)
 for k = 1:N
 U(j,:) = U(j,:) + a(k)*exp(-(pi*k/L)^2*t(j))*sin(pi*k*x/L);
 end
 plot3(x, t(j)*ones(1,length(x)), U(j,:), 'k', 'LineWidth', 1);
end
S = surf(X,T,U, -U);
set(S, 'EdgeColor', 'none');
colormap(hot);
view([147 4]);
xlabel('x', 'FontSize', 16);
ylabel('y', 'FontSize', 16);
title('Solution of the heat equation for a bar', 'FontSize', 16);