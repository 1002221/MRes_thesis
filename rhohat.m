function rhohat=rhohat(y)
syms x w 
F(w) = (2*pi)^(-.5)*fourier(exp(-x^2/2)); %finds fourier transform for 
%function f(x) = sqrt(2 pi) e^(-x^2/2), where the constant has been added
%so that rhohat(0)=1.
rhohat = double(F(y)); %convert symbolic result to number