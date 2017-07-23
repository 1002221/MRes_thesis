function rhohat2D=rhohat2D(y1,y2)
syms x w 
F(w) = (2*pi)^(-.5)*fourier(exp(-x^2/2));
rhohat2D = double(F(norm([y1,y2])));