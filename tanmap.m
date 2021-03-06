function [X,Dx,DDx] = tanmap(s,args)

%Implements a tanmap, and returns its first and second derivative at every
%point.
%S is Chebyshev points in [-1,1], and args = [Alpha, beta, A, B]
gamma = args(1);
beta = args(2);
A = args(3);
B = args(4);
% if beta<A || beta>B
%     return
% end
%[gamma, beta, A, B] = args;
scale = (B-A)/2;
mid = (B+A)/2;
beta = ( beta - mid ) / scale ;
k = atan(gamma*(1+beta))/atan(gamma*(1-beta));
s0 = (k-1)/(k+1);
lambda = atan(gamma*(1-beta))/(1-s0);
X = scale*(beta + tan(lambda*(s-s0))/gamma)+mid;
Dx = scale*(lambda/gamma)*(cos(lambda*(s-s0)).^(-2));
DDx = scale*(2*lambda^2/gamma)*(cos(lambda*(s-s0)).^(-2)).*tan(lambda*(s-s0));
return