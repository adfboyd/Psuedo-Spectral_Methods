function [x,Dx,DDx] = linmap(s, args)
%implements a linear map from s in [-1,1] to x in [A, B], where args =
%[A,B]
A = args(1);
B = args(2);
N = length(s);
scale = (B-A)/2;
mid = (B+A)/2;
x = scale*s + mid;
Dx = scale*ones(size(x));
DDx = zeros(size(x));