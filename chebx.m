function x = chebx(N)
%Returns (N+1) chebyshev points in [-1,1]

if N==0; x = 1; return, end
x = -cos(pi*(0:N)/N)'; 
