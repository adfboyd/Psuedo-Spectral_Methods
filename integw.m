function w = integw(N, Map, args)
%Returns integration weights for a given set of points.
[D,s] = cheb(N);
[x, Dx, ~] = Map(s, args);

i = 2:N+1;
Dx = Dx(i);
x = x(i);
Di = inv(D(i, i));
w = Dx'.*Di(N, :);


