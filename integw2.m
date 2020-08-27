function w = integw2(N, Map, args)
%Returns the integration weights for a given set of points
[D,s] = cheb(N);
[x, Dx, ~] = Map(s, args);

i = 2:N+1;
Dx = Dx(i);
x = x(i);
Di = inv(D(i, i));
w = Dx'.*Di(N, :);
w = [0,w];
