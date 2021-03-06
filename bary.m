%An illustration of the Interpolation without the matrix, as seen within
%the Barycentric Interpolation paper by Lloyd N. Trefethen.
n = 50;
fun = inline('abs(x) + .5*x- x.^2');
x = cos(pi*(0:n)'/n);
f = fun(x);
c = [1/2; ones(n-1,1); 1/2].*(-1).^((0:n)');
xx = linspace(-1,1,500)';
numer = zeros(size(xx));
denom = zeros(size(xx));
exact = zeros(size(xx));
for j = 1:n+1
    xdiff = xx-x(j);
    temp = c(j)./xdiff;
    numer = numer + temp*f(j);
    denom = denom + temp;
    exact(xdiff==0) = 1;
end
ff = numer./denom;
jj = find(exact);
ff(jj) = f(exact(jj));
plot(x,f,'.',xx,ff,'-')

