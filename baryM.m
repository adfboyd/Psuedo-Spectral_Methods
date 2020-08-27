function A = baryM(x,xx)
%Computes the interpolation matrix from initial grid x to final grid xx
n = length(x);
nn = length(xx);

w = (-1).^(0:(n-1))';
w(1) = w(1)/2;
w(n) = w(n)/2;

[sorted, indices] = sort([x(2:n);xx(1:nn)]);
sorted(indices) = cumsum(indices<=n-1)+1;

index = sorted(n:n+nn-1);

maska = (x(index)==xx);
maskb = (maska == 0);
xxx = xx(maskb);

denom = zeros(length(xxx),1);
A     = zeros(nn,n);

if(~isempty(xxx))
    if(length(xxx)==1)
        temp = w./(xxx-x);
        A(maskb,:) = temp;
        denom = denom + sum(temp);
    else
        for i=1:n    %go through all points in x
            temp = w(i)./(xxx-x(i));        % DNS: Paper p507: temp = c(j)./xdiff
            A(maskb,i) = temp;
            denom = denom + temp;
        end
    end
end
% DNS: This loop is like in the paper mentioned above, p506-507
% c = [1/2; ones(n-1,1); 1/2].*(-1).^((0:n)'); --- which is basically the
% same as w given above, but with n=n-1 in c. Then xx is defined, then
% numer = zeros(size(xx)); --- similar
% denom = zeros(size(xx)); --- similar
% for j = 1:n+1; xdiff = xx-x(j); temp = c(j)./xdiff --- Same essentially
% numer = numer + temp*f(j); --- this line isn't there
% denom = demon + temp; end; --- this is

A(maskb,:)  =  diag(1./denom) *  A(maskb,:) ;

for i=1:length(maska)
    if(maska(i)==1)
        A(i,index(i)) = 1;
    end
end
        
        