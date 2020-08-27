function DDm = ddiffm(N, map, args)
%Computes an appropriately transformed 2nd order differentiation matrix.
[D, s] = cheb(N);
[~, Dx, DDx] = map(s, args);
Dm = diag(1./Dx)*D;
D2 = D^2;
DDm = (D2 - D*diag(1./Dx)*diag(DDx))*diag(1./Dx.^(2));
