function Dm = diffm(N, map, args)
%Computes an appropriately transformed Differentiation matrix
[D,s] = cheb(N);
[x, Dx, DDx] = map(s, args);
Dm = diag(1./Dx)*D;
