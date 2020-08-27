function M_conv = ConvM(N, Map, args, R)
%Computes a convolution matrix as discussed in use for the OpinionDynamics
%eqn.
argL = length(args);
left = args(argL-1);
right = args(argL);

% R = 0.2;

ChebPts = chebx(N);
y = Map(ChebPts, args);

M_conv = zeros(N+1, N+1);
f = @(x)x;

for i = 1:N+1
    
    y0 = y(i);
    
    yMin = max(left, y0-R);
    yMax = min(right, y0+R);
    
    subPts = linmap(chebx(N), [yMin, yMax]);
    
    IP_i = baryM(y, subPts);
    
    Int_i = integw2(N, @linmap, [yMin, yMax]);
    
    fP = f(y0 - subPts);
    
    M_conv(i, :) = (Int_i.*fP')*IP_i;
    
end