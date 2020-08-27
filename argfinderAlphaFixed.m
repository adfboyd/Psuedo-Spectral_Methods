function args = argfinderAlphaFixed(X, fx, N, nn, start_args)
%Calculates optimal new arguments for the tanmap given that the value for
%alpha should not change.
function y = interpy(x, fx, xx)
    IntM = baryM(x, xx);
    y  = IntM*fx;
        
end
        
f = @(z)interpy(X, fx, z);

left = start_args(3);
right = start_args(4);
chebpoints = chebx(N);
xx = linspace(left, right, nn)';

[s, ~, ~] = linmap(chebpoints, [left right]);
fs = f(s);
D = diffm(N, @linmap, [left right]);
df = D * fs;



interp = baryM(s, xx);
dfinterp = interp*df;

[~, dfmaxi] = max(abs(dfinterp));
xmax = xx(dfmaxi);
beta = xmax;
% init_alpha = 5; % Initial alpha guess
% alpha = init_alpha;

alpha = start_args(1);
% beta = start_args(2);

ab = beta;


ab_min = fminsearch(@abfind, ab); %This built in routine should find optimal alpha and beta(roughly) given initial guess.
% ab_min = fminbnd(@abfind, start_args(3), start_args(4));
beta = ab_min;

%     disp(['Alpha = ' num2str(alpha)])
%     disp(['Beta = ' num2str(beta)])

args = [alpha beta -1 1];


function Sob = abfind(ab)
%    disp(num2str(alpha))
    beta = ab(1);
    args = [alpha beta -1 1];
    [x, ~, ~] = tanmap(chebpoints, args);
    Fx = f(x);
    
    Sob = sobnorm(Fx, @tanmap, args);
%     if alpha>=20
%         disp(['with Old alpha = ', num2str(Sob)]);
%         
%         Sob = 20000;
%         disp(['With new alpha = ', num2str(Sob)]);
%     end

end
end
