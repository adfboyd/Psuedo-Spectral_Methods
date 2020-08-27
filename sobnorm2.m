function I = sobnorm2(y0, Map, Args, Steepness)
%Calculates the Sobolev norm, with an additional input of Steepness,
%allowing for 'tuning' of the constants A,B and C



N = length(y0)-1;

% arglength = length(Args);
% left = Args(arglength-1);
% right = Args(arglength);
% 
% nn = 1001;
% xx = linspace(left,right, nn)';
% alpha = alphabeta(1);
% beta = alphabeta(2);
% args = [alpha, beta, -1,  1];
%args = [-1 1];
args = Args;

%Steepness = 1/0.03;

[x, ~, ~] = Map(chebx(N), args);
A = 0.25/Steepness^4;
B = 1/Steepness^2; %Set up weights, not sure what to put here,
C = 1;          % in the case of tanh(Nx), DDx increases N^2, so tried to account for that.

% IntM = baryM(x, xx);

%rho = tanh(Steepness*x);
%  u  = rho;
u = y0;
% uInt = IntM*u;

function w = weightfn(x)
%    w = (1- x.^2).^(-1/2);
    w = 1;
end

w = weightfn(x);
% wInt = weightfn(xx);

D = diffm(N, Map, args);
%D2 = ddiffm(N, Map, args);
D2 = D^2;

u_s = (D*u);
% u_sInt = IntM*u_s;

u_ss = (D2*u);
% u_ssInt = IntM*u_ss;


% plot(x, u,'.', x, B*u_s.^2, x, A*u_ss.^2)
Integrand = (A*(u_ss).^2 + B*u_s.^2 + C*u.^2);
Integrandw = w.*Integrand;


integrator = integw(N, Map, args);
%disp(Integrandw)
I = integrator*Integrandw(2:N+1);


% Integrand2 = (A*(u_ssInt).^2 + B*u_sInt.^2 + C*uInt.^2);
% Integrandw2 = wInt.*Integrand2;

% plot(xx, Integrandw2,'.', xx,A*(u_ssInt).^2,'.' , xx, B*u_sInt.^2, xx, C*uInt.^2);
% width = (right - left)/(nn-1);

% I2 = sum(Integrandw2(1:nn-1))*width;

end

