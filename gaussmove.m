function [E, argmax] = gaussmove()
%Simulates a spreading out gaussian, and forces the function to regrid
%appropriately.

N = 50; 
nn = 1000;
xx = linspace(-1,1,nn)';


c = 0;  %mean of Gaussian

initMap = @tanmap;
initArgs = [40 c -1 1];

args = initArgs;

ChebPoints = chebx(N);
BaseInterp = baryM(ChebPoints,xx);

x = initMap(ChebPoints, args);

sig = 0.005;
t0 = sig;
D0 = 1;%Change to small number eg 0.01 and make upper bound of linspace large at start of for loop to swap directions
Gaussian = @(x, sig) exp(-((x-c)/sig).^2);

y = Gaussian(x, sig);
y = ExactSolution(x, t0);

tol = 10^-1;
c1  = 1-tol;
c2 = 1+tol;

initNorm = sobnorm2(y, initMap, args, 1/sig);
ratios = [];
norms = [];
argstored = [];
errors = [];
errors2 = [];
xxInterp = baryM(x, xx);
count = 0;
sigrange = linspace(sig, 1, 501);
figure('position', [50, 200, 700, 500]);
for i = sigrange
    sig = i;
    GaussC = @(x)Gaussian(x, sig);
    y1 = GaussC(x);
    Norm = sobnorm2(y1, @tanmap, args, 1/sig);
    
    ratio = Norm/initNorm;
    ratios = [ratios;ratio];
    
    if ratio < c1 || ratio > c2
        count = count +1;
        args(1) = 0.8*args(1);
        
%        argnew = argfinderwint(x, y1', N, nn, 2, [-1 1], [1 10]);
        
%        argnew = argfinderwint2(x, y1, N, nn, 2, [-1 1], args);
        
        argnew = argfinderBetaFixed(x, y1, N, nn, args, 1/sig);
%         Similarly to in tanhmove, this is a version of argfinder that uses the initial beta
%         and doesn't change it. This seems to perform best. Again good
%         initial guesses are important, but are easier in this case.
        args = argnew;
        args(1) =  abs(args(1)); % When going from wide to steep it seems to obtain negative values of alpha for some reason
        %This seems to fix it
        
        argstored = [argstored;args];
        disp('New args found:');
        disp(args)
        disp(['When sig = ', num2str(sig)]);
        xnew = tanmap(ChebPoints, args);
        InterpM = baryM(x, xnew);
        ynew = InterpM*y1;
        x = xnew;
        y1 = ynew;
        xxInterp = baryM(x, xx);
        Normnew = sobnorm2(y1, @tanmap, args, 1/sig);
        initNorm = Normnew;
        norms = [norms;Normnew];
    else
        norms = [norms;Norm];
    end
    
    y2 = GaussC(ChebPoints);
    
    
    hold off
    plot(x, y1, '.', xx, xxInterp*y1)
    hold on
    plot(ChebPoints, y2, '.', xx, BaseInterp*y2);
    ylim([-1 1]);
    xlim([-1 1]);
    pause(0.001);
    error = norm(GaussC(xx) - xxInterp*y1, 1);
    errors = [errors;error];
    
    error2 = norm(GaussC(xx) - BaseInterp*y2, 1);
    errors2 = [errors2;error2];
%     pause;
end
[E, argmax]  = max(errors);
%disp(argmax);
figure('position', [750 200 750 500]);
hold off
loglog(sigrange,errors);
hold on 
grid on
semilogy(sigrange, errors2);
legend('Adaptive Grid', 'Fixed Grid');
title('Error in interpolation as Gaussian spreads out');
xlabel('Variance');
ylabel('Error');
disp(['Number of Grid changes = ', num2str(count)]);



  function rho = ExactSolution(y,t)
        % exact solution to diffusion equation on unbounded domain
        rho = sqrt(t0/t)*exp(-y.^2/(4*D0*t));
    end


end