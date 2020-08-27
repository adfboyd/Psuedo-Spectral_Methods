function E  =  tanhmove()

%Simulates a moving tanhwave and adapts the grid as it moves. 

close all;  

N = 50;
nn = 1000;
xx = linspace(-1,1, nn)';


A = 40; %Steepness of tanh map


initMap = @tanmap;

c = -0.2; %Initial wave position
initArgs = [A/2 c -1 1];

args = initArgs;

ChebPoints = chebx(N);
BaseInterp = baryM(ChebPoints, xx);

x = initMap(ChebPoints, args);


tnh = @(x,c) tanh(A*(x-c));
y = tnh(x, c);

tol = 10^-4;  % If this is smaller than e-5 it is very unstable unless it's right in the middle
c1 = 1-tol;
c2 = 1+tol;  %These seem to need to be much smaller in this case, in comparison to the paper, and the gaussian.

initNorm = sobnorm(y, initMap, args);
ratios = [];
norms = [];
argstored = [];
errors = [];
errors2 = [];
xxInterp = baryM(x, xx);
count = 0;
crange = linspace(c, -c, 401);
figure('position', [50 200 700 500]);

for i = crange
    c = i;
    tnhc = @(x) tnh(x,c);
    y1 = tnhc(x);
    Norm = sobnorm(y1, @tanmap, args);
    
    ratio = Norm/initNorm;
    ratios = [ratios;ratio];
%    disp(['Ratio = ', num2str(ratio)]);
    if ratio < c1 || ratio > c2
        count = count+1;
%       
%         %This
%         is the version with the built in optimiser. From some investigation it seems very sensitive to initial
%         guesses, so I have left in a bit where it finds the max of the
%         derivative to guess beta. This helps
%         I have also included a new line in abfind() within argfinder2,
%         this attempts to penalise too large alpha appearing. It weirdly
%         dodoesn't seem to work very well.
        argnew = argfinderAlphaFixed(x, y1, N, nn, args);
%       This is a version with Alpha fixed. This does seem to work well but
%       is of course limited to this kind of problem, and a good first
%       guess. It takes in the old args as it's first guess for the new
%       ones, except it tries to find the max of the derivative to estimate
%       beta.
        args = argnew;
        argstored = [argstored;args];
        disp('New args found:');
        disp(args)
        disp(['When c = ', num2str(c)]);
        xnew = tanmap(ChebPoints, args);
        InterpM = baryM(x, xnew);
        ynew = InterpM*y1;
        x = xnew;
        y1 = ynew;
        xxInterp = baryM(x, xx);
        Normnew = sobnorm(y1, @tanmap, args);
        initNorm = Normnew;
        norms = [norms;Normnew];
    else
        norms = [norms;Norm];
    end
    
    y2 = tnhc(ChebPoints);
    
    hold off
    plot(x, y1, 'b.', xx, xxInterp*y1, 'b--');%, x(N/2+1), y1(N/2 +1), 'ro', args(2), 0, 'bo');
    %Also plots midpoint, which should be beta
    hold on
    plot(ChebPoints, y2 , 'r.', xx, BaseInterp*y2, 'r--');
    plot(xx, tnhc(xx));
    legend('Adaptive Points', 'Adaptive Interpolation', 'Fixed Points', 'Fixed Interpolation', 'Exact Solution');
    grid on
    ylim([-1.2,1.2]);
    xlim([-1,1]);
    text(0.1, 0, ['C = ', num2str(c)]);
    xlabel('x');
    ylabel('f(x)');
    title('Moving Tanh Interpolation');
    pause(0.001);
    error = norm(tnhc(xx) - xxInterp*y1, 1);
    error2 = norm(tnhc(xx) - BaseInterp*y2, 1);
    
    errors = [errors;error];
    errors2 = [errors2;error2];
    
end
E = max(errors);
figure('position', [750 200 700 500], 'Colormap',[1 0 0; 0 0 1; 1 1 0]);
% figure('Colormap',[1 0 1; 0 0 1; 1 1 0])
hold off
semilogy(crange,errors);
hold on
grid on
semilogy(crange, errors2);
xlabel('C');
ylabel('Error');
title(['Error in interpolation of tanh(', num2str(A), '(x-C)) as C varies. N = ', num2str(N)]);
legend('Adaptive Grid','Fixed Grid');
disp(['Number of Grid changes = ', num2str(count)]);
end
