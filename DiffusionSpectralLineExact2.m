function E = DiffusionSpectralLineExact2()
%Runs an adaptive version of a diffusion equation with exact BCs.
    close all;
    
    N = 50;
    nn = 1000;
    xx = linspace(-1,1,nn)';
    mult = 0.85;
    
    
    initMap = @tanmap;
    initArgs = [20 0 -1 1];
    
    
    tol = 10^-0.5;
    c1 = 1-tol;
    c2 = 1+tol;
    
    args = initArgs;
    
    ChebPoints = chebx(N);
%    [y,Int] = ChebyshevPoints(N-1);
    x = initMap(ChebPoints, args);
    
        
    D0 = 1; % diffusion coefficient
    t0 = 0.0001; % initial time
    initT = t0;
    rho_ic = ExactSolution(x,t0);   
    
%     plot(x, rho_ic, '.');
%     figure();
    y0 = rho_ic;
    
    argnew = argfinderBetaFixed(x, y0, N, nn, args, 1/t0);
    args = argnew;
    disp(['Found starting args = ', num2str(args)]);
    x = initMap(ChebPoints, args);
    Dy = diffm(N, initMap, args);
    rho_ic = ExactSolution(x, t0);
    
    y0 = rho_ic;
    
    initNorm = sobnorm(y0, initMap, args);
    disp(['Init norm = ', num2str(initNorm)]);
    
    
    tMax = 1;
    nTimes = 80;

    mM    = ones(size(x));
%     mM(1) = 0; % make the equation a DAE with algebraic ones at 1 and N
    mM([1,N+1]) = 0; % make the equation a DAE with algebraic ones at 1 and N
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM), 'event', @event);
%    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    

    %ADAPTIVE GRID--------------------------------------
    
    Rho_out = [];
    Times_out = [];
    teout = [];
    yeout = [];
    ieout = [];
    xout=[];
    count= 0;
    
    timestep = (tMax-t0)/nTimes;
    outTimes = t0:timestep:tMax;
   
    tic
    while t0 < tMax-timestep
        count = count+1;
        
        [t,Rho_t, te, ye, ie] = ode15s(@rhs,outTimes,y0,opts);
        
        nt = length(t);
%         disp(['Nt = ', num2str(nt)]);

        Rho_out = [Rho_out;Rho_t(1:nt-1,:)];
        Times_out = [Times_out;t(1:nt-1)];

        teout = [teout;te];
        yeout = [yeout;ye];
        ieout = [ieout; ie];
        
        X = repmat(x', nt-1, 1);
        
        xout = [xout;X]; %Store the grids used for each timestep, for plotting purposes.
        
 %       disp(X);
 
        
        y1 = Rho_t(nt,:)';
%         disp(['Y = ']);  disp([num2str(y1)]);
        args(1) = mult*args(1);
%         argnew = argfinderBetaFixed(x, y1, N, nn, 1, args,1/t0);
% %         argnew = argfinderBetaFixed(x, y1, N, nn, 1, args);
%         args = argnew;
        
        disp(['Break ', num2str(count), ', Time reached = ', num2str(t(nt)), ', New args found = ', num2str(args)]);
        args(1) = abs(args(1));
        xnew = tanmap(ChebPoints, args);
        InterpM = baryM(x, xnew);
        ynew = InterpM*y1;
        initNorm = sobnorm(ynew, @tanmap, args);
        x = xnew;
        y0 = ynew;
        Dy = diffm(N, @tanmap, args);
        t0 = t(nt);
        outTimes = t0:timestep:tMax;
    end
    toc
    
    %FIXED GRID---------------------------------
%     disp(size(Times_out));

    Nfixed = N;
    Dy2 = diffm(Nfixed, @linmap, [-1,1]);
    x2 = chebx(Nfixed);
    rho_ic2 = ExactSolution(x2, initT);
   
    mM    = ones(size(x2));
    mM([1,Nfixed+1]) = 0; % make the equation a DAE with algebraic ones at 1 and N
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    
    
    tic
    [t2, Rho_t2] = ode15s(@rhs2, Times_out, rho_ic2, opts);
    toc
    
%     disp([size(Times_out);size(t2)])
    
    %%%PLOTTING--------------------------------
    
    errors1=[];
    errors2=[];
    
    errors1i = [];
    errors2i = [];
    
    figure('Position', [50, 200, 600, 500])
    IntFix = baryM(x2, xx);
    for iTimes = 1:length(t2)
        Time = Times_out(iTimes);
        xT = xout(iTimes, :)';
        rhot = Rho_out(iTimes,:)';
        rhot2 = Rho_t2(iTimes,:)';
        %rhoE = ExactSolution(xT,outTimesorig(iTimes));
        
        IntM = baryM(xT, xx);
        
        rhoE = ExactSolution(xT,Time);
        rhoE2 = ExactSolution(x2, Time);
        
        yy = ExactSolution(xx, Time);
        
        hold off
        plot(xT,rhot,'b.');
        hold on
        plot(xx, IntM*rhot, 'b');
        plot(x2, rhot2, 'r.');
        plot(xx, IntFix*rhot2, 'r');
        plot(xx,yy,'--k');
%         err(iTimes) = norm((abs(rhot - rhoE)) / (abs(rhoE)),1);
        legend('Adaptive Solution', 'Adaptive Interpolation', 'Fixed Solution', 'Fixed Interpolation', 'Exact Solution');
        e1 = norm(rhot - rhoE, 1);
        errors1 = [errors1;e1];
        
        e2 = norm(rhot2 - rhoE2, 1);
        errors2 = [errors2;e2];
        
        e1i = norm(IntM*rhot - yy, 1);
        errors1i = [errors1i;e1i];
        
        e2i = norm(IntFix*rhot2 - yy,1);
        errors2i = [errors2i;e2i];
        
        ylim([-0.2,1.1]);
        xlim([-1,1]);
        grid on;
        xlabel('x');
        ylabel('u(x)');
        title('Solution of Diffusion Equation');
        text(-0.8, 0.3, ['Time = ', num2str(Time), 's']);
        pause(initT/Time)
        
    end
    
    
   
    figure('Position', [800 200 600 500]);
    
    loglog(t2, errors1, 'b');
    hold on
    loglog(t2, errors2, 'r');
    loglog(t2, errors1i, 'b--');
    loglog(t2, errors2i, 'r--');
    grid on
    legend('Adaptive Grid', 'Fixed Grid', 'Adaptive Interpolation', 'Fixed Interpolation');
    xlabel('t');
    ylabel('Error');
    title('Error in diffusion equation')
%     disp(Times_out);
    E = log10(max(errors1i));
    
    %----------------------------------------------------------------------
    
    function dydt = rhs(t,rho)
        
        
        flux = getFlux(rho); % -d_x \rho
        dydt =  - Dy*flux; % d_xx \rho
        
        % Neumann problem

        dydt(N+1)  = rho(N+1) - ExactSolution(x(N+1),t); % match to exact boundary condition
        dydt(1)  = rho(1) - ExactSolution(x(1),t);  % match to exact boundary condition
        
        
%         dydt(N+1)  = rho(N+1) - ExactSolution(x(N+1),t-delay); % match to exact boundary condition
%         dydt(1)  = rho(1) - ExactSolution(x(1),t-delay);  % match to exact boundary condition

        dydt = dydt(:);
                
    end

    function f = getFlux(rho)
        f = - D0*Dy*rho;
    end

    function dydt = rhs2(t, rho)
        
        flux = getFlux2(rho);
        dydt = -Dy2*flux;
        
        % Neumann problem

        dydt(Nfixed+1)  = rho(Nfixed+1) - ExactSolution(x2(Nfixed+1), t); % make d_x \rho = 0 on right endpoint
        dydt(1)   = rho(1) - ExactSolution(x2(1), t);  % make d_x \rho = 0 on left endpoint


        % Dirichlet problem
%         dydt(N)  = rho(N) - 0; % make \rho = 0 on right endpoint
%         dydt(1)   = rho(1) - 0;  % make \rho = 0 on left endpoint
    end

    function f = getFlux2(rho)
        f = -D0*Dy2*rho;
    end
    
    function [value, isterminal, direction] = event(~,y)
        Norm  = sobnorm(y, @tanmap, args);
        ratio = Norm/initNorm;
 
        value = (ratio>c1) & (ratio<c2);
        isterminal = 1;
        direction = 0;
    end

    function rho = ExactSolution(y,t)
        % exact solution to diffusion equation on unbounded domain
        rho = sqrt(initT/t)*exp(-y.^2/(4*D0*t));
    end


end