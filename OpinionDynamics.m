function OpinionDynamics()
% Implements the opinions dynamics equation studied in the project
    close all;
    
    N = 100;
    nn = 1001;
    
    left = 0;
    right = 1;
    xx = linspace(left,right,nn)';
    
    mult = 1.1;
    sig = 0.05;
    sigma = sig/sqrt(2*pi);
    ymax = 30;
    R = 0.2;
    
    initMap = @tanmap;
    initArgs = [2 (right-left)/2 left right];
    
    
    tol = 10^-0.5;
    c1 = 1-tol;
    c2 = 1+tol;
    
    args = initArgs;
    
    ChebPoints = chebx(N);
%    [y,Int] = ChebyshevPoints(N-1);
    x = initMap(ChebPoints, args);
    
    
    rho_ic_UN = UnnormIC(x);
    
   
    
    
%     plot(x, rho_ic, '.');
%     figure();
    y0 = rho_ic_UN;
    
    argnew = argfinderBetaFixed(x, y0, N, nn, args, 20);
    args = argnew;
    disp(['Found starting args = ', num2str(args)]);
    x = initMap(ChebPoints, args);
    Dy = diffm(N, initMap, args);
    rho_ic_UN = UnnormIC(x);
    M_Conv = ConvM(N, initMap, args, R);
    
    w = integw2(N, initMap, args);
    C = w*rho_ic_UN;
    
    rho_ic = rho_ic_UN/C;
    
    
    y0 = rho_ic;
    
    
    initNorm = sobnorm(y0, initMap, args);
    disp(['Init norm = ', num2str(initNorm)]);
    
    t0 = 0;
    tMax = 2*pi*1;
    nTimes = 100;

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
        argnew = argfinderBetaFixed(x, y1, N, nn, args,5);
        if argnew(1)>=args(1)/2
            
            args = argnew;
        end
        
        disp(['Break ', num2str(count), ', Time reached = ', num2str(t(nt)/(2*pi)), ', New args found = ', num2str(args)]);
        args(1) = abs(args(1));
        xnew = tanmap(ChebPoints, args);
        InterpM = baryM(x, xnew);
        ynew = InterpM*y1;
        initNorm = sobnorm(ynew, @tanmap, args);
        x = xnew;
        y0 = ynew;
        Dy = diffm(N, @tanmap, args);
        M_Conv = ConvM(N, @tanmap, args, R);
        t0 = t(nt);
        outTimes = t0:timestep:tMax;
    end
    toc
    
    %FIXED GRID---------------------------------
%     disp(size(Times_out));
    tic 
    
    Nfixed = 500;
    Dy2 = diffm(Nfixed, @linmap, [left,right]);
    x2 = linmap(chebx(Nfixed), [left,right]);
    M_Conv2 = ConvM(Nfixed, @linmap, [left,right], R);
    
    rho_ic2 = UnnormIC(x2)/C;
    
    mM    = ones(size(x2));
    mM([1,Nfixed+1]) = 0; % make the equation a DAE with algebraic ones at 1 and N
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    
    
    [t2, Rho_t2] = ode15s(@rhs2, Times_out, rho_ic2, opts);
    toc
    
    %FIXED GRID SMALL N---------------------------
    
    
    tic 
    
    Nfixed = N;
    Dy2 = diffm(Nfixed, @linmap, [left,right]);
    x3 = linmap(chebx(Nfixed), [left,right]);
    M_Conv2 = ConvM(Nfixed, @linmap, [left,right], R);
    
    rho_ic2 = UnnormIC(x3)/C;
    
    mM    = ones(size(x3));
    mM([1,Nfixed+1]) = 0; % make the equation a DAE with algebraic ones at 1 and N
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    
    
    [t2, Rho_t3] = ode15s(@rhs2, Times_out, rho_ic2, opts);
    toc
    
%     disp([size(Times_out);size(t2)])
    
    %%%PLOTTING--------------------------------
    
%     errors1=[];
%     errors2=[];
    
    figure('Position', [50, 200, 700, 500])
    intFix = baryM(x2, xx);
    intFix2 = baryM(x3, xx);
    for iTimes = 1:length(t2)
        Time = Times_out(iTimes);
        xT = xout(iTimes, :)';
        rhot = Rho_out(iTimes,:)';
        rhot2 = Rho_t2(iTimes,:)';
        rhot3 = Rho_t3(iTimes,:)';
        %rhoE = ExactSolution(xT,outTimesorig(iTimes));
        
        IntM = baryM(xT, xx);
%         rhoE = ExactSolution(xT,Time);
%         rhoE2 = ExactSolution(x2, Time);
%         
%         yy = ExactSolution(xx, Time);
%         
        hold off
        plot(xT,rhot,'b.', 'MarkerSize', 20);
        hold on
        plot(xx, IntM*rhot, 'b', 'Linewidth', 2);
        plot(x3, rhot3, 'g.', 'MarkerSize', 20);
        plot(xx, intFix2*rhot3, 'g', 'Linewidth', 2);
%         plot(x2, rhot2, 'r.');
        plot(xx, intFix*rhot2, 'r--', 'Linewidth', 2);
        legend('Adaptive','Adaptive Interpolation', 'Fixed', 'Fixed Interpolation', 'Exact');%, 'Exact Interpolation');
        text(0.2, ymax/2, ['Time = ', num2str(round(Time/(2*pi),2)), 's']);
%         plot(xx,yy,'--k');
% %         err(iTimes) = norm((abs(rhot - rhoE)) / (abs(rhoE)),1);
%         
%         e1 = norm(rhot - rhoE, 1);
%         errors1 = [errors1;e1];
%         
%         e2 = norm(rhot2 - rhoE2, 1);
%         errors2 = [errors2;e2];
        
        ylim([-ymax/5,ymax]);
        xlim([left,right]);
        xlabel('x');
        ylabel('\rho (x)');
        title(['Solution of OD equation, \sigma = ', num2str(sig)]);
        pause(0.0001)
        
        
    end
    
    
   
%     figure('Position', [800 200 600 500]);
%     
%     semilogy(t2, errors1, 'b');
%     hold on
%     loglog(t2, errors2, 'r');
%     grid on
%     legend('Adaptive Grid', 'Fixed Grid');
%     xlabel('t');
%     ylabel('Error');
%     title('Error in diffusion equation')
%     disp(Times_out);
    
    
    %----------------------------------------------------------------------
    
    function dydt = rhs(t,rho)
        
        
        flux = getFlux(rho); % -d_x \rho
        dydt =  Dy*flux; % d_xx \rho
        
        % Neumann problem

        dydt(N+1)  = flux(N+1) - 0; % match to exact boundary condition
        dydt(1)  = flux(1) - 0;  % match to exact boundary condition
        
        
%         dydt(N+1)  = rho(N+1) - ExactSolution(x(N+1),t-delay); % match to exact boundary condition
%         dydt(1)  = rho(1) - ExactSolution(x(1),t-delay);  % match to exact boundary condition

        dydt = dydt(:);
                
    end

    function f = getFlux(rho)
%         M_Conv = ConvM(N, @tanmap, args);
%         disp([size(rho);size(M_Conv);size(Dy);size(sigma^2/2)]);

        f1 = rho .*( M_Conv * rho ) ;
        f2 = (sigma^2/2)*Dy*rho;
        f = f1 + f2;
    end

    function dydt = rhs2(t, rho)
        
        flux = getFlux2(rho);
        dydt = Dy2*flux;
        
        % Neumann problem

        dydt(Nfixed+1)  = flux(Nfixed+1) - 0; % make d_x \rho = 0 on right endpoint
        dydt(1)   = flux(1) - 0;  % make d_x \rho = 0 on left endpoint


        % Dirichlet problem
%         dydt(N)  = rho(Nfixed+1) - 0; % make \rho = 0 on right endpoint
%         dydt(1)   = rho(1) - 0;  % make \rho = 0 on left endpoint
    end

    function f = getFlux2(rho)
        f = rho .* (M_Conv2 * rho) + sigma^2/2*Dy2*rho;
    end
    
    function [value, isterminal, direction] = event(~,y)
        Norm  = sobnorm(y, @tanmap, args);
        ratio = Norm/initNorm;
 
        value = (ratio>c1) & (ratio<c2);
        isterminal = 1;
        direction = 0;
    end
    
    function rho = UnnormIC(y)
        rho = exp(-20*(y-(right-left)/2).^2);
    end

    


end