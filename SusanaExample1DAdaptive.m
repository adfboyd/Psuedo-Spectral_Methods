function SusanaExample1DAdaptive()
%Runs the 1D Pedestrian Dynamics Equation from the project.
    close all;
    
    N = 50;
    nn = 2001;
    xx = linspace(-1,1,nn)';
    
    
    
    initMap = @tanmap;
    ibeta = -0.5;
    initArgs = [5 ibeta -1 1];
    
    
    tol = 10^-1;
    c1 = 1-tol;
    c2 = 1+tol;
    
    args = initArgs;
    
    ChebPoints = chebx(N);
%    [y,Int] = ChebyshevPoints(N-1);
    x = initMap(ChebPoints, args);
    
    
    vMax = 1.5;
    Sigma = (0.05)^2;
    
    a = 0.5;
    b = 0.6;


%     a = 0.9;
%     b = 0.975;

%     a = 0.5;
%     b = 0.5;
    
    rho_ic = InitialCondition(x);
    y0 = rho_ic;
    
    
    IndLeft = 1;
    IndRight = N+1;
   
    
    argnew = argfinderBetaFixed(x, y0, N, nn, args, 20);
    args = argnew;
    disp(['Found starting args = ', num2str(args)]);
    x = initMap(ChebPoints, args);
    Dy = diffm(N, initMap, args);
    div = Dy;
    grad = Dy;
    rho_ic = InitialCondition(x);
    x0 = x;
    y0 = rho_ic;
    
    initNorm = sobnorm(y0, initMap, args);
    disp(['Init norm = ', num2str(initNorm)]);
    
    t0 = 0;
    tMax = 0.6;
    nTimes = 20;

    mM    = ones(size(x));
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
    argsout = [];
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
        
        A = repmat(args, nt-1, 1);
        
        argsout = [argsout;A];
 %       disp(X);
 
        
        y1 = Rho_t(nt,:)';
%         disp(['Y = ']);  disp([num2str(y1)]);
        argnew = argfinderAlphaFixed(x, y1, N, nn, args);
        args = argnew;
        
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

    Nfixed = 10*N;
    x2 = chebx(Nfixed);
    
    Dy2 = diffm(Nfixed, @linmap, [-1,1]);
    div2 = Dy2;
    grad2 = Dy2;
    
    IndRight2 = Nfixed+1;
    IndLeft2 = 1;
    
    rho_ic2 = InitialCondition(x2);
   
    
    mM    = ones(size(x2));
    mM([1,Nfixed+1]) = 0; % make the equation a DAE with algebraic ones at 1 and N
    opts = odeset('RelTol',10^-6,'AbsTol',10^-6,'Mass',diag(mM));
    
    
    tic
    [~, Rho_t2] = ode15s(@rhs2, Times_out, rho_ic2, opts);
    toc
    
%     disp([size(Times_out);size(t2)])
    
    %%%PLOTTING--------------------------------
    
%     figure();
%     plot(x0, rho_ic,'.', x2, rho_ic2,'o')
    
%     pause;
%     figure();
%     plot(Times_out);
%     pause;
    
    
    figure('Position', [50, 200, 700, 500])
    intFix = baryM(x2, xx);
%     disp([size(Times_out);size(t2)]);
    for iTimes = 1:length(Times_out)
        Time = Times_out(iTimes);
        xT = xout(iTimes, :)';
        argsT = argsout(iTimes,:);
        beta = argsT(2);
        rhot = Rho_out(iTimes,:)';
        rhot2 = Rho_t2(iTimes,:)';
        %rhoE = ExactSolution(xT,outTimesorig(iTimes));
        
        intM = baryM(xT, xx);
        
%         rhoE = ExactSolution(xT,Time);
%         rhoE2 = ExactSolution(x2, Time);
%         
%         yy = ExactSolution(xx, Time);
        
        hold off
        plot(xx, intFix*rhot2, 'r')%, 'LineWidth', 3);
        hold on
        plot(xx, intM*rhot, 'b--')% 'LineWidth', 3);
        plot(xT,rhot,'b.', 'MarkerSize', 20);
%         plot(x2, rhot2, 'r.');
        
        plot(beta, 0.5, 'ko', 'LineWidth', 5);
%         plot(xx,yy,'--k');
% %         err(iTimes) = norm((abs(rhot - rhoE)) / (abs(rhoE)),1);
        
%         e1 = norm(rhot - rhoE, 1);
%         errors1 = [errors1;e1];
%         
%         e2 = norm(rhot2 - rhoE2, 1);
%         errors2 = [errors2;e2];
        
        ylim([-0.1,0.6]);
        xlim([-1,1]);
        xlabel('x');
        ylabel('\rho (x)');
        title('Solution of Pedestrian Dynamics Equation');
        legend('Fixed Grid Interpolation, N = 500', 'Adaptive Interpolation', 'Adaptive, N = 50', '\beta estimate');
        text(-0.8, 0.1,['Time = ', num2str(Time)]);
        pause(0.5)
        
    end
    
    
    %ERRORS------------------------------------(N/A)
%     figure('Position', [800 200 600 500]);
%     
%     loglog(t2, errors1, 'b');
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
        
        flux = getFlux(rho);
        
        dydt =  - div*flux;

        dydt(IndRight)  = flux(IndRight) - b*rho(IndRight);
        dydt(IndLeft)   = -flux(IndLeft) + a*(1-rho(IndLeft));
        
        dydt = dydt(:);
                
    end

    function f = getFlux(rho)
        f = - ( Sigma*grad*rho - vMax .* rho.*(1-rho) );
        %f = - ( Sigma*grad*rho - vMax .* rho.*exp(-5*rho) );
    end

    function dydt = rhs2(t,rho)
        
        flux = getFlux2(rho);
        
        dydt =  - div2*flux;

        dydt(IndRight2)  = flux(IndRight2) - b*rho(IndRight2);
        dydt(IndLeft2)   = -flux(IndLeft2) + a*(1-rho(IndLeft2));
        
        dydt = dydt(:);
                
    end

    function f = getFlux2(rho)
        f = - ( Sigma*grad2*rho - vMax .* rho.*(1-rho) );
        %f = - ( Sigma*grad*rho - vMax .* rho.*exp(-5*rho) );
    end


    

    
    
    function [value, isterminal, direction] = event(t,y)
        Norm  = sobnorm(y, @tanmap, args);
        ratio = Norm/initNorm;
   %     disp(ratio);
        value = (ratio>c1) & (ratio<c2);
        
        isterminal = 1;
        direction = 0;
    end

    function rho0 = InitialCondition(y)
   
        rho0 = 0.33*(-tanh(50*(y-ibeta))+1)/2;
    end

end