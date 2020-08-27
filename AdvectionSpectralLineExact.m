function RESULT = AdvectionSpectralLineExact()
%Runs the Advection Equation with steep initial condition using an adaptive
%method.
    close all;
    
    N = 100;
    nn = 1000;
    xx = linspace(-1,1,nn)';
    
    
    
    A = 20;
    InitPoint = -0.5;
    initMap = @tanmap;
    initArgs = [A*0.5 InitPoint -1 1];
    
    
    tol = 10^-8; %Best to set this as low as possible such that the solution still converges.
    c1 = 1-tol;  %tol can be reduced as points are added
    c2 = 1+tol;  %tol must be increased when IC is steeper.
    
    args = initArgs;
    
    ChebPoints = chebx(N);
%    [y,Int] = ChebyshevPoints(N-1);
    x = initMap(ChebPoints, args);
    Dy = diffm(N, initMap, args);
        
    C = 1; % wave speed
    t0 = 0; % initial time
    initT = t0;
    
    function y = InitialCondition(x)
        y = tanh(A*(x-InitPoint));
    end

    rho_ic = InitialCondition(x);
    y0 = rho_ic;
    
    initNorm = sobnorm(y0, initMap, args);
%     disp(['Init norm = ', num2str(initNorm)]);
    
    
    tMax = 0.8;
    nTimes = 100;

    mM    = ones(size(x));
    mM(1) = 0; % make the equation a DAE with algebraic one at 1 
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM), 'event', @event);
    
    Rho_out = [];
    Times_out = [];
    teout = [];
    yeout = [];
    ieout = [];
    xout=[];
    count= 0;
    argsout = [];
    
    timestep = (tMax-t0)/nTimes;
    outTimes = t0:timestep:tMax;
%     outTimesorig = outTimes; %Store original outTimes for plotting etc
%     disp(outTimes);



    %%%ADAPTIVE GRID----------------------------
    tic
    while t0 < tMax-timestep*1.5
        count = count+1;
        
        [t,Rho_t, te, ye, ie] = ode15s(@rhs,outTimes,y0,opts);
        
        nt = length(t);
        time_calc = t(nt) - t(1);
        if count < 2                        %%This condition stops the solve if the
            firsttime = time_calc;          %%Regridding starts happening too frequently.
        else 
            if time_calc/firsttime<0.005  
                break
            end
        end
%         disp(['Nt = ', num2str(nt)]);
        Rho_out = [Rho_out;Rho_t(1:nt-1,:)];
        Times_out = [Times_out;t(1:nt-1)];
%         disp(size(Rho_out))
        teout = [teout;te];
        yeout = [yeout;ye];
        ieout = [ieout; ie];
        X = repmat(x', nt-1, 1);
        xout = [xout;X]; %Store the grids used for each timestep, for plotting purposes.
        ARGS = repmat(args, nt-1, 1);
        argsout = [argsout;ARGS];
 %       disp(X);
        
        
        y1 = Rho_t(nt,:)';
%         disp(['Y = ']);  disp([num2str(y1)]);
        argnew = argfinderAlphaFixed(x, y1, N, nn, args);
%         argnew = argfinderBetaFixed(x, y1, N, nn, 1, args);
        args = argnew;
        
        disp(['Break ', num2str(count), ', Time reached = ', num2str(t(nt)), '. New args found = ', num2str(args)]);
       
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
    t1 = toc;
    %%%--------------------------------------------------
    
  %  disp(Times_out);
    
    %%%SIMPLE GRID--------------------------------
    
    Nfixed = N;
    
    Dy2 = diffm(Nfixed, @linmap, [-1 1]);
    y2 = chebx(Nfixed);
    rho_ic2 = InitialCondition(y2);
    
    mM    = ones(size(y2));
    mM(1) = 0; % make the equation a DAE with algebraic one at 1 

    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    tic 
    [~, Rho_t2] = ode15s(@rhs2, Times_out, rho_ic2, opts);
    t2 = toc;
    
    %%%----------------------------------------------
    
    
    %%%PLOTTING%%%---------------------------------------
    
    errors1 = [];
    errors2 = [];
    figure('Position', [100, 300, 600, 500])
    
    errors1i = [];
    errors2i = [];
    
    IntFix = baryM(y2, xx);
    for iTimes = 1:length(Times_out)
        xT = xout(iTimes, :)';  %Grid at this timestep
        rhot = Rho_out(iTimes,:)';  %Adaptive Solution
        rhot2 =  Rho_t2(iTimes,:)';  %Fixed Grid Solution
        iargs = argsout(iTimes, :);  %Args of adaptive grid
        Time = Times_out(iTimes);   %t
%         rhoE = ExactSolution(xT,outTimesorig(iTimes));

        rhoE = ExactSolution(xT,Time);  %Exact solutions in each grid
        rhoE2 = ExactSolution(y2, Time);
        
        yy = ExactSolution(xx, Time);
        IntM = baryM(xT, xx);
        
        hold off
        plot(xx, IntM*rhot, 'b');
        hold on
        plot(xT,rhot,'b.');   %Plot adaptive Solution

%         plot(iargs(2), 0, 'bo'); %Plot Beta
         plot(xx, IntFix*rhot2, 'r');
        plot(y2, rhot2, 'r.'); %Plot fixed grid solution
%         plot(xx, yy);   %Plot exact solution
        
        
        e1 = norm(rhot - rhoE, 1);
        errors1 = [errors1;e1];
        
        e2 = norm(rhot2 - rhoE2, 1);
        errors2 = [errors2;e2];
        
        e1i = norm(IntM*rhot - yy, 1);
        errors1i = [errors1i;e1i];
        
        e2i = norm(IntFix*rhot2 - yy, 1);
        errors2i = [errors2i;e2i];
%         plot(xx, ExactSolution(xx,Time))
        grid on
        legend('Adaptive Interpolation','Adaptive Solution', 'Fixed Grid Interpolation','Fixed Grid Solution');
        text(-0.9, 0.7, ['Time = ', num2str(Time), 's'])
        
        xlabel('x');
        ylabel('u(x)');
        title('Solution of Advection Equation');
        ylim([-1.5,1.5]);
        xlim([-1,1]);
        
        pause(0.001)
        
    end
    
%    disp(argsout);
%     if size(outTimesorig) == size(err)
%         figure()
%         semilogy(outTimesorig,err);
%     else
%         disp([size(outTimesorig);size(err)]);
%         figure();
%         semilogy(outTimesorig,err');
%     end
%     
    figure('Position', [750, 300, 600, 500])
    
    semilogy(Times_out, errors1, 'b');
    hold on
    semilogy(Times_out, errors2, 'r');
    semilogy(Times_out, errors1i, 'b--');
    semilogy(Times_out, errors2i, 'r--');
    grid on
    legend('Adaptive grid', 'Fixed Grid', 'Adaptive Grid Interpolation', 'Fixed Grid Interpolation');
    ylabel('Error');
    xlabel('Time');
    xlim([0, tMax]);
    title('Error in Solution Over Time');
    
    E = max(errors1);
    E2 = max(errors2);
    
    RESULT = [log10(E),log10(E2),t1,t2];
    %----------------------------------------------------------------------
    
    function dydt = rhs(~,rho)
        
        
        flux = getFlux(rho); % -d_x \rho
        dydt =  flux; % d_xx \rho
        
        % Neumann problem

%         dydt(N+1)  = rho(N+1) - ExactSolution(x(N+1),t); % match to exact boundary condition
        dydt(1)  = rho(1) + 1;  % match to exact boundary condition
        
        
%         dydt(N+1)  = rho(N+1) - ExactSolution(x(N+1),t-delay); % match to exact boundary condition
%         dydt(1)  = rho(1) - ExactSolution(x(1),t-delay);  % match to exact boundary condition

        dydt = dydt(:);
                
    end

    function f = getFlux(rho)
        f = - C*Dy*rho;
    end

    function dydt = rhs2(~,rho)
        
        flux = getFlux2(rho);
        dydt = flux;
        
        dydt(1)= rho(1)+1;
        
        dydt = dydt(:);
        
    end

    function f= getFlux2(rho)
        f = -C*Dy2*rho;
    end    

    function [value, isterminal, direction] = event(~,y)
        Norm  = sobnorm(y, @tanmap, args);
        ratio = Norm/initNorm;
%         disp(ratio);
        
        value = (ratio>c1) & (ratio<c2);
        isterminal = 1;
        direction = 0;
    end


    function y = ExactSolution(x, t)
        y = InitialCondition(x-C*t);
    end

    
    
    
end