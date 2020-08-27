function GeneralPDE()

    N = 20;
    y = chebx(N);
%    disp(y)
    Dy = diffm(N, @linmap, [-1,1]);


    function rho0 = InitialCondition(y)

%            rho0 = exp(-3*(y+0.3).^2);
%             for i = 1:length(y)
%                 if i>length(y)/2
%                     rho0(i) = 1;
%                 else
%                     rho0(i)= 0;
%                 end
%             end
            rho0 = tanh(20*(y+0.6)) + tanh(-20*(y-0.3));
            rho0 = rho0;
    end

    rho_ic = InitialCondition(y);
     rho_ic(1) = 0;
%     rho_ic(N) = 0;


    tMax = 4;
    outtimes = [0:tMax/100:tMax];

    mM = ones(size(y));
      mM(1) = 0;
%      mM(N+1) = 0;
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    
    tic
    [~,Rho_t] = ode15s(@rhs,outtimes,rho_ic,opts);
    toc
    
    %figure('Position', [100,100,1200,700]);
    
    %m = zeros(length(outtimes),1);
    xx = linspace(-1,1,500)';
    IntM = baryM(y, xx);
    for iTimes = 1:length(outtimes)
        rhot = Rho_t(iTimes,:)';
        
        hold off
      
        %plot(y,rhot);
        plot(xx, IntM*rhot, y, rhot,'.');
         ylim([-2.1,2.1]);
%         xlim([-1,1]);
        
        pause(0.03)
        
    end
    
    function dydt = rhs(t,rho)
        
        flux =   -0.8*Dy*rho; % -d_x \rho
        dydt =   flux; % d_xx \rho
        
        % Neumann problem
%         dydt(N+1)  = flux(N+1) - 0; % make d_x \rho = 0 on right endpoint
%          dydt(1)   = flux(1) - 0;  % make d_x \rho = 0 on left endpoint


        % Dirichlet problem
%          dydt(N+1)  =  rho(N+1) - 0; % make \rho = 0 on right endpoint
         dydt(1)   =  rho(1) + 0;  % make \rho = 0 on left endpoint
% %         
        dydt = dydt(:);
                
    end
    
    function flux = getFlux(rho)
        
        flux = -Dy*rho;
    end
end
