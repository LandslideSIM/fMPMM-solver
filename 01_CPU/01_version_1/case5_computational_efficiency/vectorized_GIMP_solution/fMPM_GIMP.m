% fMPMM-solver: Fast and efficient MATLAB-based MPM solver
%
% Copyright (C) 2020  Emmanuel Wyser, Michel Jaboyedoff, Yury Podladchikov
% -------------------------------------------------------------------------%
% version    : v1.0
% date       : february, 2020
% description: explicit mpm (GIMP) solver based on an updated Lagrangian
% frame in plane strain condition for elasto-plastic problem discretized on
% a 4-noded quadrilateral mesh
% -------------------------------------------------------------------------%
%% FANCY CALLS
clear                                                                     ;%
addpath('functions')                                                      ;%   
version    = 'GIMP_EP_'                                                    ;%
plasticity = true                                                         ;%         
typeD='double'                                                            ;%
disp('EXACT PLASTIC CORRECTION')                                          ;%
disp('------------------------')                                          ;%

set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;

numel = repmat(20,1,1)                                                    ;%
run   = zeros(length(numel),6)                                            ;%

numel = [20 40 80 160 320 640]+1
numel = 640;

meanDPIperRun = zeros(length(numel),7)                                    ;%

for sim=1:length(numel)
    disp('------------------------')                                      ;%
    disp(['Run ',num2str(sim),': nel = ',num2str(numel(sim)),''])         ;%
    disp('------------------------')                                      ;%
    
    %% NON-DIMENSIONAL CONSTANT
    ar      = 0.8                                                         ;% aspect ratio thickness/width
    nu      = 0.3                                                         ;% poisson ratio
    ni      = 3                                                           ;% number of mp in h(1) direction
    nstr    = 4                                                           ;% number of stress components
    %---------------------------------------------------------------------%
    
    %% INDEPENDANT PHYSICAL CONSTANT
    g       = 9.81                                                        ;% gravitationnal acceleration [m/s^2]
    Kc      = 0.7e6                                                       ;% bulk modulus                [Pa]
    E       = 3*Kc*(1-2*nu)                                               ;% young's modulus             [Pa]
    Gc      = E/(2*(1+nu))                                                ;% shear modulus               [Pa]
    rho0    = 2650                                                        ;% density                     [kg/m^3]
    n0      = 0.25                                                        ;% initial porosity            [-]
    yd      = sqrt(E/rho0)                                                ;% elastic wave velocity       [m/s]
    coh0    = 0.0e3                                                       ;% cohesion                    [Pa]
    phi0    = 19.8*pi/180                                                 ;% friction angle              [Rad]
    t       = 3.0                                                         ;% simulation time             [s]
    te      = 1.0                                                         ;% elastic loading             [s]
    %---------------------------------------------------------------------%
    
    %% MESH INITIALIZATION
    Lx    = 8.0                                                           ;%
    [meD,xn,yn] = meSetupCollapse(numel(sim),Lx,typeD)                    ;% - see function   
    % BOUNDARY CONDITIONS
    yinf    = 0.0                                                         ;%
    xB      = [-meD.L(1)/2 meD.L(1)/2]                                    ;%
    [row]   = find(meD.y<=yinf)                                           ;%
    BC.yi   = row                                                         ;%
    [row]   = find(meD.x<=xB(1)+0.25*meD.h(1))                            ;%
    BC.xi   = row                                                         ;%   
    [row]   = find(meD.x>=xB(2)-0.25*meD.h(1))                            ;%
    BC.xs   = row                                                         ;%   
    clear row                                                             ;%
    %---------------------------------------------------------------------%
    
    %% MPM DISCRETIZATION
    l0         = 1.0                                                      ;%
    [mpD]      = mpSetupCollapse(meD,ni,l0,xB(1),xB(2),yinf,coh0,phi0,n0,rho0,nstr,typeD);% - see function
    mpD.s(2,:) = -rho0.*g.*(l0-mpD.x(:,2));
    K0         = nu./(1-nu);
    mpD.s(1,:) = K0.*mpD.s(2,:);

    % figure(1),plot(meD.x,meD.y,'s',meD.x([BC.xi;BC.xs;BC.yi]),meD.y([BC.xi;BC.xs;BC.yi]),'gs',mpD.x(:,1),mpD.x(:,2),'x');axis equal;drawnow
    % ISOTROPIC ELASTIC MATRIX
    Del     = [ Kc+4/3*Gc,Kc-2/3*Gc,Kc-2/3*Gc,0.0;...
                Kc-2/3*Gc,Kc+4/3*Gc,Kc-2/3*Gc,0.0;...
                Kc-2/3*Gc,Kc-2/3*Gc,Kc+4/3*Gc,0.0;...
                0.0      ,0.0      ,0.0      ,Gc]                         ;%
    %---------------------------------------------------------------------%
    
    %% DISPLAY PARAMETERS AND RUNTIME INITIALIZATION
    fps  = 25                                                             ;% image per second
    %COURANT-FRIEDRICH-LEVY CONDITION
    dt   = 0.5*meD.h(1)/yd                                                ;% unconditionally stable timestep
    nit  = ceil(t/dt)                                                     ;% maximum number of interation
    nf   = max(2,ceil(round(1/dt)/fps))                                   ;% number of frame interval
    % GRAVITY INCREMENT
    dg   = linspace(0,g,round(((te))/dt))                                 ;% gravity increase
    %STORAGE DATA ALLOCATION
    runt = zeros(ceil(nit/nf),2)                                          ;% computational activity
    % RUNTIME PARAMETERS
    nc   = 0                                                              ;% initialize iteration counter                                                         
    it   = 1                                                              ;% initialize iteration
    tw   = 0.0                                                            ;% initialize time while statement
    cycle_time = zeros(nit,7)                                             ;%
    % PLOT SETTING
    [pp] = plotSetup(xB(1),xB(2),yinf,l0,meD.y,BC.yi)                     ;% - see function
    %---------------------------------------------------------------------%

    %% MPM MUSL VARIANT EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
    disp(['MPM SOLVER ON: ',num2str(meD.nN),' nodes, ',num2str(mpD.n),' material points']);
    tsolve=tic                                                            ;
    while(tw<t)% BEGIN WHILE LOOP
        time_it= tic                                                      ;%
        dpi    = time_it                                                  ;% CURRENT ITERATION TIMER BEGIN
%         if(it==1)
%         [vid1] = video('CPDI2q_slump_');
%         end
        %% LINEAR INCREASE OF GRAVITY FOR EQUILIBRIUM
        if(it<=size(dg,2))
            g = dg(end)                                                    ;% incremented gravity load
        end
        %------------------------------------------------------------------%
% %          %% ADAPTATIVE MP DOMAIN
% %         if(it*dt>te)
% %             l = 1.25;
% %             [mpD] = f_mpDomain(l,mpD);
% %             mpD.n = length(mpD.x);
% %         end
% %         %------------------------------------------------------------------%

        %% TRACK MATERIAL POINTS (p) IN ELEMENTS (E)
        p2.e = (floor((mpD.x(:,2)-min(meD.y))./meD.h(2))+1)+...
               (meD.nEy).*floor((mpD.x(:,1)-min(meD.x))./meD.h(1))        ;% index mp to element
        p2.N = meD.e2N(p2.e,:)                                            ;% index mp (to element) to node 
        l2g  = [meD.DoF*p2.N-1;meD.DoF*p2.N]                              ;% local to global node index list [x_I,y_I]
        neon = size(unique(p2.e),1)                                       ;% number of active element
        %------------------------------------------------------------------%
        
        %% BASIS FUNCTIONS
        tic;
        [mpD] = SdS(meD,mpD,p2.N)                                         ;% - see function 
        cycle_time(it,1)=toc;
        %------------------------------------------------------------------%
        %% MAPPING FROM MATERIAL POINTS (p) TO NODES (N)
        tic
        [meD] = p2Nsolve(meD,mpD,g,dt,l2g,p2.N,BC)                        ;% - see function
        cycle_time(it,2)=toc; 
        %------------------------------------------------------------------%

        %% MAPPING FROM NODES (N) TO MATERIAL POINTS (p)
        tic
        [meD,mpD] = mapN2p(meD,mpD,dt,l2g,p2.N,BC)                        ;% - see function
        cycle_time(it,3)=toc;
        %------------------------------------------------------------------%
        
        %% UPDATE INCREMENTAL DEFORMATION & STRAIN 
        tic
        [mpD] = DefUpdate(meD,mpD,l2g)                                    ;% - see function
        cycle_time(it,4)=toc;
        %------------------------------------------------------------------%   
        
        %% ELASTO-PLASTIC RELATION: ELASTIC PREDICTOR - PLASTIC CORRECTOR
        tic
        [mpD] = constitutive(mpD,Del,Gc,Kc,plasticity,dt,it,te)           ;%
        cycle_time(it,5)=toc                                              ;%
        dpi=toc(dpi)                                                      ;% CURRENT ITERATION TIMER END
        cycle_time(it,6)=dpi                                              ;%
        cycle_time(it,7)=length(unique(p2.N))                             ;%
        %------------------------------------------------------------------%
        
        %% TERMINAL DISPLAY
        if(mod(it,nf)==1)
            nc    = nc+1                                                  ;%
            clocktimer(((nit-it)*toc(time_it)),'Remaining estimated time:',1/dpi);%
            runt(nc,:) = [1/dpi mean(runt(1:nc,1))]                       ;%
%             if(dt*it>te)
%                 figure(2)
%                 clf
%                 mpD.P    = (mpD.s'*[1;1;1;0])./3;
%                 pp.cbchar='$\epsilon_{\mathrm{eqv}}^p$';
%                 pp.cbpos =[0.42 0.5 0.2 0.05];
%                 pos            = pp.cbpos;
%                 pp.cblpos=[pos(1) pos(2)+2];
%                 %plotprop.caxis =(0.25*E)/1e3;
%                 pp.caxis =8;
%                 pp.ps    =3.0;
%                 pp.tit   =['time: ',num2str(it*dt-te,'%.2f'),' (s), $g=',num2str(g,'%.2f'),'$ (m/s$^2$)'];
%                 pp.cbclass = 40;
%                 dd = mpD.epII;
%                 dis(dd,mpD.x(:,1),mpD.x(:,2),it*dt,pp);
%                 set(gca,'xtick',[xB(1):0.125*(xB(2)-xB(1)):xB(2)]);
%                 grid on;
% %                writeVideo(vid1,getframe(gcf));
%             end
        end
        %------------------------------------------------------------------%        
        
        %% CONDITION: NO ACTIVE ELEMENT
        if(isempty(neon)==1 || sum(isnan(mpD.v(:,1))+isnan(mpD.v(:,2)))>0.0)
            break                                                         ;%
        end
        %------------------------------------------------------------------%
        
        %% ITERATION INCREMENT
        tw=tw+dt                                                          ;%
        it=it+1                                                           ;%
        %------------------------------------------------------------------%
        
    end% END WHILE LOOP
    tsolve = toc(tsolve)                                                  ;%
    clocktimer(tsolve,'Runtime MPM solver:',mean(runt(:,2)))              ;%

    fig2=figure(54363);
    clf
    set(fig2,'Units','pixels','Position',[200 200 364 353]);
    xs=mpD.x(:,1)+[-1 -1 1 1 -1].*mpD.l(:,1);
    ys=mpD.x(:,2)+[-1  1 1 -1 -1].*mpD.l(:,2);
    hold on
    ax1=plot(xs',ys','b-',mpD.x(:,1),mpD.x(:,2),'b.');axis equal;axis tight;
    xN = reshape(meD.x,meD.nNy,meD.nNx);
    yN = reshape(meD.y,meD.nNy,meD.nNx);
    hold off
    xlim([xB(1) xB(1)+5]);
    ylim([0 5])
    xlabel('$x$ (m)');ylabel('$y$ (m)');
    box on;
    title(['$t=',num2str(dt*nit,'%.2f'),'$ (s)'])
    set(gca,'FontSize',15,'TickLabelInterpreter','latex');
    axis tight
    
    name=['.\data\GIMP_time_vectorized_' num2str(sim) '.mat'];
    save(name,'runt','cycle_time','tsolve');
    name=['.\data\GIMP_vectorized_data_' num2str(sim) '.mat'];
    save(name,'mpD','meD','ni','rho0','E','nu','dt','nit');
end


