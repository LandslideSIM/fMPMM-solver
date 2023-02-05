% fMPMM-solver: Fast and efficient MATLAB-based MPM solver
%
% Copyright (C) 2020  Emmanuel Wyser, Michel Jaboyedoff, Yury Podladchikov
% -------------------------------------------------------------------------%
% version    : v1.0
% date       : february, 2020
% description: explicit mpm (CPDI2q) solver based on an updated Lagrangian
% frame in plane strain condition for elasto-plastic problem discretized on
% a 4-noded quadrilateral mesh
% -------------------------------------------------------------------------%
%% FANCY CALLS
clear                                                                     ;%
addpath('functions')                                                      ;%
version    = 'CPDI_E_'                                                   ;%
plasticity = false                                                         ;%
typeD='double'                                                            ;%
disp('EXACT PLASTIC CORRECTION')                                          ;%
disp('------------------------')                                          ;%

set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;

numel = repmat(20,1,1)                                                    ;%
run   = zeros(length(numel),6)                                            ;%

for sim=1:length(numel)
    disp('------------------------')                                      ;%
    disp(['Run ',num2str(sim),': nel = ',num2str(numel(sim)),''])         ;%
    disp('------------------------')                                      ;%
    
    %% NON-DIMENSIONAL CONSTANT
    ar      = 0.8                                                         ;% aspect ratio thickness/width
    nu      = 0.3                                                         ;% poisson ratio
    ni      = 4                                                           ;% number of mp in h(1) direction
    nstr    = 4                                                           ;% number of stress components
    %---------------------------------------------------------------------%
    
    %% INDEPENDANT PHYSICAL CONSTANT
    g       = 0.0                                                         ;% gravitationnal acceleration [m/s^2]
    E       = 1.0e4                                                       ;% young's modulus             [Pa]
    Gc      = E/(2*(1+nu))                                                ;% shear modulus               [Pa]
    Kc      = E/(3*(1-2*nu))                                              ;% bulk modulus                [Pa]
    rho0    = 1000                                                        ;% density                     [kg/m^3]
    n0      = 0.25                                                        ;% initial porosity            [-]
    yd      = sqrt(E/rho0)                                                ;% elastic wave velocity       [m/s]
    coh0    = 20.0e3                                                      ;% cohesion                    [Pa]
    phi0    = 30.0*pi/180                                                 ;% friction angle              [Rad]
    H       = -60e3                                                       ;% softening modulus           [Pa]
    cohr    =  4.0e3                                                      ;% residual cohesion           [Pa]
    phir    = 7.5*pi/180                                                  ;% residual friction angle     [Rad]
    t       = 1.13                                                        ;% simulation time             [s]
    te      = 0.0                                                         ;% elastic loading             [s]
    %---------------------------------------------------------------------%
    
    %% MESH INITIALIZATION
    L     = 1;
    [meD] = meSetup(numel(sim),L,typeD)                                     ;% - see function
    % BOUNDARY CONDITIONS
    yinf    = 0.0                                                         ;%
    xB      = [-meD.L(1)/2 meD.L(1)/2]                                    ;%
    [row]   = find(meD.y<=yinf)                                           ;%
    BC.yi   = []                                                          ;%
    [row]   = find(meD.x<=xB(1))                                          ;%
    BC.xi   = []                                                         ;%
    [row]   = find(meD.x>=xB(2))                                          ;%
    BC.xs   = []                                                         ;%
    clear row                                                             ;%
    %---------------------------------------------------------------------%
    
    %% MPM DISCRETIZATION
    lx      = (ar*meD.L(1))                                               ;% layer length [m]
    [mpD]   = mpSetup(meD,ni,L,xB(1),xB(2),yinf,coh0,cohr,phi0,phir,n0,rho0,nstr,typeD);% - see function
    % ISOTROPIC ELASTIC MATRIX
    Del     = [ Kc+4/3*Gc,Kc-2/3*Gc,Kc-2/3*Gc,0.0;...
        Kc-2/3*Gc,Kc+4/3*Gc,Kc-2/3*Gc,0.0;...
        Kc-2/3*Gc,Kc-2/3*Gc,Kc+4/3*Gc,0.0;...
        0.0      ,0.0      ,0.0      ,Gc]                         ;%
    Hp      = H*meD.h(1)                                                  ;%
    %---------------------------------------------------------------------%
    
    %% DISPLAY PARAMETERS AND RUNTIME INITIALIZATION
    fps  = 100                                                             ;% image per second
    %COURANT-FRIEDRICH-LEVY CONDITION
    dt   = 0.1*meD.h(1)/yd                                                ;% unconditionally stable timestep
    nit  = ceil(t/dt)                                                     ;% maximum number of interation
    nf   = max(2,ceil(round(1/dt)/fps))                                   ;% number of frame interval
    % GRAVITY INCREMENT
    dg   = linspace(0,g,round(((te)/1.5)/dt))                             ;% gravity increase
    %STORAGE DATA ALLOCATION
    runt = zeros(ceil(nit/nf)+1,2)                                         ;% computational activity
    % RUNTIME PARAMETERS
    nc   = 0                                                              ;% initialize iteration counter
    it   = 1                                                              ;% initialize iteration
    tw   = 0.0                                                            ;% initialize time while statement
    cycle_time = zeros(nit,7)                                             ;%
    % PLOT SETTING
    %[pp] = plotSetup(xB(1),xB(2),yinf,ly,meD.y,BC.yi)                     ;% - see function
    %---------------------------------------------------------------------%
    
    %% MPM MUSL VARIANT EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
    disp(['MPM SOLVER ON: ',num2str(meD.nN),' nodes, ',num2str(mpD.n),' material points']);
    tsolve=tic                                                            ;
    while(tw<t)% BEGIN WHILE LOOP
        time_it= tic                                                      ;%
        dpi    = time_it                                                  ;% CURRENT ITERATION TIMER BEGIN
        %------------------------------------------------------------------%
        
        %% TRACK MATERIAL POINTS CORNER (C) IN ELEMENTS (E)
        c2e  = reshape((permute((floor((mpD.yc-min(meD.y))./meD.h(2))+1)+...
            (meD.nEy).*floor((mpD.xc-min(meD.x))./meD.h(1)),[2 1])),mpD.n*4,1);%
        neon = length(unique(c2e));
        c2N  = reshape((meD.e2N(c2e,:)'),meD.nNp,mpD.n)';
        %------------------------------------------------------------------%
        
        %% BASIS FUNCTIONS
        tic;
        [mpD,mp,no] = SdS(mpD,meD,c2N)                                    ;% - see function
        cycle_time(it,1)=toc;
        %------------------------------------------------------------------%
        %% PROJECTION FROM MATERIAL POINTS (p) TO NODES (N)
        tic
        [meD] = p2Nsolve(meD,mpD,g,dt,mp,no,BC)                           ;% - see function
        cycle_time(it,2)=toc;
        %------------------------------------------------------------------%
        
        %% INTERPOLATION FROM NODAL SOLUTIONS (N) TO MATERIAL POINTS (p)
        tic
        [meD,mpD] = mapN2p(meD,mpD,dt,mp,no,BC)                           ;% - see function
        cycle_time(it,3)=toc;
        %------------------------------------------------------------------%
        
        %% UPDATE INCREMENTAL DEFORMATION & STRAIN
        tic
        [mpD] = DefUpdate(meD,mpD,mp,no)                                  ;% - see function
        cycle_time(it,4)=toc;
        %------------------------------------------------------------------%
        
        %% ELASTO-PLASTIC RELATION: ELASTIC PREDICTOR - PLASTIC CORRECTOR
        tic
        [mpD] = constitutive(mpD,Del,Hp,cohr,plasticity,dt,it,te)         ;%
        cycle_time(it,5)=toc                                              ;%
        dpi=toc(dpi)                                                      ;% CURRENT ITERATION TIMER END
        cycle_time(it,6)=dpi                                              ;%
        %------------------------------------------------------------------%
        
        %% TERMINAL DISPLAY
        if(mod(it,nf)==1)
            nc    = nc+1                                                  ;%
            clocktimer(((nit-it)*toc(time_it)),'Remaining estimated time:',1/dpi);%
            runt(nc,:) = [1/dpi mean(runt(1:nc,1))]                       ;%
            
%             fig2=figure(54363);
%             clf
%             set(fig2,'Units','pixels','Position',[69.6667 41.6667 621.3333 599.3333]);
%             plot(mpD.x(:,1),mpD.x(:,2),'o');
%             axis equal;axis tight;
%             xlabel('$x$ (m)');ylabel('$y$ (m)');
%             xlim([-L/2 L/2])
%             ylim([0 L])
%             set(gca,'FontSize',15,'TickLabelInterpreter','latex');
        end

        %------------------------------------------------------------------%
        
        %% CONDITION: NO ACTIVE ELEMENT
        if(sum(isnan(mpD.v(:,1))+isnan(mpD.v(:,2)))>0.0)
            break                                                         ;%
        end
        %------------------------------------------------------------------%
 
        K(it)  = 0.5.*sum(sqrt(mpD.v(:,1).^2+mpD.v(:,2).^2).*mpD.m,1);
        kappa  = 3-4*nu;

        dsigma = (((kappa+1)/4).*(mpD.s(1,:).^2+mpD.s(2,:).^2)-2.*(mpD.s(1,:).*mpD.s(2,:)-mpD.s(4,:).^2));
        U(it)  = sum(mpD.V'.*(1./(4.*Gc)).*dsigma,2);

        %% ITERATION INCREMENT
        tw=tw+dt                                                          ;%
        it=it+1                                                           ;%
        %------------------------------------------------------------------%
        
    end% END WHILE LOOP
    tsolve = toc(tsolve)                                                  ;%
    clocktimer(tsolve,'Runtime MPM solver:',mean(1./cycle_time(it,6)))    ;%

    fig2=figure(54363);
    clf
    set(fig2,'Units','pixels','Position',[69.6667 41.6667 621.3333 599.3333]);
    xs=[mpD.xc mpD.xc(:,1)];
    ys=[mpD.yc mpD.yc(:,1)];
    plot(xs',ys','k-');axis equal;axis tight;
    xlabel('$x$ (m)');ylabel('$y$ (m)');
    set(gca,'FontSize',15,'TickLabelInterpreter','latex');
    
    figure(2)
    plot(cycle_time(:,1:5))
    % %     fig2=figure(54363);
    % %     clf
    % %     set(fig2,'Units','pixels','Position',[1882 647 1127 283]);
    % %     xs=mpD.x(:,1)+[-1 -1 1 1 -1].*mpD.l(:,1);
    % %     ys=mpD.x(:,2)+[-1  1 1 -1 -1].*mpD.l(:,2);
    % %     plot(xs',ys','k-');axis equal;axis tight;
    % %     xlabel('$x$ (m)');ylabel('$y$ (m)');
    % %     set(gca,'FontSize',15,'TickLabelInterpreter','latex');
    % %     print(fig2,['./data/' version '_fem_exact_domain' '_',num2str(sim),'' ],'-dpng','-r0');
    % %
    %     figure(3)
    %     plot(1:length(runt),runt(:,1),'r--',1:length(runt),runt(:,2),'k-')
    
    name=['.\data\CPDI_time_vectorized_' num2str(sim) '.mat'];
    save(name,'runt','cycle_time');
end
