% fMPMM-solver: Fast and efficient MATLAB-based MPM solver
%
% Copyright (C) 2020  Emmanuel Wyser, Yury Alkhimenkov, Michel Jaboyedoff, Yury Podladchikov
% -------------------------------------------------------------------------%
% version    : v1.1
% date       : october, 2020
% description: explicit mpm solver based on an updated Lagrangian
% frame in plane strain condition for elasto-plastic problem discretized on
% a 4-noded quadrilateral mesh
%% FANCY CALLS
clear                                                                     ;%
addpath('functions')                                                      ;%   
version    = 'CPDI_EP_'                                                   ;%
plasticity = true                                                         ;%         
typeD='double'                                                            ;%
disp('EXACT PLASTIC CORRECTION')                                          ;%
disp('------------------------')                                          ;%

set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;

numel = repmat(110,1,1)                                                 ;%                                          ;%
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
    E       = 70.0e6                                                       ;% young's modulus             [Pa]
    Gc      = E/(2*(1+nu))                                                ;% shear modulus               [Pa]
    Kc      = E/(3*(1-2*nu))                                              ;% bulk modulus                [Pa]
    rho0    = 2100                                                        ;% density                     [kg/m^3]
    n0      = 0.25                                                        ;% initial porosity            [-]
    yd      = sqrt((Kc+4/3*Gc)/rho0)                                      ;% elastic wave speed          [m/s]
    coh0    = 10.0e3                                                      ;% cohesion                    [Pa]
    phi0    = 20.0*pi/180                                                 ;% friction angle              [Rad]
    H       = -60e3                                                       ;% softening modulus           [Pa]
    cohr    =  4.0e3                                                      ;% residual cohesion           [Pa]
    phir    = 7.5*pi/180                                                  ;% residual friction angle     [Rad]
    t       = 7.2                                                         ;% simulation time             [s]
    te      = 0.0                                                         ;% elastic loading             [s]
    %---------------------------------------------------------------------%
    
    %% MESH & MP INITIALIZATION
    [meD,bc] = meSetup(numel(sim),typeD)                                  ;% - see function   
    ly      = 35                                                       ;% layer thickness [m]
    [mpD]   = mpSetup(meD,ni,ly,coh0,phi0,n0,rho0,nstr,typeD)   ;% - see function
    % ISOTROPIC ELASTIC MATRIX
    Del     = [ Kc+4/3*Gc,Kc-2/3*Gc,Kc-2/3*Gc,0.0;...
                Kc-2/3*Gc,Kc+4/3*Gc,Kc-2/3*Gc,0.0;...
                Kc-2/3*Gc,Kc-2/3*Gc,Kc+4/3*Gc,0.0;...
                0.0      ,0.0      ,0.0      ,Gc]                         ;%
    Hp      = H*meD.h(1)                                                  ;%
    %---------------------------------------------------------------------%
    meD.vd = 0.0;
    %% DISPLAY PARAMETERS AND RUNTIME INITIALIZATION
    fps  = 25                                                             ;% image per second
    %COURANT-FRIEDRICH-LEVY CONDITION
    C    = 0.5                                                            ;%                 
    dt   = C*meD.h(1)/yd                                                  ;% unconditionally stable timestep
    nit  = ceil(t/dt)                                                     ;% maximum number of interation
    nf   = max(2,ceil(round(1/dt)/fps))                                   ;% number of frame interval
    % RUNTIME PARAMETERS
    nc   = 0                                                              ;% initialize iteration counter                                                         
    it   = 1                                                              ;% initialize iteration
    tw   = 0.0                                                            ;% initialize time while statement
    cycle_time = zeros(nit,7)                                             ;%
    % PLOT SETTING
    [pp] = plotSetup(meD.xB(1),meD.xB(2),0.0,ly,meD.y,bc.y)               ;% - see function

    %---------------------------------------------------------------------%
    %% MPM DM ALGORITHM EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
    tsolve=tic                                                            ;%
    fprintf('MPM SOLVER ON: %.0f elements \n'      ,meD.nEx*meD.nEy)      ;%
    fprintf('               %.0f nodes \n'         ,meD.nN)               ;%
    fprintf('               %.0f material point \n',mpD.n)                ;%    
    while((tw<t)||(sum(isnan(mpD.v(:,1))+isnan(mpD.v(:,2)))>0.0))          % BEGIN WHILE LOOP
        time_it= tic                                                      ;%
        dpi    = time_it                                                  ;% CURRENT ITERATION TIMER BEGIN
     
        %% CFL CONDITION & LINEAR INCREASE OF GRAVITY FOR EQUILIBRIUM
        c  = [max(yd+abs(mpD.v(:,1))),max(yd+abs(mpD.v(:,2)))]            ;%
        dt = C*min(meD.h./c)                                              ;%
        DT(it) = dt;
        %------------------------------------------------------------------%
        
        %% TRACK MATERIAL POINTS CORNER (C) IN ELEMENTS (E)
        c2e  = reshape((permute((floor((mpD.yc-min(meD.y))./meD.h(2))+1)+...
               (meD.nEy).*floor((mpD.xc-min(meD.x))./meD.h(1)),[2 1])),mpD.n*4,1);%
        neon = length(unique(c2e));
        c2N  = reshape((meD.e2N(c2e,:)'),meD.nNp,mpD.n)';
        l2g  = [meD.DoF*c2N-1;meD.DoF*c2N]                                ;% local to global node index list [x_I,y_I]
        %------------------------------------------------------------------%
        
        %% BASIS FUNCTIONS
        tic;
        [mpD] = SdS(mpD,meD,c2N)                                    ;% - see function 
        cycle_time(it,1)=toc;
        %------------------------------------------------------------------%
        %% PROJECTION FROM MATERIAL POINTS (p) TO NODES (N)
        tic
        [meD] = p2Nsolve(meD,mpD,g,dt,l2g,c2N,bc)                           ;% - see function
        cycle_time(it,2)=toc; 
        %------------------------------------------------------------------%

        %% INTERPOLATION FROM NODAL SOLUTIONS (N) TO MATERIAL POINTS (p)
        tic
        [meD,mpD] = mapN2p(meD,mpD,dt,l2g,c2N,bc)                           ;% - see function
        cycle_time(it,3)=toc;
        %------------------------------------------------------------------%
        
        %% UPDATE INCREMENTAL DEFORMATION & STRAIN 
        tic
        [mpD] = DefUpdate(meD,mpD,l2g)                                  ;% - see function
        cycle_time(it,4)=toc;
        %------------------------------------------------------------------%   
        
        %% ELASTO-PLASTIC RELATION: ELASTIC PREDICTOR - PLASTIC CORRECTOR
        tic
        [mpD] = constitutive(mpD,Del,plasticity,dt,it,te)         ;%
        cycle_time(it,5)=toc                                              ;%
        dpi=toc(dpi)                                                      ;% CURRENT ITERATION TIMER END
        cycle_time(it,6)=dpi                                              ;%
        %------------------------------------------------------------------%
        
        %% TERMINAL DISPLAY
        if(mod(it,nf)==1)
            rt    = ((nit-it)*toc(time_it))                               ;%
            dpi   = mean(1./cycle_time(1:it,6))                           ;%
            clocktimer(rt,'Remaining estimated time:',dpi)                ;%
        end
        if((tw<te)&&(tw>=te-dt))
            figure(2)
            pp.cbchar='$p$ [kPa]';
            pp.cbpos =[0.42 0.5 0.2 0.05];
            pos            = pp.cbpos;
            pp.cblpos=[pos(1) pos(2)+2];
            pp.caxis =(rho0*g*max(mpD.x(:,2)))/1e3;
            pp.ps    =3.0;
            pp.tit   =['time: ',num2str(it*dt-te,'%.2f'),' (s), $g=',num2str(g,'%.2f'),'$ (m/s$^2$)'];
            pp.cbclass = 20;
            dd = -((mpD.s(1,:)+mpD.s(2,:)+mpD.s(3,:))./3)./1e3;
            dis(dd,mpD.x(:,1),mpD.x(:,2),it*dt,pp);
            print(figure(2),['./data/' version 'pressure_equilibrium' '_',num2str(sim),'' ],'-dpng');
            mpD.u(:,:) = 0.0;
        end
        %------------------------------------------------------------------%        
        
        %% ITERATION INCREMENT
        tw=tw+dt                                                          ;%
        it=it+1                                                           ;%
        %------------------------------------------------------------------%
        
    end% END WHILE LOOP
    tsolve = toc(tsolve)                                                  ;%
    clocktimer(tsolve,'Runtime MPM solver:',mean(1./cycle_time(:,6)))     ;%
    %----------------------------------------------------------------------%
    
% %     %% SAVE DATA
% %     name=['.\data\CPDI_time_vectorized_' num2str(sim) '.mat'];
% %     save(name,'cycle_time','tsolve');
% %     name=['.\data\data_' num2str(sim) '.mat'];
% %     save(name,'mpD','meD','ni','rho0','E','nu','dt','nit','tsolve');
% %     %----------------------------------------------------------------------%
    
end

fig3=figure(3);
set(fig3,'Units','pixels','Position',[139 98.3333 504 458.6667]);
subplot(311)
xs=[mpD.xc mpD.xc(:,1)];
ys=[mpD.yc mpD.yc(:,1)];
plot(xs',ys','k-');axis equal;axis tight;
ylabel('$y$ (m)');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
xlim([0 meD.L(1)]);
ylim([0 ly]);
title('Finite deformation')
subplot(312)
ax1=scatter(mpD.x(:,1),mpD.x(:,2),pp.ps,mpD.epII,'filled');
axis equal;
box on;
ylabel('$y$ (m)');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
colormap(gca,jet);
xlim([0 meD.L(1)]);
ylim([0 ly]);
title('$\epsilon_{II}$')
title(['$\epsilon_{II}$, max($\epsilon_{II}$) = ',num2str(max(mpD.epII),'%.2f'),' [-]']);
subplot(313)
du = sqrt(mpD.u(:,1).^2+mpD.u(:,2).^2);
ax2=scatter(mpD.x(:,1),mpD.x(:,2),pp.ps,du,'filled');
axis equal;
box on;
xlabel('$x$ (m)');ylabel('$y$ (m)');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
colormap(gca,(jet));
xlim([0 meD.L(1)]);
ylim([0 ly]);
title(['$\Delta u$, max($\Delta u$) = ',num2str(max(du),'%.2f'),' [m]']);
print(fig3,['./data/' version 'summary' '_',num2str(sim),'' ],'-dpng');
path = ['C:\Users\agnes\Desktop\MPM\work\stable\variant_Difference\' version 'summary' '_',num2str(sim),''];
print(fig3,path,'-dpng');

figure(4)
X = repmat(1:length(cycle_time),5,1)';
Y = cycle_time(:,1:5);
plot(X,Y);
legend('basis function','p2n','n2p','deformation','consitutive')



  free_surf=load('C:\Users\agnes\Desktop\MPM\work\technical_note\case5_ElastoPlasticSlump_new\Huang_etal_Solution\free_surface.txt');
    fail_surf=load('C:\Users\agnes\Desktop\MPM\work\technical_note\case5_ElastoPlasticSlump_new\Huang_etal_Solution\failure_surface.txt');


    fig1=figure(1)
    set(fig1,'Units','pixels','Position',[100 287.6667 541 241.3333],'Color','white');
    pp.cbchar='$\log_{10}(\epsilon_{\mathrm{II}})$';
    pp.cbpos =[0.42 0.5 0.2 0.05];
    pos            = pp.cbpos;
    pp.cblpos=[pos(1)-(pos(3)) pos(2)+2];
    pp.caxis =log10(max(mpD.epII));
    pp.XTickLabel = {'-0.5','0.5'}
    pp.XTick = {'-0.5','0.5'}
    %pp.tit   =['time: ',num2str(it*dt-te,'%.2f'),' (s), $e_{on}=',num2str(neon),'$, $n_p \in e_i = ',num2str(ni^2),'$ '];
    pp.tit   =['$t=',num2str(it*dt-te,'%.2f'),'$ (s)'];
    D = log10(mpD.epII);
    dis(D,mpD.x(:,1),mpD.x(:,2),it*dt,pp);
    hold on
    ax1=plot(free_surf(:,1),free_surf(:,2),'r-','LineWidth',2.5)
    ax2=plot(fail_surf(:,1),fail_surf(:,2),'r:','LineWidth',2.5)
    hold off
    ylim([0 40])
    tit = {'D-P model: final geometry','D-P model: failure surface'};
    h1=legend([ax1 ax2],tit);
    set(h1,'Interpreter','latex','FontSize',12,'Position',[0.4561 0.6718 0.4413 0.1778],'NumColumns',1);
    set(fig1, 'InvertHardCopy', 'off');
    print(fig1,['./data/' version 'vectorized_plastic_strain_invariant' '_',num2str(sim),'' ],'-dpng');
    print(fig1,['./data/' version 'vectorized_plastic_strain_invariant' '_',num2str(sim),'' ],'-depsc');

%% REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
% Dunatunga, S., Kamrin, K. (2017) Continuum modeling of projectile impact
% and penetration in dry granular media. J. Mech. Phys. Solids. 100. 45-60
%-------------------------------------------------------------------------%
% Dunatunga, S., Kamrin, K. (2015) Continuum modelling and simulation of
% granular flows through their many phases. J. Fluid Mech. 779. 483-513.
%-------------------------------------------------------------------------%
% Stomakhin, A., Schroeder, C., Chai, L., Teran, C., Selle, A. (2013) A
% material point method for snow simulation. ACM Trans. Graph. 32.
%-------------------------------------------------------------------------%
% Gaume, J., Gast, T., Teran, J., van Herwijnen, A., Jiang, C. (2018)
% Dynamic anticrack propagation in snow. Nat. Com. 9
%-------------------------------------------------------------------------%
% Wang, B., Vardon, P.J., Hicks, M.A. (2016) Investigation of retrogressive
% and progressive slope failure mechanisms using the material point method.
% Comp. and Geot. 78. 88-98
%-------------------------------------------------------------------------%
% Liang, W., Zhao, J. (2018) Multiscale modeling of large deformation in
% geomechanics. Int. J. Anal. Methods Geomech. 43. 1080-1114
%-------------------------------------------------------------------------%
% Baumgarten, A.S., Kamrin, K. (2018) A general fluid-sediment mixture
% model and constitutive theory validated in many flow regimes. J. Fluid
% Mech. 861. 7211-764.
%-------------------------------------------------------------------------%
% Huang, P., Zhang, X., Ma, S., Huang, X. (2011) Contact algorithms for the
% material point method in impact and penetration simulation. Int. J. Numer.
% Meth. Engng. 85. 498-517.
%-------------------------------------------------------------------------%
% Homel, M., Herbold, E.B. (2017) Field-gradient partitioning for fracture
% and frictional contact in the material point method. Int. J. Numer. Meth.
% Engng. 109. 1013-1044
%-------------------------------------------------------------------------%
% Ma. S., Zhank, X., Qui, X.M. (2009) Comparison study of MPM and SPH in
% modeling hypervelocity impact problems. Int. J. Imp. Engng. 36. 272-282.
%-------------------------------------------------------------------------%
% Bardenhagen, S.G., Brackbill, J.U., Sulsky, D. (2000) The material-point
% method for granular materials. Comput. Methods Appl. Mech. Engrg. 187.
% 529-541.
%-------------------------------------------------------------------------%
% York, A.R., Sulsky, D., Schreyer, H.L. (1999) The material point method
% for simulation of thin membranes. Int. J. Numer. Meth. Engng. 44. 1429-1456
%-------------------------------------------------------------------------%
% Nairn, J.A., Bardenhagen, S.G., Smith, G.D. (2018) Generalized contact
% and improved frictional heating in the material point method. Comp. Part.
% Mech. 3. 285-296
%-------------------------------------------------------------------------%
% Nairn, J.A. (2013) Modeling imperfect interfaces in the material point
% method using multimaterial methods. Comp. Mod. Engrg. Sci. 92. 271-299.
%-------------------------------------------------------------------------%
% Hammerquist, C.C., Nairn, J.A. (2018) Modeling nanoindentation using the
% material point method. J. Mater. Res. 33. 1369-1381
%-------------------------------------------------------------------------%
% Hamad, F., Stolle, D., Moormann, C. (2016) Material point modelling of
% releasing geocontainers from a barge. Geotext. Geomembr. 44. 308-318.
%-------------------------------------------------------------------------%
% Bhandari, T., Hamad, F., Moormann, C., Sharma, K.G., Westrich, B. (2016)
% Comp. Geotech. 75. 126-134.
%-------------------------------------------------------------------------%
% Wang, B., Vardon, P.J., Hicks, M.A. (2016) Investigation of retrogressive
% and progressive slope failure mechanisms using the material point method.
% Comp. Geotech. 78. 88-98.
%-------------------------------------------------------------------------%
% Keller, T., May, D.A., Kaus, B.J.P. (2013) Numerical modelling of magma
% dynamics oupled to tectonic deformation of lithosphere and crust.
% Geophys. J. Int. 195. 1406-1442.