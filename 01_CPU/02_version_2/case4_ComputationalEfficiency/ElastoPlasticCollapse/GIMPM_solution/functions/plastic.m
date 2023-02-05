function [mpD,D] = plastic(mpD,Hp,cohr,Del)
%UNTITLED Summary of this function goes here
if(mpD.plast == 1)
    Q= [ 2/3 -1/3 -1/3 0;...
        -1/3  2/3 -1/3 0;...
        -1/3 -1/3  2/3 0;...
        0    0    0   2];
    corr         = 0;
    ftol         = 1e-9;
    % yield function
    coh          = mpD.coh+Hp.*mpD.epII;
    coh(coh<cohr)= cohr;
    h            = -sqrt(3);
    sig0         = sqrt(3).*coh;
    f            = sqrt(3./2.*sum((mpD.s'*Q).*mpD.s',2))-sig0';
    %-------------------------------------------------------------%
    F = f(f>0.0);
    [yld]=find(f>0.0);
    while(max(F)>ftol)% Radial Return Algorithm
        %-----------------------------------------------------%
        dfds = (3./(2.*repmat(sig0(yld),4,1))).*(Q*mpD.s(:,yld)); % De Borst etal, 2012, p.243
        dlam = F./(h+sum(dfds'*Del.*dfds',2));
        dsig = repmat(dlam',4,1).*(Del*dfds);
        dep  = Del\dsig;
        %-----------------------------------------------------%
        % correct trial stress
        mpD.s(:,yld)  = mpD.s(:,yld) - dsig;
        %-----------------------------------------------------%
        % compute yield function
        mpD.epII(yld) = mpD.epII(yld)+sqrt(2/3.*(dep(1,:).^2+dep(2,:).^2+2.*dep(4,:).^2));
        coh          = mpD.coh+Hp.*mpD.epII;
        coh(coh<cohr)= cohr;
        sig0         = sqrt(3).*coh;
        F            = sqrt(3./2.*sum((mpD.s(:,yld)'*Q).*mpD.s(:,yld)',2))-sig0(yld)';
        %-----------------------------------------------------%
        % save iteration
        corr = corr+1;
        err(corr) = corr;
        if(corr>50)
            break
        end
    end
elseif(mpD.plast == 2)
        c            = mpD.coh+(Hp*mpD.epII ); c(c<cohr)=cohr             ;% hardening if Hp \neq 0 otherwise perfectly plastic material
        %% MOHR-COULOMB YIELD FUNCTION
        ds           = mpD.s(1,:)-mpD.s(2,:)                              ;% preprocessing: sig_xx - sig_yy
        tau          = sqrt(0.25*(ds).^2 + mpD.s(4,:).^2)                 ;% deviatoric stress
        sig          = 0.5*(mpD.s(1,:)+mpD.s(2,:))                        ;% mohr-circle center: pressure
        f            = tau+sig.*sin(mpD.phi)-c.*cos(mpD.phi)              ;% yield function f<=0, if f>0 then return mapping
        mpD.sn       = mpD.s                                              ;% store stress
        %------------------------------------------------------------------%
        
        %% PLASTIC CORRECTION: ONE STEP EXACT RETURN MAPPING
        beta         = abs(c.*cos(mpD.phi)-sig.*sin(mpD.phi))./tau        ;%
        dsigA        = 0.5*beta.*(ds)                                     ;%
        dsigB        = c./tan(mpD.phi)                                    ;%
        iDA          = (sig<=dsigB & f>0)                                 ;% logical indexing
        mpD.sn(1,iDA)= sig(iDA)+dsigA(iDA)                                ;%
        mpD.sn(2,iDA)= sig(iDA)-dsigA(iDA)                                ;%
        mpD.sn(4,iDA)= beta(iDA).*mpD.s(4,iDA)                            ;%
        iDB           = (sig> dsigB & f>0)                                ;% logical indexing
        mpD.sn(1,iDB) = dsigB(iDB)                                        ;%
        mpD.sn(2,iDB) = dsigB(iDB)                                        ;%
        mpD.sn(4,iDB) = 0.0                                               ;%
        %------------------------------------------------------------------%
        
        %% COMPUTE NEW STRESS, VOLUMETRIC, DEVIATORIC AND SECOND DEVIATORIC INVARIANT
        dsig         = mpD.sn-mpD.s                                       ;%
        mpD.s        = mpD.sn                                             ;%
        mpD.ep       = Del\dsig                                           ;%
        mpD.epV      = mpD.ep'*[1;1;1;0]                                  ;%
        mpD.devep    = mpD.ep-repmat(mpD.epV',4,1).*repmat([1;1;1;0],1,mpD.n);%
        mpD.epII     = mpD.epII+sqrt(2/3*(mpD.ep(1,:).^2+mpD.ep(2,:).^2+...
                                          mpD.ep(3,:).^2+2.*mpD.ep(4,:).^2))                         ;%
        %--------------------------------------------------------------------------%
        
    end
    %% SAVE
    D            = [mpD.epII coh]                                               ;%
    %--------------------------------------------------------------------------%
    
end

