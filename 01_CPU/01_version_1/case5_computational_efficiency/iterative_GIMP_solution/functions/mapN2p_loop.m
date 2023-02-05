function [meD,mpD] = mapN2p_loop(meD,mpD,dt,l2g,p2n,BC)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for mp = 1:mpD.n
    iDx = meD.DoF*p2n(mp,:)-1                                             ;% x-component index
    iDy = iDx+1                                                           ;% y-component index
    mpD.v(mp,1) =   mpD.PF(1).*(mpD.v(mp,1)+dt*mpD.S(mp,:)*meD.a(iDx))+...
        mpD.PF(2).*(mpD.S(mp,:)*meD.v(iDx))                               ;% PIC-FLIP update
    mpD.v(mp,2) =   mpD.PF(1).*(mpD.v(mp,2)+dt*mpD.S(mp,:)*meD.a(iDy))+...
        mpD.PF(2).*(mpD.S(mp,:)*meD.v(iDy))                               ;% PIC-FLIP update
    % MP'S MOMENTUM UPDATE
    mpD.p(mp,:) = mpD.v(mp,:).*repmat(mpD.m(mp),1,meD.DoF)                ;% update mp momentum
    % MP'S COORDINATE UPDATE
    mpD.x(mp,:) = mpD.x(mp,:)+dt*[mpD.S(mp,:)*meD.v(iDx)...
                                  mpD.S(mp,:)*meD.v(iDy)]                 ;% update mp position
end
%% UPDATE NODAL MOMENTUM WITH UPDATED MP MOMENTUM
meD.p(:) = 0.0                                                            ;% initialized momentum
for mp=1:mpD.n % for any material point
    for N=1:meD.nNe % for any neighboring node of any material point
        L          = meD.DoF*p2n(mp,N)                                    ;%
        meD.p(L-1) = meD.p(L-1)+mpD.S(mp,N).*mpD.p(mp,1)                  ;%
        meD.p(L  ) = meD.p(L  )+mpD.S(mp,N).*mpD.p(mp,2)                  ;%
    end
end
meD.u = dt*meD.p./meD.mr                                                  ;% incremental displacement 
meD.u(isnan(meD.u))=0.0                                                   ;% -
%% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
BCx = [BC.xi        ;BC.xs      ]                                         ;% boundary nodes index
BCy = [BC.yi                    ]                                         ;% boundary nodes index
meD.u([meD.DoF*BCx-1;meD.DoF*BCy;meD.DoF*BCy-1])=0.0                      ;% fix BC's for displacement
clear iD BCx BCy                                                          ;% clear temporary variables
%--------------------------------------------------------------------------%

end

