function [meD] = p2Nsolve_loop(meD,mpD,g,dt,l2g,p2n,BC)
%UNTITLED Summary of this function goes here
% Input(s):
% meD - mesh structured array:
%           - nN : total number of ndoes
% mpD - material point structured array:
%           - n   : count
%           - V   : volume
%           - m   : mass
%           - p   : momentum
%           - N   : basis function
%           - B   : B matrix
%% INITIALIZATION
% NODAL VECTOR INITIALIZATION
meD.m(:) = 0.0                                                            ;% mass
meD.mr(:)= 0.0                                                            ;% repmated mass
meD.p(:) = 0.0                                                            ;% momentum
meD.f(:) = 0.0                                                            ;% force balance
meD.fi(:)= 0.0                                                            ;%
meD.d(:) = 0.0                                                            ;% damping force
meD.a(:) = 0.0                                                            ;% momentum
meD.v(:) = 0.0                                                            ;% momentum
for mp=1:mpD.n % for any material point
    for N=1:meD.nNe % for any neighboring node of any material pointv
        iDx        = N*meD.DoF-1                                          ;%
        iDy        = iDx+1                                                ;%
        
        l          = p2n(mp,N)                                            ;%
        meD.m(l)   = meD.m(l)+mpD.S(mp,N).*mpD.m(mp)                      ;%
        L          = meD.DoF*l                                            ;%
        meD.p(L-1) = meD.p(L-1)+mpD.S(mp,N).*mpD.p(mp,1)                  ;%
        meD.p(L  ) = meD.p(L  )+mpD.S(mp,N).*mpD.p(mp,2)                  ;%
        
        meD.f(L  ) = meD.f(L  )-mpD.S(mp,N).*mpD.m(mp).*g                 ;%
        meD.fi(L-1) = meD.fi(L-1)+mpD.V(mp).*(mpD.dSx(mp,N).*mpD.s(1,mp)+...
            mpD.dSy(mp,N).*mpD.s(4,mp))                                   ;%
        meD.fi(L  ) = meD.fi(L  )+mpD.V(mp).*(mpD.dSx(mp,N).*mpD.s(4,mp)+...
            mpD.dSy(mp,N).*mpD.s(2,mp))                                   ;%
    end
end
%% SOLVE EXPLICIT MOMENTUM BALANCE EQUATION
meD.f      = meD.f-meD.fi                                                 ;%
% UPDATE NODAL QUANTITITES
iDx        = 1:meD.DoF:meD.nDF(2)-1                                       ;%
iDy        = iDx+1                                                        ;%
% COMPUTE NODAL FORCE
meD.d(iDx) = sqrt(meD.f(iDx).^2+meD.f(iDy).^2)                            ;%
meD.d(iDy) = meD.d(iDx)                                                   ;%
meD.f      = meD.f - meD.vd.*meD.d.*sign(meD.p)                           ;%
% UPDATE NODAL MOMENTUM
meD.p      = meD.p + dt.*meD.f                                            ;%
% COMPUTE NODAL ACCELERATION AND VELOCITY
meD.mr     = reshape(repmat(meD.m',meD.DoF,1),meD.nDF(2),1)               ;%
meD.a      = meD.f./meD.mr                                                ;%
meD.v      = meD.p./meD.mr                                                ;%
meD.a(isnan(meD.a))=0.0                                                   ;%
meD.v(isnan(meD.v))=0.0                                                   ;%
% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
BCx =  [BC.xi;BC.xs]                                                      ;% boundary nodes index
BCy =  [BC.yi      ]                                                      ;% boundary nodes index
meD.a([meD.DoF*BCx-1;meD.DoF*BCy;meD.DoF*BCy-1])=0.0                      ;% fix BC's
meD.v([meD.DoF*BCx-1;meD.DoF*BCy;meD.DoF*BCy-1])=0.0                      ;% fix BC's
clear m p f fi iDx iDy iD BCx BCy                                         ;% clear temporary variables
%-------------------------------------------------------------------------%

end