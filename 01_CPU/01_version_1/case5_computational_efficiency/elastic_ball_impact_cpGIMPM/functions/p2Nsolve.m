function [meD] = p2Nsolve(meD,mpD,g,dt,l2g,p2N,BC)
%f_P2NINT function that integrate every contributing material points to
%their neighboring nodes, which result in nodal mass, momentum and forces.
%The latter are further used to explicitly calculate a solution of the
%momentum balance equation.
%   Detailed explanation goes here
%% INITIALIZATION
% NODAL VECTOR INITIALIZATION
meD.m(:) = 0.0 ; meD.mr(:) = 0.0 ; meD.f(:) = 0.0 ; meD.d(:) = 0.0        ;% mass | repmated mass | force balance | damping force
meD.a(:) = 0.0 ;  meD.p(:) = 0.0 ; meD.v(:) = 0.0 ; meD.u(:) = 0.0        ;% acceleration | momentum | velocity | incremental displacement
%--------------------------------------------------------------------------%

%% CONTRIBUTION TO NODES
% PREPROCESSING
m = reshape( mpD.S.*repmat(mpD.m,1,meD.nNe)      ,mpD.n*meD.nNe   ,1)     ;% preprocesing of mass global vector
p = reshape([mpD.S.*repmat(mpD.p(:,1),1,meD.nNe);...
             mpD.S.*repmat(mpD.p(:,2),1,meD.nNe)],mpD.n*meD.nDF(1),1)     ;% preprocesing of momentum global vector
f = reshape([mpD.S.*0.0                         ;...
             mpD.S.*repmat(mpD.m,1,meD.nNe).*-g ],mpD.n*meD.nDF(1),1)     ;% preprocesing of external force global vector  
fi= squeeze(sum(permute(mpD.B  ,[2 1 3]).*...
                repmat((permute(mpD.s',[3 2 1])),meD.nDF(1),1),2).*...
                repmat( permute(mpD.V ,[3 2 1]) ,meD.nDF(1),1))           ;% preprocesing of internal force global vector  
% CONTRIBUTION FROM p TO N
meD.m = accumarray(p2N(:),m,[meD.nN     1])                               ;% mass global vector
meD.p = accumarray(l2g(:),p,[meD.nDF(2) 1])                               ;% momentum global vector
meD.f = accumarray(l2g(:),f,[meD.nDF(2) 1])                               ;% external force global vector
for n = 1:meD.nNe                                                          % BEGIN ITERATION OVER meD.nNe NEIGHBORING NODES FOR ALL MP                                                               
    l = [(meD.DoF*p2N(:,n)-1);(meD.DoF*p2N(:,n))]                         ;% local to global index
    meD.f = meD.f - accumarray(l,[fi(n*meD.DoF-1,:)';...
                                  fi(n*meD.DoF  ,:)'],[meD.nDF(2) 1])     ;% external force global vector - nodal internal force global vector
end                                                                        % END ITERATION OVER meD.nNe NEIGHBORING NODES FOR ALL MP
%--------------------------------------------------------------------------%

%% SOLVE EXPLICIT MOMENTUM BALANCE EQUATION
% UPDATE GLOBAL NODAL INFORMATIONS
iDx        = 1:meD.DoF:meD.nDF(2)-1                                       ;% x component of global vector node
iDy        = iDx+1                                                        ;% y component of global vector node
% COMPUTE GLOBAL NODAL FORCE                                                  
meD.d(iDx) = sqrt(meD.f(iDx).^2+meD.f(iDy).^2)                            ;% x component damping force
meD.d(iDy) = meD.d(iDx)                                                   ;% y component damping force
meD.f      = meD.f - meD.vd*meD.d.*sign(meD.p)                            ;% nodal force - damping force
% UPDATE GLOBAL NODAL MOMENTUM
meD.p      = meD.p + dt*meD.f                                             ;% forward euler momentum
% COMPUTE GLOBAL NODAL ACCELERATION AND VELOCITY
meD.mr     = reshape(repmat(meD.m',meD.DoF,1),meD.nDF(2),1)               ;% repmat nodal mass
iD         = meD.mr==0                                                    ;% null nodal mass index
meD.a      = meD.f./meD.mr                                                ;% compute nodal acceleration
meD.v      = meD.p./meD.mr                                                ;% compute nodal velocity
meD.a(iD)  = 0.0                                                          ;% zeroed nodal acceleration
meD.v(iD)  = 0.0                                                          ;% zeroed nodal velocity
% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
BCx =  [BC.xi;BC.xs]                                                      ;% boundary nodes index
BCy =  [BC.yi      ]                                                      ;% boundary nodes index
meD.a([meD.DoF*BCx-1;meD.DoF*BCy;meD.DoF*BCy-1])=0.0                      ;% fix BC's
meD.v([meD.DoF*BCx-1;meD.DoF*BCy;meD.DoF*BCy-1])=0.0                      ;% fix BC's
clear m p f fi iDx iDy iD BCx BCy                                         ;% clear temporary variables
%-------------------------------------------------------------------------%

end
