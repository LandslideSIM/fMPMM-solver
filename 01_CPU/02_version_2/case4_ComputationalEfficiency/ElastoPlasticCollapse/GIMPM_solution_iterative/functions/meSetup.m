function [meD,bc] = meSetup(nEx,typeD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% MESH INITIALIZATION
meD.vd  = 0.10                                                            ;% viscous damping coefficient
Lx      = 0.80                                                         ;% mesh dimension (x,y)
Ly      = Lx                                                              ;% number of element in x direction
nEy     = nEx                                                             ;% number of element in y direction
meD.L   = [Lx Ly]                                                         ;% mesh length in x            [m]
meD.h   = [meD.L(1)/nEx meD.L(2)/nEy]                                     ;% [dx dy]

ratio   = 6.6667
Ly      = Lx/ratio                                                       ;% mesh dimension (x,y)
nEy     = nEx/ratio                                                              ;% number of element in y direction
meD.L   = [Lx Ly]                                                         ;% mesh length in x            [m]
meD.h   = [meD.L(1)/nEx meD.L(2)/nEy]                                     ;% [dx dy]
meD.L   = [Lx Ly]                                                         ;% mesh length in x            [m]

[xn,yn] = meshgrid(0.0-2*meD.h(1):meD.h(1):meD.L(1)+2*meD.h(1),...
                   0.0-2*meD.h(2):meD.h(2):meD.L(2)+2*meD.h(2))               ;%
xn      = flip(xn)                                                        ;%
yn      = flip(yn)                                                        ;%
               
               
meD.nNx = size(xn,2)                                                      ;% number of nodes along x
meD.nNy = size(yn,1)                                                      ;% number of nodes along y
meD.nN  = meD.nNx*meD.nNy                                                 ;% number of nodes
meD.nEx = meD.nNx-1;
meD.nEy = meD.nNy-1;
meD.nNe = 16                                                              ;% number of node per element
meD.DoF = 2                                                               ;% degree of freedom
meD.nDoF = meD.DoF.*[meD.nNe meD.nN]                                       ;% local and global number of degree of freedom
meD.x   = xn(:)                                                           ;% x coordinate
meD.y   = yn(:)                                                           ;% y coordinate
% NODAL VECTOR INITIALIZATION                                              
meD.m   = zeros(meD.nN    ,1,typeD)                                       ;% mass
meD.mr  = zeros(meD.nDoF(2),1,typeD)                                      ;% repmated mass
meD.f   = zeros(meD.nDoF(2),1,typeD)                                      ;% force balance                                                       
meD.fi  = zeros(meD.nDoF(2),1,typeD)                                      ;% internal force
meD.d   = zeros(meD.nDoF(2),1,typeD)                                      ;% damping force
meD.a   = zeros(meD.nDoF(2),1,typeD)                                      ;% acceleration
meD.p   = zeros(meD.nDoF(2),1,typeD)                                      ;% momentum
meD.v   = zeros(meD.nDoF(2),1,typeD)                                      ;% velocity
meD.u   = zeros(meD.nDoF(2),1,typeD)                                      ;% displacement
%--------------------------------------------------------------------------%

%% ELEMENT-NODE CONNECTIVITY
[meD.e2N] = e2N(meD.nNy,meD.nNx,meD.nEx,meD.nEy,meD.nNe)                  ;% element to node topology
%--------------------------------------------------------------------------%

% BOUNDARY CONDITIONS
meD.xB = [min(meD.x)+2*meD.h(1) max(meD.x)-2*meD.h(1) 0.0 max(meD.y)-2*meD.h(2)]                ;%
[row]  = find(meD.y<=meD.xB(3))                                           ;%
bc.y   = row                                                              ;%
[row]  = find(meD.x<=meD.xB(1))                                           ;%
bc.x   = row                                                              ;%
[row]  = find(meD.x>=meD.xB(2))                                           ;%
bc.x   = [bc.x;row]                                                       ;%
                                                
bc.y   = meD.DoF*bc.y(:,1)                                                ;%
bc.y   = [bc.y zeros(size(bc.y))]                                         ;%
bc.x   = meD.DoF*bc.x(:,1)-1                                              ;%  
bc.x   = [bc.x zeros(size(bc.x))]                                         ;%
end
function [g_num] = e2N(nny,nnx,nelx,nely,nnpe)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gnumbers= flip(reshape(1:(nnx*nny),nny ,nnx ))                            ;%
g_num   = zeros(nelx*nely,nnpe,'int32')                                   ;%
iel     = 1                                                               ;%
disp('------------------------')                                          ;%
disp('16 nodes per quadrilaterial element')                               ;%
disp('------------------------')                                          ;%
for i = 1:nelx
    for j = 1:nely
        if(i>1 && i<nelx && j>1 && j<nely)
            g_num(iel,1 ) = gnumbers(j-1,i-1)                         ;%
            g_num(iel,2 ) = gnumbers(j-0,i-1)                         ;%
            g_num(iel,3 ) = gnumbers(j+1,i-1)                         ;%
            g_num(iel,4 ) = gnumbers(j+2,i-1)                         ;%
            
            g_num(iel,5 ) = gnumbers(j-1,i  )                         ;%
            g_num(iel,6 ) = gnumbers(j-0,i  )                         ;%
            g_num(iel,7 ) = gnumbers(j+1,i  )                         ;%
            g_num(iel,8 ) = gnumbers(j+2,i  )                         ;%
            
            g_num(iel,9 ) = gnumbers(j-1,i+1)                         ;%
            g_num(iel,10) = gnumbers(j-0,i+1)                         ;%
            g_num(iel,11) = gnumbers(j+1,i+1)                         ;%
            g_num(iel,12) = gnumbers(j+2,i+1)                         ;%
            
            g_num(iel,13) = gnumbers(j-1,i+2)                         ;%
            g_num(iel,14) = gnumbers(j-0,i+2)                         ;%
            g_num(iel,15) = gnumbers(j+1,i+2)                         ;%
            g_num(iel,16) = gnumbers(j+2,i+2)                         ;%
        end
        iel = iel+1;
    end
end
end
