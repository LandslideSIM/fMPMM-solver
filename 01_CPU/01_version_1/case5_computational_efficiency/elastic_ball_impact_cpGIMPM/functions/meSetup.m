function [meD] = meSetup(nEx,L,typeD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% MESH INITIALIZATION
meD.vd  = 0.0                                                             ;% viscous damping coefficient
Lx      = L                                                               ;% mesh dimension (x,y)
Ly      = Lx;
nEy     = nEx                                                 ;% number of element in y direction
meD.L   = [Lx Ly]                                                    ;% mesh length in x            [m]
meD.h   = [meD.L(1)/nEx meD.L(2)/nEy]                                     ;% [dx dy]
[xn,yn] = meshgrid(-0.5*meD.L(1)-meD.h(1):meD.h(1):0.5*meD.L(1)+meD.h(1),...
                   (1.0*meD.L(2)+meD.h(2):-meD.h(2):0-meD.h(2)))                   ;%

meD.nNx = size(xn,2)                                                      ;% number of nodes along x
meD.nNy = size(yn,1)                                                      ;% number of nodes along y
meD.nN  = meD.nNx*meD.nNy                                                 ;% number of nodes
meD.nEx = meD.nNx-1;
meD.nEy = meD.nNy-1;
meD.nNe = 16                                                              ;% number of node per element
meD.DoF = 2                                                               ;% degree of freedom
meD.nDF = meD.DoF.*[meD.nNe meD.nN]                                       ;% local and global number of degree of freedom
meD.x   = xn(:)                                                           ;% x coordinate
meD.y   = yn(:)                                                           ;% y coordinate
% NODAL VECTOR INITIALIZATION                                              
meD.m   = zeros(meD.nN    ,1,typeD)                                       ;% mass
meD.mr  = zeros(meD.nDF(2),1,typeD)                                       ;% repmated mass
meD.f   = zeros(meD.nDF(2),1,typeD)                                       ;% force balance                                                       
meD.fi  = zeros(meD.nDF(2),1,typeD)                                       ;% internal force
meD.d   = zeros(meD.nDF(2),1,typeD)                                       ;% damping force
meD.a   = zeros(meD.nDF(2),1,typeD)                                       ;% acceleration
meD.p   = zeros(meD.nDF(2),1,typeD)                                       ;% momentum
meD.v   = zeros(meD.nDF(2),1,typeD)                                       ;% velocity
meD.u   = zeros(meD.nDF(2),1,typeD)                                       ;% displacement
%-------------------------------------------------------------------------%

%% ELEMENT-NODE CONNECTIVITY
[meD.e2N] = e2N(meD.nNy,meD.nNx,meD.nEx,meD.nEy,meD.nNe)                        ;% element to node topology
%-------------------------------------------------------------------------%

end

