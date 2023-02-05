function [N,dNx,dNy] = f_SdS_linear(dX,dY,dx,dy)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% COMPUTE X-BASIS FUNCTION
[Sx,dSx] = f_S(dX,dx)                                                     ;%
%% COMPUTE Y-BASIS FUNCTION
[Sy,dSy] = f_S(dY,dy)                                                     ;%
%% CONVOLUTION OF BASIS FUNCTIONS
N    =  Sx.* Sy                                                           ;%  shape function
dNx  = dSx.* Sy                                                           ;% x-gradient shape function
dNy  =  Sx.*dSy                                                           ;% y-gradient shape function
end
function [S,dS]=f_S(dXi,hi,lp)
%% INITIALIZE BASIS FUNCTION MATRIX
S  = zeros(size(dXi));
dS = zeros(size(dXi));

%% COMPUTE BASIS FUNCTIONS
c1 = abs(dXi)<=hi                                                         ;%
% BASIS FUNCTION
S(c1==1)=(1-abs(dXi(c1==1))./hi)                                                 ;%
% BASIS FUNCTION GRADIENT
dS(c1==1)=(-sign(dXi(c1==1))./hi);
end
