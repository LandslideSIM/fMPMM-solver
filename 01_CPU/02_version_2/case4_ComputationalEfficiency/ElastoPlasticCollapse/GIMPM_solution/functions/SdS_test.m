function [mpD] = SdS(meD,mpD,p2N)
% SdS : function which calculate the basis function between
% material points and nodes through the array p2N, which is a topological
% indexing between a material point and all the nodes it integrates to
%   Detailed explanation goes here
%% COMPUTE (X,Y)-BASIS FUNCTION
D        = (repmat(mpD.x(:,1),1,meD.nNe) - meD.x(p2N))                    ;% distance x from point to node
[Sx,dSx] = NdN(D,meD.h(1),repmat(mpD.l(:,1),1,meD.nNe))                   ;% - see function      
D        = (repmat(mpD.x(:,2),1,meD.nNe) - meD.y(p2N))                    ;% distance y from point to node
[Sy,dSy] = NdN(D,meD.h(2),repmat(mpD.l(:,2),1,meD.nNe))                   ;% - see function 
%-------------------------------------------------------------------------%

%% CONVOLUTION OF BASIS FUNCTIONS
mpD.S    =  Sx.* Sy                                                       ;% basis 
mpD.dSx  = dSx.* Sy                                                       ;% x-gradient basis 
mpD.dSy  =  Sx.*dSy                                                       ;% y-gradient basis  
%-------------------------------------------------------------------------%

% xc = accumarray(p2e,mpD.x(:,1),[meD.nNy*meD.nNx 1]);
% yc = accumarray(p2e,mpD.x(:,2),[meD.nNy*meD.nNx 1]);
% [el,~,ic] = unique(p2e);
% el_counts = accumarray(ic,1);
% xc(el)=xc(el)./el_counts;
% yc(el)=yc(el)./el_counts;
% 
% figure(2524624)
% plot(meD.x,meD.y,'s',xc,yc,'x')
% axis equal
% drawnow
% 
xc         = mean(meD.x(p2N),2);
yc         = mean(meD.y(p2N),2);
D          = (repmat(xc,1,meD.nNe) - meD.x(p2N))                    ;% distance x from point to node
[Sxc,dSxc] = NdN(D,meD.h(1),repmat(2*mpD.l(:,2),1,meD.nNe))                                              ;% - see function      
D          = (repmat(yc,1,meD.nNe) - meD.y(p2N))                    ;% distance y from point to node
[Syc,dSyc] = NdN(D,meD.h(2),repmat(2*mpD.l(:,2),1,meD.nNe))                                              ;% - see function 
dSxc       = dSxc.* Syc                                                       ;% x-gradient basis 
dSyc       =  Sxc.*dSyc                                                       ;% y-gradient basis
%-------------------------------------------------------------------------%

%% B MATRIX ASSEMBLY
iDx         = 1:meD.DoF:meD.nDF(1)-1                                   ;% x component global node index
iDy         = iDx+1                                                    ;% y component global node index
B1          = zeros(4,meD.nDF(1),mpD.n);
B1(1,iDx,:) = mpD.dSx'                                                 ;% -
B1(2,iDy,:) = mpD.dSy'                                                 ;% -
B1(4,iDx,:) = mpD.dSy'                                                 ;% -
B1(4,iDy,:) = mpD.dSx'                                                 ;% -
%-------------------------------------------------------------------------%

%% B MATRIX ASSEMBLY
iDx         = 1:meD.DoF:meD.nDF(1)-1                                   ;% x component global node index
iDy         = iDx+1                                                    ;% y component global node index
B2          = zeros(4,meD.nDF(1),mpD.n);
B2(1,iDx,:) = dSxc';
B2(1,iDy,:) = dSyc';
B2(2,iDx,:) = dSxc';
B2(2,iDy,:) = dSyc';
B2(4,iDx,:) = dSy'                ;
B2(4,iDy,:) = dSx'                ;

B3          = zeros(4,meD.nDF(1),mpD.n);
B3(1,iDx,:) = 0.5*(mpD.dSx +    dSxc)';
B3(1,iDy,:) = 0.5*(    dSyc-mpD.dSy )';
B3(2,iDx,:) = 0.5*(    dSxc-mpD.dSx )';
B3(2,iDy,:) = 0.5*(mpD.dSy +    dSyc)';
B3(4,iDx,:) = mpD.dSy'                ;
B*(4,iDy,:) = mpD.dSx'                ;


end
function [N,dN]=NdN(dX,h,lp)
%% COMPUTE BASIS FUNCTIONS
lp = 2*lp                                                                 ;% length of mp domain
c1 = ( abs(dX)< (  0.5*lp)                        )                       ;% logical array 1
c2 = ((abs(dX)>=(  0.5*lp)) & (abs(dX)<(h-0.5*lp)))                       ;% logical array 2
c3 = ((abs(dX)>=(h-0.5*lp)) & (abs(dX)<(h+0.5*lp)))                       ;% logical array 3
% BASIS FUNCTION
N1 = 1-((4*dX.^2+lp.^2)./(4*h.*lp))                                       ;% basis function according to c1
N2 = 1-(abs(dX)./h)                                                       ;% basis function according to c2
N3 = ((h+0.5*lp-abs(dX)).^2)./(2*h.*lp)                                   ;% basis function according to c3
N  = c1.*N1+c2.*N2+c3.*N3                                                 ;% basis function
% BASIS FUNCTION GRADIENT
dN1= -((8*dX)./(4*h.*lp))                                                 ;% gradient basis function according to c1
dN2= sign(dX).*(-1/h)                                                     ;% gradient basis function according to c2
dN3=-sign(dX).*(h+0.5*lp-abs(dX))./(h*lp)                                 ;% gradient basis function according to c3
dN = c1.*dN1+c2.*dN2+c3.*dN3                                              ;% gradient basis function
%-------------------------------------------------------------------------%

% % %% COMPUTE BASIS FUNCTIONS
% % c1 = (((-h-lp)<dX) & (dX<=(-h+lp)))                                       ;% logical array 1
% % c2 = (((-h+lp)<dX) & (dX<=(  -lp)))                                       ;% logical array 2
% % c3 = ((    -lp<dX) & (dX<=    lp) )                                       ;% logical array 3
% % c4 = ((     lp<dX) & (dX<=( h-lp)))                                       ;% logical array 4
% % c5 = ((( h-lp)<dX) & (dX<=( h+lp)))                                       ;% logical array 5
% % % BASIS FUNCTION
% % N1 = ((h+lp+dX).^2)./(4*h*lp)                                             ;% basis function according to c1
% % N2 = 1+(dX./h)                                                            ;% basis function according to c2
% % N3 = 1-((dX.^2+lp.^2)./(2*h*lp))                                          ;% basis function according to c3
% % N4 = 1-(dX./h)                                                            ;% basis function according to c4
% % N5 = ((h+lp-dX).^2)./(4*h*lp)                                             ;% basis function according to c5
% % N  = c1.*N1+c2.*N2+c3.*N3+c4.*N4+c5.*N5                                   ;% basis function
% % % BASIS FUNCTION GRADIENT
% % dN1 = (h+lp+dX)./(2*h*lp)                                                 ;% gradient basis function according to c1
% % dN2 = 1/h.*ones(size(dX))                                                 ;% gradient basis function according to c2
% % dN3 = -dX./(h*lp)                                                         ;% gradient basis function according to c3
% % dN4 = -1/h.*ones(size(dX))                                                ;% gradient basis function according to c4
% % dN5 = -(h+lp-dX)./(2*h*lp)                                                ;% gradient basis function according to c5
% % dN  = c1.*dN1+c2.*dN2+c3.*dN3+c4.*dN4+c5.*dN5                             ;% gradient basis function
% % %-------------------------------------------------------------------------%

end