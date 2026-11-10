% Used for: York regression used in Extended Data Fig. 8 and Supplementary Fig. 4
% Original: Manuscript\MLD\manuscript\nature\NCC\R1\code for revision\mld_npp_relationship_grid\york_reg.m (data paths updated to I:\Prj; computation unchanged)
function [a,b,sa,sb] = york_reg(X,Y,sX,sY)
% YORK_REG  York (1966) bivariate weighted-errors-in-both-variables regression.
%   [a,b,sa,sb] = york_reg(X,Y,sX,sY) fits Y = a + b*X given per-point
%   standard errors sX, sY on both variables. Iterative solution following
%   Cael et al. (2026) replication in I:\Prj3\Cael_2026\cael_replication_and_correction.m.
X=X(:); Y=Y(:); sX=sX(:); sY=sY(:);
wX=1./sX.^2; wY=1./sY.^2;
p=polyfit(X,Y,1); b=p(1);
for iter=1:200
    W=wX.*wY./(wX+b^2.*wY);
    Xb=sum(W.*X)/sum(W); Yb=sum(W.*Y)/sum(W);
    U=X-Xb; V=Y-Yb;
    beta=W.*(U./wY+b.*V./wX);
    den=sum(W.*beta.*U); if den==0, break; end
    bn=sum(W.*beta.*V)/den;
    if abs(bn-b)<1e-12, b=bn; break; end
    b=bn;
end
a=Yb-b*Xb; x=Xb+beta; xb=sum(W.*x)/sum(W); u=x-xb;
sb=sqrt(1/sum(W.*u.^2)); sa=sqrt(1/sum(W)+xb^2*sb^2);
end
