function [xiv,etav,phiv,xi,eta,phi]=rtpv2xietaphiv(isProlate,conk,rv,tv,pv,r,t,p);
%rtpv2xietaphiv - vector field transformation for spherical to spheroidal.
%
%Usage:
%
% [xiv,etav,phiv,xi,eta,phi]=rtpv2xietaphiv(isProlate,conk,rv,tv,pv,r,t,p);
%
% conk - coordinate interfocal distance (real>0).
% xiv, etav, phiv - spheroidal vector field.
% xi, eta, phi - spheroidal coordinates.
% rv, tv, pv, - spherical vector field.
% r, t, p, - spherical coordinates.
% 
% [xietaphiv,xietaphi]=rtpv2xietaphiv(isProlate,conk,rv,tv,pv,r,t,p);
%
% xietaphiv - [N x 3] spheroidal vector field.
% xietaphi - [N x 3] spheroidal coordinates.
%
% [xietaphiv,xietaphi]=rtpv2xietaphiv(isProlate,conk,rtpv,rtp);
%
% rtpv - [N x 3] spherical vector field.
% rtp - [N x 3] spherical coordinates.
%
% PACKAGE INFO

if nargin < 6
   r = tv(:,1);
   t = tv(:,2);
   p = tv(:,3);
   pv = rv(:,3);
   tv = rv(:,2);
   rv = rv(:,1);
end

sigma=2*isProlate-1;

% Convert points to spherical coordinates
[xi,eta,phi] = rtp2xietaphi(isProlate,conk,r,t,p);

eta=eta(:);
xi=xi(:);
phi=phi(:);

rv=rv(:);
tv=tv(:);
pv=pv(:);

etafac=sqrt(1-eta.^2);
xifac=sqrt(-sigma+xi.^2);

a=xi.*sqrt(-sigma*(etafac.^2)./(eta.^2-sigma*xi.^2));
b=eta.*sqrt((xifac.^2)./(-sigma*eta.^2+xi.^2));

%special cases on interval or disc planes:
% a(abs(eta)==0)=1;
% a(abs(eta)==1)=0;
% b(abs(eta)==0)=0;
% b(abs(eta)==1)=1; %maybe sqrt(sigma)

%Normalised jacobian of prolate to polar. Mathematica generated.
J=[ [a.*sin(t)+b.*cos(t),a.*cos(t)-b.*sin(t),0*phi];
    [a.*cos(t)-b.*sin(t),-a.*sin(t)-b.*cos(t),0*phi];
    [0*phi,0*phi,ones(size(eta))]];

%Pack the spherical three vectors
rtpv=[rv,tv,pv];

%Separate the Jacobian and multiply for each unit vector.
xiv = dot(J(1:length(xi),:),rtpv,2);
etav = dot(J(length(xi)+1:2*length(xi),:),rtpv,2);
phiv = dot(J(2*length(xi)+1:3*length(xi),:),rtpv,2);

if nargout < 3
   xiv = [ xiv etav phiv ];
   etav = [ xi eta phi ];
end

return
