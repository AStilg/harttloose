function [rv,thetav,phiv,r,t,p]=xietaphiv2rtpv(isProlate,conk,xiv,etav,phiv,xi,eta,p);
%xietaphiv2rtpv - vector field transformation for spheroidal to sph. polar.
%
%Usage:
%
% [rv,thetav,phiv,r,theta,phi]=xietaphiv2rtpv(isProlate,conk,xiv,etav,phiv,xi,eta,phi);
%
% conk - coordinate interfocal distance (real>0).
% xiv, etav, phiv - spheroidal vector field.
% xi, eta, phi - spheroidal coordinates.
% rv, thetav, phiv, - sph. polar vector field.
% x, y, z, - sph. polar coordinates.
% 
% [rtpv,rtp]=xietaphiv2rtpv(isProlate,conk,xiv,etav,phiv,xi,eta,phi);
%
% rtpv - [N x 3] sph. polar vector field.
% rtp - [N x 3] sph. polar coordinates.
%
% [rtpv,rtp]=xietaphiv2rtpv(isProlate,conk,xietaphiv,xietaphi);
%
% xietaphiv - [N x 3] spheroidal vector field.
% xietaphi - [N x 3] spheroidal coordinates.
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

if nargin<8
    p=etav(:,3);
    eta=etav(:,2);
    xi=etav(:,1);
    phiv=xiv(:,3);
    etav=xiv(:,2);
    xiv=xiv(:,1);
end
eta=eta(:);
xi=xi(:);
p=p(:);
etav=etav(:);
xiv=xiv(:);
phiv=phiv(:);

sigma=2*isProlate-1;

[r,t,p] = xietaphi2rtp(isProlate,conk,xi,eta,p);

etafac=sqrt(1-eta.^2);
xifac=sqrt(-sigma+xi.^2);

a=xi.*sqrt((etafac.^2)./(-sigma*eta.^2+xi.^2));
b=eta.*sqrt((xifac.^2)./(-sigma*eta.^2+xi.^2));

%special cases on interval or disc planes:
a(abs(eta)==0)=1;
a(abs(eta)==1)=0;
b(abs(eta)==0)=0;
b(abs(eta)==1)=1; %maybe sqrt(sigma)

%Normalised jacobian of prolate to polar. Mathematica generated.
J=[ [a.*sin(t)+b.*cos(t),a.*cos(t)-b.*sin(t),0*p];
    [a.*cos(t)-b.*sin(t),-a.*sin(t)-b.*cos(t),0*p];
    [0*p,0*p,ones(size(eta))]];

%Pack the spherical three vectors
rtpv=[xiv,etav,phiv];

%Separate the Jacobian and multiply for each unit vector.
rv = dot(J(1:length(xi),:),rtpv,2);
thetav = dot(J(length(xi)+1:2*length(xi),:),rtpv,2);
phiv = dot(J(2*length(xi)+1:3*length(xi),:),rtpv,2);

if nargout < 3
   rv = [ rv thetav phiv ];
   thetav = [ r t p ];
end

return
