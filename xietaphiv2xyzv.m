function [xv,yv,zv,x,y,z]=xietaphiv2xyzv(isProlate,conk,xiv,etav,phiv,xi,eta,phi);
% xietaphiv2xyzv - vector field transformation for spheroidal to cartesean.
%
% Usage:
%
% [xv,yv,zv,x,y,z]=xietaphiv2xyzv(isProlate,conk,xiv,etav,phiv,xi,eta,phi);
%
% conk - coordinate interfocal distance (real>0).
% xiv, etav, phiv - spheroidal vector field.
% xi, eta, phi - spheroidal coordinates.
% xv, yv, zv, - cartesean vector field.
% x, y, z, - cartesean coordinates.
% 
% [xyzv,xyz]=xietaphiv2xyzv(isProlate,conk,xiv,etav,phiv,xi,eta,phi);
%
% xyzv - [N x 3] cartesean vector field.
% xyz - [N x 3] cartesean coordinates.
%
% [xyzv,xyz]=xietaphiv2xyzv(isProlate,conk,xietaphiv,xietaphi);
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
    phi=etav(:,3);
    eta=etav(:,2);
    xi=etav(:,1);
    phiv=xiv(:,3);
    etav=xiv(:,2);
    xiv=xiv(:,1);
end
eta=eta(:);
xi=xi(:);
phi=phi(:);
etav=etav(:);
xiv=xiv(:);
phiv=phiv(:);

sigma=2*isProlate-1;

[x,y,z] = xietaphi2xyz(isProlate,conk,xi,eta,phi);


etafac=(1-eta.^2);
xifac=(-sigma+xi.^2);

a=xi.*sqrt((etafac)./(-sigma*eta.^2+xi.^2));
b=eta.*sqrt((xifac)./(-sigma*eta.^2+xi.^2));


%special cases on interval or disc planes:
a(abs(eta)==0)=1;
a(abs(eta)==1)=0;
b(abs(eta)==0)=0;
b(abs(eta)==1)=1; %maybe sqrt(sigma)

%Normalised jacobian of cartesian to prolate. Mathematica generated.
J=[ [a.*cos(phi),-b.*cos(phi),-sin(phi)];
    [a.*sin(phi),-b.*sin(phi),cos(phi)];
    [b,a,zeros(size(eta))]];

%Pack the spherical three vectors
rtpv=[xiv,etav,phiv];

%Separate the Jacobian and multiply for each unit vector.
xv = dot(J(1:length(xi),:),rtpv,2);
yv = dot(J(length(xi)+1:2*length(xi),:),rtpv,2);
zv = dot(J(2*length(xi)+1:3*length(xi),:),rtpv,2);

if nargout < 3
   xv = [ xv yv zv ];
   yv = [ x y z ];
end

return
