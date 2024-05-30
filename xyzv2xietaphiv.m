function [xiv,etav,phiv,xi,eta,phi]=xyzv2xietaphiv(isProlate,conk,xv,yv,zv,x,y,z);
%xietaphiv2xyzv - vector field transformation for cartesean to spheroidal.
%
%Usage:
%
% [xiv,etav,phiv,xi,eta,phi]=xyzv2xietaphiv(isProlate,conk,xv,yv,zv,x,y,z);
%
% conk - coordinate interfocal distance (real>0).
% xiv, etav, phiv - spheroidal vector field.
% xi, eta, phi - spheroidal coordinates.
% xv, yv, zv, - cartesean vector field.
% x, y, z, - cartesean coordinates.
% 
% [xietaphiv,xietaphi]=xyzv2xietaphiv(isProlate,conk,xv,yv,zv,x,y,z);
%
% xietaphiv - [N x 3] spheroidal vector field.
% xietaphi - [N x 3] spheroidal coordinates.
%
% [xietaphiv,xietaphi]=xyzv2xietaphiv(isProlate,conk,xyzv,xyz);
%
% xyzv - [N x 3] cartesean vector field.
% xyz - [N x 3] cartesean coordinates.
%
% PACKAGE INFO

if nargin < 8
   z = yv(:,3);
   y = yv(:,2);
   x = yv(:,1);
   zv = xv(:,3);
   yv = xv(:,2);
   xv = xv(:,1);
end

sigma=2*isProlate-1;

% Convert points to spherical coordinates
[xi,eta,phi] = xyz2xietaphi(isProlate,conk,x,y,z);

eta=eta(:);
xi=xi(:);
phi=phi(:);

xv=xv(:);
yv=yv(:);
zv=zv(:);

etafac=(1-eta.^2);
xifac=(-sigma+xi.^2);

a=xi.*sqrt((etafac)./(-sigma*eta.^2+xi.^2));
b=eta.*sqrt((xifac)./(-sigma*eta.^2+xi.^2));

% fa=a./((a-b).*(a+b));
% fb=b./((a-b).*(a+b));

%special cases on interval or disc planes:
% a(abs(eta)==0)=1;
% a(abs(eta)==1)=0;
% b(abs(eta)==0)=0;
% b(abs(eta)==1)=1; %maybe sqrt(sigma)

%Normalised jacobian of cartesian to prolate. Mathematica generated.
% J=[ [-fa.*cos(phi),-fa.*sin(phi),-fb];
    % [ fb.*cos(phi), fb.*sin(phi),fa];
    % [-sin(phi),cos(phi),zeros(size(eta))]];

% Normalised jacobian of cartesian to prolate. Mathematica generated.
J=[ [a.*cos(phi),a.*sin(phi),b];
    [-b.*cos(phi),-b.*sin(phi),a];
    [-sin(phi),cos(phi),zeros(size(eta))]];

%Pack the spherical three vectors
rtpv=[xv,yv,zv];

%Separate the Jacobian and multiply for each unit vector.
xiv = dot(J(1:length(xi),:),rtpv,2);
etav = dot(J(length(xi)+1:2*length(xi),:),rtpv,2);
phiv = dot(J(2*length(xi)+1:3*length(xi),:),rtpv,2);

if nargout < 3
   xiv = [ xiv etav phiv ];
   etav = [ xi eta phi ];
end

return
