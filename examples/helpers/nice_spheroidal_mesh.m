function [externalMesh,internalMesh]=nice_spheroidal_mesh(ac,nxi,neta,width,height);
% nice_spheroidal_mesh - design a pretty XZ mesh for graphics. Note: there
% can be two regions: an internal mesh and external mesh.
%
% [externalMesh,internalMesh]=nice_spheroidal_mesh(ac,nxi,neta,xwidth,zwidth);
%
% ac -- the [x,z] coordinate of the internal-to-external mesh interface.
% width -- the total width of the box. Includes origin.
% height -- the total height of the box. Includes origin.
% nxi -- the number of xi between the centre of the space and the outermost
%        needed range.
% neta -- the number of eta values around the disc/line.
%
% externalMesh -- structure of X,Y,Z for external mesh.
% internalMesh -- structure of X,Y,Z for internal mesh.
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

% ac=[0.5,1];
% maxr=max(ac);
% xrange=1.05*[-maxr,maxr];
% zrange=xrange;
% nxi=10;
% neta=21;

xrange=[-.5,.5]*width;
zrange=[-.5,.5]*height;

[isProlate,c]=aspect_ratio_to_conk(ac);

[XR,ZR]=meshgrid(xrange,zrange);

xep=xyz2xietaphi(isProlate,c,XR(:),0*XR(:),ZR(:)); %get outer box
xi_max=max(xep(:,1));

%create positive grid seeds:
[x_max,~,~]=xietaphi2xyz(isProlate,c,xi_max,0,0);

x_ext=linspace(ac(1),x_max,nxi).';
x_int=linspace(0,ac(1),nxi).';

[xi_ext,~,~]=xyz2xietaphi(isProlate,c,x_ext,0*x_ext,0*x_ext);
[xi_int,~,~]=xyz2xietaphi(isProlate,c,x_int,0*x_int,0*x_int);

%make eta nice on ac(1); we're lazy and are going to use a spline.
etad=(linspace(-1,1,1000)');
[xr,~,zr]=xietaphi2xyz(isProlate,c,xi_ext(1).*ones(size(etad)),etad,0*etad);
s=cumsum(sqrt(diff(xr).^2+diff(zr).^2));
s=[0;s]; % this is numerically accurate as trap.

eta=spline(s,etad(1:end),linspace(s(1),s(end),neta));
%see, we used a spline. I hate elliptic integrals!

%create grids duplicate for phi
[Xi_ext,Eta_ext]=meshgrid(xi_ext,eta);
Phi_ext=0*Xi_ext;
Xi_ext=[flipud(Xi_ext);Xi_ext];
Eta_ext=[flipud(Eta_ext);Eta_ext];
Phi_ext=[Phi_ext+pi;Phi_ext];

[X_ext,Y_ext,Z_ext]=xietaphi2xyz(isProlate,c,Xi_ext(:),Eta_ext(:),Phi_ext(:));
X_ext=reshape(X_ext,size(Xi_ext));
Y_ext=reshape(Y_ext,size(Xi_ext));
Z_ext=reshape(Z_ext,size(Xi_ext));

[Xi_int,Eta_int]=meshgrid(xi_int,eta);
Phi_int=0*Xi_int;
Xi_int=[flipud(Xi_int);Xi_int];
Eta_int=[flipud(Eta_int);Eta_int];
Phi_int=[Phi_int+pi;Phi_int];

[X_int,Y_int,Z_int]=xietaphi2xyz(isProlate,c,Xi_int(:),Eta_int(:),Phi_int(:));
X_int=reshape(X_int,size(Xi_int));
Y_int=reshape(Y_int,size(Xi_int));
Z_int=reshape(Z_int,size(Xi_int));
% deal with phi now

externalMesh.X=X_ext;
externalMesh.XI=Xi_ext;
externalMesh.Y=Y_ext;
externalMesh.ETA=Eta_ext;
externalMesh.Z=Z_ext;
externalMesh.PHI=Phi_ext;

internalMesh.X=X_int;
internalMesh.XI=Xi_int;
internalMesh.Y=Y_int;
internalMesh.ETA=Eta_int;
internalMesh.Z=Z_int;
internalMesh.PHI=Phi_int;

% surf(X_ext,Z_ext,Y_ext)
% hold on
% surf(X_int,Z_int,Y_int)
% hold off
% view(2)
% axis equal
% xlim([-1,1])
% ylim([-1,1])