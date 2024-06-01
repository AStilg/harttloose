function [xi,eta,phi,xifac,etafac]=rtp2xietaphi(isProlate,conk,r,theta,phi)
% xyz2xietaphi - computes the coordinate transformation from sph. polar to
% 					spheroidal coordinates.
%
% For accurate values near singular points we must also compute 
% sqrt(xi.^2-sigma) and sqrt(1-eta.^2).
%
% usage:
%
% [xi,eta,phi]=rtp2xietaphi(isProlate,conk,r,theta,phi);
%
% or 
% 
% [xi,eta,phi,xifac,etafac]=rtp2xietaphi(isProlate,conk,r,theta,phi);
%
% or 
%
% [xietaphi,xifac_etafac]=rtp2xietaphi(isProlate,conk,r,theta,phi);
% 
% or
% 
% [allparameters]=rtp2xietaphi(isProlate,conk,rtp);
%
% NOTE: [Nx3]: rtp can always be be used.
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

if nargin <5

   theta = r(:,2);
   phi = r(:,3);
   r = r(:,1);
   
end
if conk==0
    conk=1e-100;
end

phi=phi(:);

conk=(conk);
rho=r(:).*sin(theta(:));
z=r(:).*cos(theta(:));

if isProlate
    unifactor=(sqrt(rho.^2+(z-conk).^2)+sqrt(rho.^2+(z+conk).^2))/2;
    unifactor(unifactor<conk)=conk; %cannot be less than conk.
else
    transverseFactor=((sqrt((rho+conk).^2+z.^2)+sqrt((rho-conk).^2+z.^2))/2);
    transverseFactor(transverseFactor<conk)=conk; %cannot be less than conk.
    unifactor=sqrt((transverseFactor-conk).*(transverseFactor+conk));
end

% eta=2.*z./unifactor;
xi=unifactor./conk;
if isProlate
    xi(xi<1)=1;
end
eta=z./conk./xi;

if isProlate %some of the below are redundant
    eta(xi>conk)=z(xi>conk)./conk./xi(xi>conk);
    eta(rho==0&abs(z)>=conk)=sign(z(rho==0&abs(z)>=conk));
    xifac=sqrt((xi-1).*(xi+1));
    etafac=sqrt((1-eta).*(1+eta));
    etafac(abs(z)>conk)=rho(abs(z)>conk)./conk./xifac(abs(z)>conk);
    xifac(abs(z)<conk)=rho(abs(z)<conk)./conk./etafac(abs(z)<conk);
else
    eta(unifactor==0)=sqrt(conk.^2-rho(unifactor==0).^2)/conk;
    xi(unifactor==0)=z(unifactor==0)./conk./eta(unifactor==0);
    xi(z==0&rho<=conk)=0;
    xifac=sqrt(1+xi.^2);
    etafac=rho./conk./xifac;
end

if nargout==1
    xi=[xi,eta,phi,xifac,etafac];
end

if nargout==2
    xi=[xi,eta,phi];
    eta=[xifac,etafac];
end


