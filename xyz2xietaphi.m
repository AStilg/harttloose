function [xi,eta,phi,xifac,etafac]=xyz2xietaphi(isProlate,conk,x,y,z)
% xyz2xietaphi - computes the coordinate transformation from cartesian to
% 				spheroidal coordinates.
% 				For accurate values near singular points we must also compute 
% 				sqrt(xi.^2-sigma) and sqrt(1-eta.^2).
%
% usage:
%
% [xi,eta,phi]=xyz2xietaphi(isProlate,conk,x,y,z);
%
% or 
% 
% [xi,eta,phi,xifac,etafac]=xyz2xietaphi(isProlate,conk,x,y,z);
%
% or 
%
% [xietaphi,xifac_etafac]=xyz2xietaphi(isProlate,conk,x,y,z);
% 
% or
% 
% [allparameters]=xyz2xietaphi(isProlate,conk,xyz);
%
% NOTE: [Nx3]: xyz can always be be used.
%
% PACKAGE INFO

if nargin <5

   z = x(:,3);
   y = x(:,2);
   x = x(:,1);
   
end
if conk==0
    conk=1e-100;
end

conk=(conk);
phi=atan2((y),(x));
rho2=(x.^2+y.^2);
rho=sqrt(rho2);
r2=rho2+z.^2;

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


