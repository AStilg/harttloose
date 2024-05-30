function [x,y,z]=xietaphi2xyz(isProlate,conk,xi,eta,phi,xifac,etafac);
% xietaphi2xyz - compute the transformation of spheroidal coordinates back
% 				to cartesian. For accurate conversion near singular points 
% 				sqrt(xi.^2-sigma) and sqrt(1-eta.^2) should also be provided.
%
% usage:
% 
% [x,y,z] = xietaphi2xyz(isProlate,conk,xi,eta,phi,xifac,etafac); %(5)
%
% or 
% 
% [x,y,z] = xietaphi2xyz(isProlate,conk,xi,eta,phi,xifac); %(4)
% 
% or
% 
% [x,y,z] = xietaphi2xyz(isProlate,conk,xi,eta,phi); %(3)
%
% or
%
% [allparams] = xietaphi2xyz(isProlate,conk,xietaphi);
%
% NOTE: xietaphi can be [N x (3,4,5)] of the three (n) cases of the above.
% NOTE: allparams can always be used as an [N x 3] output.
%
% PACKAGE INFO

if nargin==3
    if size(xi,2)==3
        phi=xi(:,3);
        eta=xi(:,2);
        xi=xi(:,1);
    elseif size(xi,2)==4
        xifac=xi(:,4);
        phi=xi(:,3);
        eta=xi(:,2);
        xi=xi(:,1);
    elseif size(xi,2)==5
        etafac=xi(:,5);
        xifac=xi(:,4);
        phi=xi(:,3);
        eta=xi(:,2);
        xi=xi(:,1);
    end
end
if conk==0
    conk=1e-100;
end

if ~exist('xifac','var')
    xifac=sqrt(xi.^2-(2*isProlate-1));
end
if ~exist('etafac','var')
    etafac=sqrt((1-eta).*(1+eta));
end

x=conk.*etafac.*xifac;
y=x.*sin(phi);
x=x.*cos(phi);
z=conk.*eta.*xi;

if nargout<2
    x=[x(:),y(:),z(:)];
end