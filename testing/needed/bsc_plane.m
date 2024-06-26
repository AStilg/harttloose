function [nn,mm,a,b] = bsc_plane(nmax,lambda,theta,phi,Etheta,Ephi)
% Finds VSWF representation of plane wave, output BSCs producing a plane
% with amplitude determined by the field input. If a vector input then a
% and b are matrices.
%
% Usage:
% [n,m,a,b] = bsc_plane(nmax,lambda,theta,phi,Etheta,Ephi)
%
% or
%
% [a,b] = bsc_plane(nmax,lambda,theta,phi,Etheta,Ephi)
%
% theta and phi (in radians) give the direction of propagation
% of the plane wave. +z direction is theta = 0, phi = any value
%
% NOTE: lambda will be leaving this code in a future release.
%
% You must obtain/distribute a copy of the LICENSE with any derivations of this 
% file along with the following message.
%
% Author: Alexander Stilgoe (a.stilgoe@uq.edu.au)
% Copyright (C) The University of Queensland.
% This file is part of HARTTloose.
% The LICENSE can be obtained at: https://github.com/AStilg/harttloose/blob/main/LICENSE

%special case ... linear polarisation
if nargin==5
    szE=size(Etheta);
    if szE(2)==2
        if szE(1)==1
            Etheta=repmat(Etheta,[length(phi),1]);
        end
        Ephi=-sin(phi).*Etheta(:,1)+cos(phi).*Etheta(:,2);
        Etheta=sign(cos(real(theta))).*(cos(phi).*Etheta(:,1)+sin(phi).*Etheta(:,2));
    else
        error('ott:bsc_plane:badargs','Most likely: Invalid uniform polarisation. Otherwise: Invalid argument.');
    end
end

theta=theta(:);
phi=phi(:);
Etheta=1i*Etheta(:);
Ephi=1i*Ephi(:);

if any(theta>pi)
    warning('ott:bsc_plane:thetatoobig','Theta is larger than PI, treating all elements as angles in degrees.')
    theta=theta/180*pi;
end

%calculate the mode indices we are going to find.
[nn,mm]=combined_index([1:nmax*(nmax+2)].');

a = zeros(length(nn),length(Etheta));
b = zeros(length(nn),length(Etheta));

for n = 1:nmax
    
    iter=[(n-1)*(n+1)+1:n*(n+2)];
    leniter=2*n+1;
    
    %expand theta and phi components of field to match spherical harmonics
    ET=repmat(Etheta,[1,leniter]);
    EP=repmat(Ephi,[1,leniter]);
    
    %power normalisation.
    Nn = 1/sqrt(n*(n+1));
    
    %Generate the farfield components of the VSWFs
    [~,dtY,dpY] = spharm(n,[-n:n],theta,phi);
    
    %equivalent to dot((1i)^(n+1)*C,E);
    a(iter,:) = 4*pi*Nn*(-1i)^(n+1)*(conj(dpY).*ET - conj(dtY).*EP).';
    %equivalent to dot((1i)^(n)*B,E);
    b(iter,:) = 4*pi*Nn*(-1i)^(n)  *(conj(dtY).*ET + conj(dpY).*EP).';
    
end

p=abs(a).^2+abs(b).^2;
binaryvector=(p>1e-15*max(p));

if nargout>2
    nn=nn(any(binaryvector,2),:);
    mm=mm(any(binaryvector,2),:);
    a=a(any(binaryvector,2),:);
    b=b(any(binaryvector,2),:);
end

if nargout==2
    a(~binaryvector)=0;
    b(~binaryvector)=0;
    nn=sparse(a);
    mm=sparse(b);
end
return
